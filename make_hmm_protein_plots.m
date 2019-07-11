% script to examine results of HMM inference
clear
close all
addpath('utilities')
% set folder paths
project = 'Dl-Ven_snaBAC-mCh';
dropboxFolder =  'E:\Nick\LivemRNA\Dropbox\';
dataPath = [dropboxFolder 'ProcessedEnrichmentData\' project '/'];
% specify inference parameters
K = 3;
w = 7;
load([dataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'],'hmm_input_output')
resultsPath = [dataPath 'hmm_inference_protein\w' num2str(w) '_K' num2str(K) '\'];
% specify figure path
figPath = [dropboxFolder '\LocalEnrichmentFigures\' project '\hmm_protein_inference_w' num2str(w) '_K' num2str(K) '\'];
mkdir(figPath);
% read in inference results
inf_files = dir([resultsPath '*.mat']);
inference_results = struct;
iter = 1;
for i = 1:numel(inf_files)
    load([resultsPath inf_files(i).name])
    if ~output.skip_flag
        fn_list = fieldnames(output);
        for f = 1:numel(fn_list)
            inference_results(iter).(fn_list{f}) = output.(fn_list{f});
        end
        iter = iter + 1;
    end
end
% make indexing arrays
particle_index = [hmm_input_output.ParticleID];
mf_protein_index = NaN(size(particle_index));
for i = 1:numel(hmm_input_output)
    mf_protein = hmm_input_output(i).mf_protein;
    mf_protein_index(i) = nanmean(mf_protein);
end
dT = inference_results(1).deltaT;
% make alpha kernel
alpha_kernel = zeros(1,w);
alpha = inference_results(1).alpha;
for i = 1:w
    if i-1 > alpha
        alpha_kernel(i) = 1;
    elseif i > alpha
        alpha_kernel(i) = i-alpha + (i-1)/alpha * (alpha-i+1) + (1-(i-1)/alpha)*.5*(alpha-i+1);
    elseif i == 1
        alpha_kernel(i) = .5*i/alpha;
    else
        error('uh oh')
    end
end
        
% Calculate burst freq, burst dur, and r for each trace
trace_inf_struct = struct;
iter = 1;
for i  = 1:numel(inference_results)
    particle_ids = inference_results(i).particle_ids;
    p_z_cell = inference_results(i).soft_struct.p_z_log_soft;
    p_zz_cell = inference_results(i).soft_struct.p_zz_log_soft;  
    r_avg = inference_results(i).r;
    [r_sorted, ri] = sort(r_avg);
    protein_bin = inference_results(i).protein_bin;
    % iterate
    for j = 1:numel(particle_ids)
        ptID = particle_ids(j);
        tr_ft = particle_index == ptID;
        if ~any(tr_ft)
            continue
        end
        mf_pt = mf_protein_index(tr_ft);
        fluo = hmm_input_output(tr_ft).fluo;
        pz_array = exp(p_z_cell{j}(ri,:));
        pzz_mat = sum(exp(p_zz_cell{j}(ri,ri,:)),3);
        
        % solve for r
        pz2 = vertcat(pz_array(1,:), sum(pz_array(2:3,:)));
        pz1_counts = conv(pz2(1,:),alpha_kernel);
        pz1_counts = pz1_counts(1:numel(fluo));
        pz2_counts = conv(pz2(2,:),alpha_kernel);
        pz2_counts = pz2_counts(1:numel(fluo));
        % define objective fun
        r_ob_fun = @(params)  pz2_counts*params(1) - fluo;
        r_fit = lsqnonlin(r_ob_fun,[1],[0],[Inf]);
        
        % solve for 2 state burst dur and burst freq
        pzz2 = NaN(2,2);
        pzz2(1,1) = pzz_mat(1,1);
        pzz2(1,2) = sum(pzz_mat(1,2:3));
        pzz2(2,1) = sum(pzz_mat(2:3,1));
        pzz2(2,2) = sum(sum(pzz_mat(2:3,2:3)));
        a_mat = pzz2 ./ sum(pzz2);
        r_mat = logm(a_mat) / dT;
        % solve
        burst_dur = -1/r_mat(2,2);
        burst_sep = -1/r_mat(1,1);
        
        % record
        trace_inf_struct(iter).mf_protein = mf_pt;
        trace_inf_struct(iter).ParticleID = ptID;
        trace_inf_struct(iter).burst_duration = burst_dur;
        trace_inf_struct(iter).burst_separation = burst_sep;
        trace_inf_struct(iter).burst_separation = burst_sep;
        trace_inf_struct(iter).r = r_fit;
        trace_inf_struct(iter).mf_id = protein_bin;        
        % increment
        iter = iter + 1;        
    end         
end

%% plot results
close all

mf_protein_index = [trace_inf_struct.mf_protein];
mf_axis_vec = linspace(prctile(mf_protein_index,1),prctile(mf_protein_index,99));
mf_sigma = 7.5;
burst_amp_full = [trace_inf_struct.r];
burst_dur_full = [trace_inf_struct.burst_duration];
burst_sep_full = [trace_inf_struct.burst_separation];
% initialize mean trend arrays
burst_dur_vec = NaN(size(mf_axis_vec));
burst_sep_vec = NaN(size(mf_axis_vec));
burst_amp_vec = NaN(size(mf_axis_vec));

for i = 1:numel(mf_axis_vec)
    dm_vec = exp(-(mf_axis_vec(i)-mf_protein_index).^2 / 2 / mf_sigma^2);
    burst_dur_vec(i) = nansum(burst_dur_full.*dm_vec) / nansum(dm_vec);
    burst_sep_vec(i) = nansum(burst_sep_full.*dm_vec) / nansum(dm_vec);
    burst_amp_vec(i) = nansum(burst_amp_full.*dm_vec) / nansum(dm_vec);
end

r_scatter = figure;
hm_cm = flipud(brewermap([],'Spectral'));
colormap(hm_cm);
hold on
scatter([trace_inf_struct.mf_protein],[trace_inf_struct.r],30,'o','MarkerFaceColor',hm_cm(end-10,:),'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.3);
plot(mf_axis_vec,burst_amp_vec,'Color',hm_cm(end-10,:)/1.2,'LineWidth',2)
grid on
xlim([100 300])
ylim([15 40])
xlabel('average Dl concentration (au)')
ylabel('transcription burst amplitude (au)')
set(gca,'Fontsize',14)
saveas(r_scatter,[figPath,'burst_amp_mf_pt.tif'])


dur_scatter = figure;
hm_cm = flipud(brewermap([],'Spectral'));
colormap(hm_cm);
hold on
scatter(mf_protein_index,burst_dur_full,30,'o','MarkerFaceColor',hm_cm(20,:),'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.3);
plot(mf_axis_vec,burst_dur_vec,'Color',hm_cm(20,:)/1.2,'LineWidth',2)
grid on
xlim([100 300])
ylim([0 150])
xlabel('average Dl concentration (au)')
ylabel('transcription burst duration (au)')
set(gca,'Fontsize',14)
saveas(dur_scatter,[figPath,'burst_dur_mf_pt.tif'])

sep_scatter = figure;
hm_cm = flipud(brewermap([],'Spectral'));
colormap(hm_cm);
hold on
scatter(mf_protein_index,burst_sep_full,30,'o','MarkerFaceColor',hm_cm(5,:),'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.3);
plot(mf_axis_vec,burst_sep_vec,'Color',hm_cm(5,:)/1.2,'LineWidth',2)
grid on
xlim([100 300])
ylim([0 175])
xlabel('average Dl concentration (au)')
ylabel('transcription separation duration (au)')
set(gca,'Fontsize',14)
saveas(sep_scatter,[figPath,'burst_sep_mf_pt.tif'])