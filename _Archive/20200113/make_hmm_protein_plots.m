% script to examine results of HMM inference
clear
close all
addpath('utilities')
% set folder paths
project = 'Dl-Ven_snaBAC-mCh_v3';
DropboxFolder =  'E:\Nick\LivemRNA\Dropbox (Personal)\';
[~, DataPath, FigRoot] =   header_function(DropboxFolder, project);

% specify inference parameters
K = 3;
w = 7;
load([DataPath 'nucleus_struct.mat'])
load([DataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'],'hmm_input_output')
resultsPath = [DataPath 'hmm_inference_mf\w' num2str(w) '_K' num2str(K) '\'];
% specify figure path
FigPath = [FigRoot '\' project '\hmm_figs\'];
mkdir(FigPath);
tic
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
%     protein_bin = inference_results(i).protein_bin;
    % iterate
    for j = 1:numel(particle_ids)
        ptID = particle_ids(j);
        tr_ft = find(particle_index == ptID);
        if isempty(tr_ft)
            continue
        end
        mf_pt = mf_protein_index(tr_ft(1));
        fluo = hmm_input_output(tr_ft(1)).fluo;
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
        options = optimoptions(@lsqnonlin,'Display','off');
        r_fit = lsqnonlin(r_ob_fun,[1],[0],[Inf],options);
        r_fit = r_fit / dT;
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
        trace_inf_struct(iter).r = r_fit;
%         trace_inf_struct(iter).mf_id = protein_bin;        
        % increment
        iter = iter + 1;        
    end         
end

% plot results
close all
n_boots = 100;
mf_protein_index = [trace_inf_struct.mf_protein];
sample_vec = 1:numel(mf_protein_index)/5;
mf_axis_vec = linspace(prctile(mf_protein_index,1),prctile(mf_protein_index,99),20);
mf_sigma = .05;
burst_amp_full = [trace_inf_struct.r];
burst_dur_full = [trace_inf_struct.burst_duration]/60;
burst_freq_full = 60./[trace_inf_struct.burst_separation];
% initialize mean trend arrays
burst_dur_array = NaN(n_boots,numel(mf_axis_vec));
burst_freq_array = NaN(n_boots,numel(mf_axis_vec));
burst_amp_array = NaN(n_boots,numel(mf_axis_vec));

for n = 1:n_boots
    boot_indices = randsample(sample_vec,numel(sample_vec),true);
    burst_dur_boot = burst_dur_full(boot_indices);
    burst_freq_boot = burst_freq_full(boot_indices);
    burst_amp_boot = burst_amp_full(boot_indices);
    mf_protein_boot = mf_protein_index(boot_indices);
    for i = 1:numel(mf_axis_vec)
        dm_vec = exp(-(mf_axis_vec(i)-mf_protein_boot).^2 / 2 / mf_sigma^2);
        burst_dur_array(n,i) = nansum(burst_dur_boot.*dm_vec) / nansum(dm_vec);
        burst_freq_array(n,i) = nansum(burst_freq_boot.*dm_vec) / nansum(dm_vec);
        burst_amp_array(n,i) = nansum(burst_amp_boot.*dm_vec) / nansum(dm_vec);
    end
end
% calculate bootstrap average and ste
burst_dur_mean = nanmean(burst_dur_array);
burst_dur_ste = nanstd(burst_dur_array);

burst_freq_mean = nanmean(burst_freq_array);
burst_freq_ste = nanstd(burst_freq_array);

burst_amp_mean = nanmean(burst_amp_array);
burst_amp_ste = nanstd(burst_amp_array);
toc

% calculate Voxel Size
PixelSize = nucleus_struct(1).PixelSize;
zStep = nucleus_struct(1).zStep;
VoxelSize = PixelSize^2 * zStep;
mf_axis_vec = mf_axis_vec * VoxelSize;

MarkerSize = 75;
blue = [115 143 193]/256;
purple = [171 133 172]/256;
red = [213 108 85]/256;
close all
r_trend = figure;
hm_cm = flipud(brewermap([],'Spectral'));
colormap(hm_cm);
hold on
e = errorbar(mf_axis_vec,burst_amp_mean,burst_amp_ste,'Color','black','LineWidth',1.5);
e.CapSize = 0;
scatter(mf_axis_vec,burst_amp_mean,MarkerSize,'o','MarkerFaceColor',red,'MarkerEdgeAlpha',0);
p = plot(0,0);
grid on
xlim([0.54 1.8])
ylim([0 1.7])
xlabel('Dl concentration (au)')
ylabel('burst amplitude (au)')
% set(gca,'Fontsize',14)
StandardFigure(p,gca)
box on
saveas(r_trend,[FigPath,'burst_amp_mf_pt.tif'])
saveas(r_trend,[FigPath,'burst_amp_mf_pt.pdf'])

dur_trend = figure;
hold on
e = errorbar(mf_axis_vec,burst_dur_mean,burst_dur_ste,'Color','black','LineWidth',1.5);
e.CapSize = 0;
scatter(mf_axis_vec,burst_dur_mean,MarkerSize,'o','MarkerFaceColor',blue,'MarkerEdgeAlpha',0);
grid on
xlim([100 300])
p = plot(0,0);
grid on
xlim([0.54 1.8])
ylim([0 2])
xlabel('Dl concentration (au)')
ylabel('burst duration (min)')
% set(gca,'Fontsize',14)
StandardFigure(p,gca)
box on
saveas(dur_trend,[FigPath,'burst_dur_mf_pt.tif'])
saveas(dur_trend,[FigPath,'burst_dur_mf_pt.pdf'])


freq_trend = figure;
hold on
e = errorbar(mf_axis_vec,1./burst_freq_mean,burst_freq_ste./(burst_freq_mean.^2),'Color','black','LineWidth',1.5);
e.CapSize = 0;
scatter(mf_axis_vec,1./burst_freq_mean,MarkerSize,'o','MarkerFaceColor',purple,'MarkerEdgeAlpha',0);
grid on
xlim([0.54 1.8])
ylim([0 1.2])
xlabel('Dl concentration (au)')
ylabel('burst separation (min)')
p = plot(0,0);
StandardFigure(p,gca)
box on
saveas(freq_trend,[FigPath,'burst_sep_mf_pt.tif'])
saveas(freq_trend,[FigPath,'burst_sep_mf_pt.pdf'])