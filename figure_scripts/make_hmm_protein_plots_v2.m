% script to examine results of HMM inference
clear
close all
addpath('../utilities')
% set folder paths
project = '2xDl-Ven_snaBAC-mCh_v4';
DropboxFolder = 'S:\Nick\Dropbox\';
[~, DataPath, FigRoot] =   header_function(DropboxFolder, project);

% specify inference parameters
K = 3;
w = 7;
fluo = 2;
load([DataPath 'nucleus_struct_protein.mat'])

resultsPath = [DataPath 'hmm_inference_protein\w' num2str(w) '_K' num2str(K)  '_f' num2str(fluo) 'D\'];
% specify figure path
FigPath = [FigRoot '\' project '\hmm_figs\w' num2str(w) '_K' num2str(K)  '_f' num2str(fluo) 'D\'];
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
Tres = inference_results(1).deltaT;
% get list of mf dorsal values used for inference binsq
mf_axis_vec = inference_results(1).protein_bin_edges(1:end-1) + diff(inference_results(1).protein_bin_edges)/2;
dorsal_bins = 1:numel(inference_results(1).protein_bin_list);
bin_id_vec = [inference_results.protein_bin];

% check that bin assignmnets were performed correctly
pt_id_vec = [];
mf_id_vec = [];
for i = 1:numel(inference_results)
    pt_id_vec = [pt_id_vec inference_results(i).particle_ids];
    mf_id_vec = [mf_id_vec repelem(inference_results(i).protein_bin,numel(inference_results(i).particle_ids))];
end

pt_mf_key = unique([pt_id_vec' mf_id_vec'],'rows');
hmm_pt_vec = [nucleus_struct_protein.ParticleID];
hmm_pt_vec = hmm_pt_vec(~isnan(hmm_pt_vec));

protein_check_vec = NaN(1,size(pt_mf_key,1));
for i = 1:size(pt_mf_key,1)
    indices = find(hmm_pt_vec==pt_mf_key(i,1));
    if ~isempty(indices)
        protein_check_vec(i) = nanmean(nucleus_struct_protein(indices(1)).mf_null_protein_vec);
    end
end

pt_check_fig = figure;
scatter(pt_mf_key(:,2),protein_check_vec,5,'filled')
xlabel('Dorsal bin')
ylabel('average Dorsal')
grid on
set(gca,'Fontsize',14)
saveas(pt_check_fig,[FigPath 'protein_binning_check.png'])

%% calculate average initiation rate, burst freq, and burst duration
init_vec_mean = NaN(size(dorsal_bins));
init_vec_ste = NaN(size(dorsal_bins));
freq_vec_mean = NaN(size(dorsal_bins));
freq_vec_ste = NaN(size(dorsal_bins));
dur_vec_mean = NaN(size(dorsal_bins));
dur_vec_ste = NaN(size(dorsal_bins));
fluo_mean = NaN(size(dorsal_bins));
for d = dorsal_bins
    d_ids = find(bin_id_vec==d);
    init_array = NaN(size(d_ids));
    freq_array = NaN(size(d_ids));
    dur_array = NaN(size(d_ids));
    fluo_array = NaN(size(d_ids));
    for i = 1:numel(d_ids)
        [r,ri] = sort(inference_results(d_ids(i)).r);
        A = inference_results(d_ids(i)).A_mat(ri,ri);
        
        R = logm(A) / Tres;
        if ~isreal(R) || sum(R(:)<0) > K
            out = prob_to_rate_fit_sym(A, Tres, 'gen', .005, 1);            
            R = out.R_out;     
        end
        [V,D] = eig(R);
        [~,di] = max(diag(D));
        ss_vec = V(:,di) / sum(V(:,di));
        % initiation rate
        init_array(i) = (r(2) * ss_vec(2) + r(3) * ss_vec(3)) / (ss_vec(2)+ss_vec(3));   
        if sum(r.*ss_vec) > 1.1*init_array(i)
            error('nontrivial "off" state initiation')
        end            
        % burst freq
        freq_array(i) = -R(1,1);%_eff(2,1);
        % burst dur
        dur_array(i) = -1/R(1,1) *(1/ss_vec(1) - 1);%1/R_eff(1,2);
        
        fluo_array(i) = nanmean([inference_results(d_ids(i)).fluo_data{:}]);
    end
    % calculate average and ste
    init_vec_mean(d) = nanmean(init_array)*60;
    init_vec_ste(d) = nanstd(init_array)*60;
    
    freq_vec_mean(d) = nanmean(freq_array)*60;
    freq_vec_ste(d) = nanstd(freq_array)*60;
    
    dur_vec_mean(d) = nanmean(dur_array)/60;
    dur_vec_ste(d) = nanstd(dur_array)/60;
    
    fluo_mean(d) = nanmean(fluo_array);
end              
%%
MarkerSize = 50;
blue = [115 143 193]/256;
purple = [171 133 172]/256;
red = [213 108 85]/256;
ind_list = [2:numel(mf_axis_vec)-1];
mf_axis_long = linspace(min(mf_axis_vec(ind_list)),max(mf_axis_vec(ind_list)));
% set x axis
x_lim = [min(mf_axis_long)-.05 max(mf_axis_long)+.05];
% x_lim = [.85 2];
% mf_axis_long = linspace(min(mf_axis_vec(ind_list)),max(mf_axis_vec(ind_list)));

p_init = polyfit(mf_axis_vec(ind_list),init_vec_mean(ind_list),2);
p_trend_init = polyval(p_init,mf_axis_long);


mf_sigma = 5*median(diff(mf_axis_long));
init_trend_sm = NaN(size(mf_axis_long));
% calculate smoothed trend
for i = 1:numel(mf_axis_long)
    mf = mf_axis_long(i);
    sigma_wt_vec = exp(-.5*((mf_axis_vec(ind_list)-mf)/mf_sigma).^2);
    init_trend_sm(i) = nansum(init_vec_mean(ind_list).*sigma_wt_vec) / nansum(sigma_wt_vec);
end
    
close all
r_trend = figure;
hm_cm = flipud(brewermap([],'Spectral'));
colormap(hm_cm);
hold on
p1 = plot(mf_axis_long,p_trend_init,'--','Color','black','LineWidth',1.5);
e = errorbar(mf_axis_vec(ind_list),init_vec_mean(ind_list),init_vec_ste(ind_list),'o','Color','black','LineWidth',1);
e.CapSize = 0;
s = scatter(mf_axis_vec(ind_list),init_vec_mean(ind_list),MarkerSize,'o','MarkerFaceColor',red,'MarkerEdgeColor','black');
p = plot(0,0);
% grid on
xlim(x_lim)
% ylim([50 95])
xlabel('Dorsal concentration (au)')
ylabel('burst amplitude (au/min)')
% legend([s p1],'raw HMM results','trend','Location','southeast')
% set(gca,'Fontsize',14)
StandardFigure(p,gca)
box on
saveas(r_trend,[FigPath,'burst_amp_mf_pt.tif'])
saveas(r_trend,[FigPath,'burst_amp_mf_pt.pdf'])

%

p_dur = polyfit(mf_axis_vec(ind_list),dur_vec_mean(ind_list),1);
p_trend_dur = polyval(p_dur,mf_axis_long);

dur_trend = figure;
hold on
p1 = plot(mf_axis_long,p_trend_dur,'--','Color','black','LineWidth',1.5);
e = errorbar(mf_axis_vec(ind_list),dur_vec_mean(ind_list),dur_vec_ste(ind_list),'o','Color','black','LineWidth',1);
e.CapSize = 0;
s = scatter(mf_axis_vec(ind_list),dur_vec_mean(ind_list),MarkerSize,'o','MarkerFaceColor',blue,'MarkerEdgeColor','black');
% grid on
xlim(x_lim)
p = plot(0,0);
% grid on
xlim(x_lim)
% ylim([.5 3])
xlabel('Dorsal concentration (au)')
ylabel('burst duration (min)')
% legend([s p1],'raw HMM results','trend','Location','northeast')
% set(gca,'Fontsize',14)
StandardFigure(p,gca)
box on
saveas(dur_trend,[FigPath,'burst_dur_mf_pt.tif'])
saveas(dur_trend,[FigPath,'burst_dur_mf_pt.pdf'])

%
p_freq = polyfit(mf_axis_vec(ind_list),freq_vec_mean(ind_list),1);
p_trend_freq = polyval(p_freq,mf_axis_long);

freq_trend = figure;
hold on
p1 = plot(mf_axis_long,p_trend_freq,'--','Color','black','LineWidth',1.5);
e = errorbar(mf_axis_vec(ind_list),freq_vec_mean(ind_list),freq_vec_ste(ind_list),'o','Color','black','LineWidth',1);
e.CapSize = 0;
s = scatter(mf_axis_vec(ind_list),freq_vec_mean(ind_list),MarkerSize,'o','MarkerFaceColor',purple,'MarkerEdgeColor','black');
% grid on
xlim(x_lim)
% ylim([.5 1.7])
xlabel('Dorsal concentration (au)')
ylabel('burst frequency (1/min)')
p = plot(0,0);
% legend([s p1],'raw HMM results','trend','Location','southeast')

% legend([s p1],'raw HMM results','trend','Location','northeast')
StandardFigure(p,gca)
box on
saveas(freq_trend,[FigPath,'burst_sep_mf_pt.tif'])
saveas(freq_trend,[FigPath,'burst_sep_mf_pt.pdf'])