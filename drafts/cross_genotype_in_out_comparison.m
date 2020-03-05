% Script to generate figures establishing presence of enrichment bursts at
% start of transcription bursts
clear 
% close all
addpath('utilities')
% set ID variables
project1 = 'Dl-Ven_snaBAC-mCh_F-F-F_v1';
project2 = 'Dl-Ven_snaBAC-mCh_v3';
DropboxFolder = 'S:\Nick\Dropbox\';

% Params
fluo_dim = 3;
K = 3;
w = 7;

% set write paths
[~, DataPath1, FigureRoot] =   header_function(DropboxFolder, project1); 
[~, DataPath2, ~] =   header_function(DropboxFolder, project2); 

FigPath = [FigureRoot '\' project1 '\input_comparisons\'];
mkdir(FigPath)

% load data
load([DataPath1 'hmm_input_output_results_w' num2str(w) '_K' num2str(K) '_f' num2str(fluo_dim) '.mat'])
results_struct1 = results_struct;
load([DataPath2 'hmm_input_output_results_w' num2str(w) '_K' num2str(K) '_f' num2str(fluo_dim) '.mat'])
results_struct2 = results_struct;
clear results_struct;

Tres = 20; % seconds
% extract relevant arrays from project 1
lag_dur_vec_1 = results_struct1.lag_dur_vec;
lead_dur_vec_1 = results_struct1.lead_dur_vec;
hmm_array_1 = results_struct1.hmm_array;
spot_array_dt_1 = results_struct1.spot_array_dt;
spot_array_dm_1 = results_struct1.spot_array_dm;
virtual_array_dt_1 = results_struct1.virtual_array_dt;
virtual_array_dm_1 = results_struct1.virtual_array_dm;
feature_sign_vec_1 = results_struct1.feature_sign_vec;
% now project 2
lag_dur_vec_2 = results_struct2.lag_dur_vec;
lead_dur_vec_2 = results_struct2.lead_dur_vec;
hmm_array_2 = results_struct2.hmm_array;
spot_array_dt_2 = results_struct2.spot_array_dt;
spot_array_dm_2 = results_struct2.spot_array_dm;
virtual_array_dt_2 = results_struct2.virtual_array_dt;
virtual_array_dm_2 = results_struct2.virtual_array_dm;
feature_sign_vec_2 = results_struct2.feature_sign_vec;

%  determine snip size
n_col = size(spot_array_dt_1,2);
window_size = floor(n_col/2);
time_axis = (-window_size:window_size)*Tres/60;

%% set basic analyisis parameters
nBoots = 100; % number of bootstrap samples to use
% min_pause_len = 5; % minimum length of preceding OFF period (in time steps)
% max_pause_len = 10;
min_pause_len = 6; % minimum length of preceding OFF period (in time steps)
max_pause_len = 1000;
min_burst_len = 3;
max_burst_len = 1000;
% max_burst_len = 12;

% generate basic filter for target locus and computational controls
burst_ft_1 = feature_sign_vec_1 == 1&lead_dur_vec_1>=min_pause_len&lead_dur_vec_1<=max_pause_len...
    &lag_dur_vec_1>=min_burst_len&lag_dur_vec_1<=max_burst_len;%
burst_ft_2 = feature_sign_vec_2 == 1&lead_dur_vec_2>=min_pause_len&lead_dur_vec_2<=max_pause_len...
    &lag_dur_vec_2>=min_burst_len&lag_dur_vec_2<=max_burst_len;%
% sampling vectors
sample_options_1 = find(burst_ft_1);
sample_options_2 = find(burst_ft_2);


% (1) make de-trended input-output figure with controls
burst_rise_spot_array_dm_1 = NaN(nBoots,n_col);
burst_rise_virt_array_dm_1 = NaN(nBoots,n_col);
burst_rise_spot_array_dm_2 = NaN(nBoots,n_col);
burst_rise_virt_array_dm_2 = NaN(nBoots,n_col);
% take bootstrap samples
for n = 1:nBoots
    % project 1
    s_ids_1 = randsample(sample_options_1,numel(sample_options_1),true);        
    burst_rise_spot_array_dm_1(n,:) = nanmean(spot_array_dm_1(s_ids_1,:));
    burst_rise_virt_array_dm_1(n,:) = nanmean(virtual_array_dm_1(s_ids_1,:));
    % project 2
    s_ids_2 = randsample(sample_options_2,numel(sample_options_2),true);        
    burst_rise_spot_array_dm_2(n,:) = nanmean(spot_array_dm_2(s_ids_2,:));
    burst_rise_virt_array_dm_2(n,:) = nanmean(virtual_array_dm_2(s_ids_2,:));
end

% project #1
% calculate mean and standard error for spot
burst_rise_spot_mean_1 = nanmean(burst_rise_spot_array_dm_1);
burst_rise_spot_ste_1 = nanstd(burst_rise_spot_array_dm_1);

% calculate mean and standard error for virtual spot
burst_rise_virt_mean_1 = nanmean(burst_rise_virt_array_dm_1);
burst_rise_virt_ste_1 = nanstd(burst_rise_virt_array_dm_1);

% project #2
% calculate mean and standard error for spot
burst_rise_spot_mean_2 = nanmean(burst_rise_spot_array_dm_2);
burst_rise_spot_ste_2 = nanstd(burst_rise_spot_array_dm_2);

% calculate mean and standard error for virtual spot
burst_rise_virt_mean_2 = nanmean(burst_rise_virt_array_dm_2);
burst_rise_virt_ste_2 = nanstd(burst_rise_virt_array_dm_2);


%% make figure
cmap1 = brewermap([],'Set2');

burst_dt_fig_virt = figure;
% Dorsal activity
hold on

% locus (project 1)
p1 = plot(time_axis,burst_rise_spot_mean_1,'-','Color',cmap1(2,:),'LineWidth',2);
ylabel('relative Dl concentration (au)')
% locus (project 2)
p2 = plot(time_axis,burst_rise_spot_mean_2,'--','Color',cmap1(2,:),'LineWidth',2);
ylabel('relative Dl concentration (au)')


p = plot(0,0);
ax = gca;
ax.YColor = 'black';%cmap1(2,:);
% grid on
xlabel('offset (minutes)')
legend('Dl at {\it snail} locus (FFF)','Dl at {\it snail} locus (OG)', 'Location','northwest');

set(gca,'Fontsize',14,'xtick',-4:2:4)
chH = get(gca,'Children');
set(gca,'Children',flipud(chH));
set(gca,    'Box','off',...
            'Color',[228,221,209]/255,...            
            'TickLength',[0.02,0.05])    
burst_dt_fig_virt.Color = 'white';        
burst_dt_fig_virt.InvertHardcopy = 'off';
% save
saveas(burst_dt_fig_virt,[FigPath 'locus_trend_comparisons.tif'])
saveas(burst_dt_fig_virt,[FigPath 'locus_trend_comparisons.pdf'])