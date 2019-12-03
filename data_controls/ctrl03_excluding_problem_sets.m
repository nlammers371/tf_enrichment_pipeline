% Script to generate figures establishing presence of enrichment bursts at
% start of transcription bursts
clear 
% close all
addpath('utilities')
% set ID variables
targetProject = 'Dl-Ven_snaBAC-mCh';
controlProject = 'Dl-Ven_hbP2P-mCh';
DropboxFolder = 'E:\Nick\LivemRNA\Dropbox (Personal)\';
[~, DataPathTarget, FigureRoot] =   header_function(DropboxFolder, targetProject); 
[~, DataPathControl, ~] =   header_function(DropboxFolder, controlProject); 


% load data
load([DataPathTarget 'hmm_input_output_results.mat'])
target_results_struct = results_struct;
load([DataPathControl 'hmm_input_output_results.mat'])
control_results_struct = results_struct;


FigPath = [FigureRoot 'data_controls\'];
mkdir(FigPath)
Tres = 20; % seconds
% extract relevant arrays from target project 

lag_dur_vec_target = target_results_struct.lag_dur_vec;
lead_dur_vec_target = target_results_struct.lead_dur_vec;

spot_array_dt = target_results_struct.spot_array_dt;
target_mf_vec = target_results_struct.mf_protein_vec;
feature_sign_vec_target = target_results_struct.feature_sign_vec;
target_set_vec = floor(target_results_struct.particle_id_vec);
% extract arrays from control project
hmm_array = control_results_struct.hmm_array;
lag_dur_vec_control = control_results_struct.lag_dur_vec;
lead_dur_vec_control = control_results_struct.lead_dur_vec;
biocontrol_spot_array_dt = control_results_struct.spot_array_dt;
biocontrol_mf_vec = control_results_struct.mf_protein_vec;
feature_sign_vec_control = control_results_struct.feature_sign_vec;
control_set_vec = floor(control_results_struct.particle_id_vec);
%  determine snip size
n_col = size(spot_array_dt,2);
window_size = floor(n_col/2);
time_axis = (-window_size:window_size)*Tres/60;

% set basic analyisis parameters
nBoots = 100; % number of bootstrap samples to use
min_pause_len = 6; % minimum length of preceding OFF period (in time steps)
min_burst_len = 2;

% generate basic filter for target locus and computational controls
burst_ft_primary = feature_sign_vec_target == 1&lead_dur_vec_target>=min_pause_len&lag_dur_vec_target>min_burst_len & target_set_vec~=4; % filter for rise events
burst_ft_control = feature_sign_vec_control == 1&lead_dur_vec_control>=min_pause_len&lag_dur_vec_control>min_burst_len & control_set_vec~=3; % filter for rise events


%%% check for potential effects from disparate sampel sizes
close all
burst_rise_hmm_array = NaN(nBoots,n_col);
burst_rise_spot_array = NaN(nBoots,n_col);
burst_rise_bio_array = NaN(nBoots,n_col);

sample_options_target = find(burst_ft_primary);
sample_options_control = find(burst_ft_control);
% take bootstrap samples
for n = 1:nBoots
    s_ids_target = randsample(sample_options_target,numel(sample_options_control),true);    
    s_ids_control = randsample(sample_options_control,numel(sample_options_control),true);
    %
    burst_rise_hmm_array(n,:) = nanmean(hmm_array(s_ids_control,:));    
    burst_rise_spot_array(n,:) = nanmean(spot_array_dt(s_ids_target,:));
    burst_rise_bio_array(n,:) = nanmean(biocontrol_spot_array_dt(s_ids_control,:));
end
% calculate mean and standard error
burst_rise_hmm_mean = nanmean(burst_rise_hmm_array);
burst_rise_spot_mean = nanmean(burst_rise_spot_array);
burst_rise_spot_ste = nanstd(burst_rise_spot_array);
burst_rise_bio_mean = nanmean(burst_rise_bio_array);
burst_rise_bio_ste = nanstd(burst_rise_bio_array);


%% make figure
burst_trend_fig = figure;
cmap1 = brewermap([],'Set2');

% snail activity
yyaxis right
p1 = area(time_axis,burst_rise_hmm_mean,'FaceColor',cmap1(end,:),'LineWidth',1.5,'FaceAlpha',.4);
ylabel('snail transcription (au)')
% set(gca,'ytick',.2:.1:1.2)
% ylim([.2 1.2])
ax = gca;
ax.YColor = 'black';

% Dorsal activity
yyaxis left
hold on
% fill([time_axis fliplr(time_axis)],[br_spot_ub fliplr(br_spot_lb)],cmap1(2,:),'FaceAlpha',.5,'EdgeAlpha',0)
p2 = plot(time_axis,burst_rise_spot_mean,'-','Color',cmap1(2,:),'LineWidth',2);
p3 = plot(time_axis,burst_rise_bio_mean,'-','Color',cmap1(6,:),'LineWidth',2);
ylabel('relative Dl concentration (au)')
% set(gca,'ytick',-20:4:20)
ax = gca;
ax.YColor = 'black';

grid on
xlabel('offset (minutes)')
legend([p1 p2 p3],'snail transcription','Dl (target)','Dl (control)','Location','northwest')
set(gca,'Fontsize',12,'xtick',-4:2:4)
chH = get(gca,'Children');
set(gca,'Children',flipud(chH));

% % save
% saveas(burst_trend_fig,[FigPath 'downsampling_locus_control.tif'])
% saveas(burst_trend_fig,[FigPath 'downsampling_locus_control.pdf'])
