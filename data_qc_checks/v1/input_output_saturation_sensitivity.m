% Script to generate figures establishing presence of enrichment bursts at
% start of transcription bursts
clear 
% close all
addpath('utilities')
% set ID variables
DropboxFolder = 'E:\Nick\LivemRNA\Dropbox (Personal)\';
project = 'Dl-Ven_snaBAC-mCh';
[~, DataPath, FigureRoot] =   header_function(DropboxFolder, project); 

FigPath = [FigureRoot '\' project '\data_qc_checks\'];
mkdir(FigPath)

% load data
load([DataPath 'hmm_input_output_results.mat'])
target_results_struct = results_struct;

% list of sets that are at saturating MCP levels
sat_sets = [1 2 4 6];

Tres = 20; % seconds
% extract relevant arrays from target project 
set_id_vec = floor(target_results_struct.particle_id_vec);
lag_dur_vec_target = target_results_struct.lag_dur_vec;
lead_dur_vec_target = target_results_struct.lead_dur_vec;
hmm_array = target_results_struct.hmm_array;
spot_array_dt = target_results_struct.mf_array;
virtual_array_dt = target_results_struct.virtual_array_dt;
feature_sign_vec_target = target_results_struct.feature_sign_vec;
%  determine snip size
n_col = size(spot_array_dt,2);
window_size = floor(n_col/2);
time_axis = (-window_size:window_size)*Tres/60;

% set basic analyisis parameters
nBoots = 100; % number of bootstrap samples to use
min_pause_len = 6; % minimum length of preceding OFF period (in time steps)
min_burst_len = 2;

%%% (1) make basic input-output figure
close all
burst_rise_hmm_array = NaN(nBoots,n_col);
burst_rise_spot_array_full = NaN(nBoots,n_col);
burst_rise_spot_array_sat = NaN(nBoots,n_col);
burst_rise_spot_array_nsat = NaN(nBoots,n_col);
% generate basic filter for target locus and computational controls
burst_ft_primary = feature_sign_vec_target == 1&lead_dur_vec_target>=min_pause_len&lag_dur_vec_target>min_burst_len; % filter for rise events
sample_opts = find(burst_ft_primary);
sample_opts_sat = find(burst_ft_primary & ismember(set_id_vec,sat_sets));
sample_opts_nsat = find(burst_ft_primary & ~ismember(set_id_vec,sat_sets));
% take bootstrap samples
for n = 1:nBoots
    % sample full
    s_ids = randsample(sample_opts,numel(sample_opts),true);    
    burst_rise_spot_array_full(n,:) = nanmean(spot_array_dt(s_ids,:));
    % saturating only
    s_ids_sat = randsample(sample_opts_sat,numel(sample_opts_sat),true);
    burst_rise_spot_array_sat(n,:) = nanmean(spot_array_dt(s_ids_sat,:));
    % nonsaturating only
    s_ids_nsat = randsample(sample_opts_nsat,numel(sample_opts_nsat),true);
    burst_rise_spot_array_nsat(n,:) = nanmean(spot_array_dt(s_ids_nsat,:));
    burst_rise_hmm_array(n,:) = nanmean(hmm_array(s_ids_sat,:));
end
% calculate mean and standard error
burst_rise_hmm_mean = nanmean(burst_rise_hmm_array);
burst_rise_spot_full = nanmean(burst_rise_spot_array_full);
burst_rise_spot_sat = nanmean(burst_rise_spot_array_sat);
burst_rise_spot_nsat = nanmean(burst_rise_spot_array_nsat);
% burst_rise_spot_ste = nanstd(burst_rise_spot_array_full);
% calculate upper and lower bound vectors
% br_spot_ub = burst_rise_spot_dt + burst_rise_spot_ste;
% br_spot_lb = burst_rise_spot_dt - burst_rise_spot_ste;

% make figure
burst_trend_fig = figure;
cmap1 = brewermap([],'Set2');

% snail activity
yyaxis right
p1 = area(time_axis,burst_rise_hmm_mean,'FaceColor',cmap1(end,:),'LineWidth',1.5,'FaceAlpha',.4);
ylabel('snail transcription (au)')
set(gca,'ytick',.2:.1:1.2)
ylim([.2 1.2])
ax = gca;
ax.YColor = 'black';

% Dorsal activity
yyaxis left
hold on
% fill([time_axis fliplr(time_axis)],[br_spot_ub fliplr(br_spot_lb)],cmap1(2,:),'FaceAlpha',.5,'EdgeAlpha',0)
p2 = plot(time_axis,burst_rise_spot_full,'-','Color',cmap1(2,:),'LineWidth',2);
p3 = plot(time_axis,burst_rise_spot_sat,'-','Color',cmap1(1,:),'LineWidth',2);
p4 = plot(time_axis,burst_rise_spot_nsat,'-','Color',cmap1(3,:),'LineWidth',2);
ylabel('relative Dl concentration (au)')
% set(gca,'ytick',-20:4:20)
ax = gca;
ax.YColor = 'black';

grid on
xlabel('offset (minutes)')
legend([p1 p2 p3 p4],'snail transcription','combined','saturating','non-saturating','Location','northwest')
set(gca,'Fontsize',12,'xtick',-4:2:4)
chH = get(gca,'Children');
set(gca,'Children',flipud(chH));
% save
saveas(burst_trend_fig,[FigPath 'de-meaned_locus_trend.tif'])
saveas(burst_trend_fig,[FigPath 'de-meaned_locus_trend.pdf'])

% % (2) make de-trended input-output figure with controls
% 
% burst_rise_spot_array_dt = NaN(nBoots,n_col);
% burst_rise_swap_array_dt = NaN(nBoots,n_col);
% burst_rise_virt_array_dt = NaN(nBoots,n_col);
% burst_rise_bio_array_dt = NaN(nBoots,n_col);
% % take bootstrap samples
% for n = 1:nBoots
%     % primary
%     s_ids_target = randsample(sample_options_target,numel(sample_options_target),true);    
%     burst_rise_spot_array_dt(n,:) = nanmean(spot_array_dt(s_ids_target,:));
%     burst_rise_swap_array_dt(n,:) = nanmean(swap_array_dt(s_ids_target,:));
%     burst_rise_virt_array_dt(n,:) = nanmean(virtual_array_dt(s_ids_target,:));
%     % biological control
%     s_ids_control = randsample(sample_options_control,numel(sample_options_control),true);    
%     burst_rise_bio_array_dt(n,:) = nanmean(biocontrol_array_dt(s_ids_control,:));
% end
% 
% % calculate mean and standard error for spot
% burst_rise_spot_mean = nanmean(burst_rise_spot_array_dt);
% burst_rise_spot_ste = nanstd(burst_rise_spot_array_dt);
% br_spot_ub = burst_rise_spot_mean + burst_rise_spot_ste;
% br_spot_lb = burst_rise_spot_mean - burst_rise_spot_ste;
% % calculate mean and standard error for nn swap
% burst_rise_swap_mean = nanmean(burst_rise_swap_array_dt);
% burst_rise_swap_ste = nanstd(burst_rise_swap_array_dt);
% br_swap_ub = burst_rise_swap_mean + burst_rise_swap_ste;
% br_swap_lb = burst_rise_swap_mean - burst_rise_swap_ste;
% % calculate mean and standard error for virtual spot
% burst_rise_virt_mean = nanmean(burst_rise_virt_array_dt);
% burst_rise_virt_ste = nanstd(burst_rise_virt_array_dt);
% br_virt_ub = burst_rise_virt_mean + burst_rise_virt_ste;
% br_virt_lb = burst_rise_virt_mean - burst_rise_virt_ste;
% % calculate mean and standard error for virtual spot
% burst_rise_bio_mean = nanmean(burst_rise_bio_array_dt);
% burst_rise_bio_ste = nanstd(burst_rise_bio_array_dt);
% br_bio_ub = burst_rise_bio_mean + burst_rise_bio_ste;
% br_bio_lb = burst_rise_bio_mean - burst_rise_bio_ste;
% 
% % make figure
% burst_dt_fig = figure;
% % snail activity
% yyaxis right
% p1 = area(time_axis,burst_rise_hmm_mean,'FaceColor',cmap1(end,:),'LineWidth',1,'FaceAlpha',.25);
% ylabel('snail transcription (au)')
% set(gca,'ytick',.2:.1:1.2)
% ylim([.2 1.2])
% ax = gca;
% ax.YColor = 'black';
% % Dorsal activity
% yyaxis left
% hold on
% % virtual control
% fill([time_axis fliplr(time_axis)],[br_virt_ub fliplr(br_virt_lb)],cmap1(3,:),'FaceAlpha',.15,'EdgeAlpha',0)
% p2 = plot(time_axis,burst_rise_virt_mean,'-','Color',cmap1(3,:),'LineWidth',1.5);
% % swap control
% fill([time_axis fliplr(time_axis)],[br_swap_ub fliplr(br_swap_lb)],cmap1(5,:),'FaceAlpha',.15,'EdgeAlpha',0)
% p3 = plot(time_axis,burst_rise_swap_mean,'-','Color',cmap1(5,:),'LineWidth',1.5);
% %biological control
% fill([time_axis fliplr(time_axis)],[br_bio_ub fliplr(br_bio_lb)],cmap1(6,:),'FaceAlpha',.15,'EdgeAlpha',0)
% p4 = plot(time_axis,burst_rise_bio_mean,'-','Color',cmap1(6,:),'LineWidth',1.5);
% %locus
% fill([time_axis fliplr(time_axis)],[br_spot_ub fliplr(br_spot_lb)],cmap1(2,:),'FaceAlpha',.5,'EdgeAlpha',0)
% p5 = plot(time_axis,burst_rise_spot_mean,'-','Color',cmap1(2,:),'LineWidth',2);
% ylabel('relative Dl concentration (au)')
% % set(gca,'ytick',-12:3:24)
% % ylim([-12 18])
% ax = gca;
% ax.YColor = 'black';
% grid on
% xlabel('offset (minutes)')
% legend([p1 p2 p3 p4 p5],'snail transcription','virtual spot',...
%     'nearest neighbor','off-target ({\ithbP2P})','target ({\itsnail})',...
%     'Location','northwest')
% set(gca,'Fontsize',12,'xtick',-4:2:4)
% chH = get(gca,'Children');
% set(gca,'Children',flipud(chH));
% % save
% saveas(burst_dt_fig,[FigPath 'locus_trend_w_controls.tif'])
% saveas(burst_dt_fig,[FigPath 'locus_trend_w_controls.pdf'])