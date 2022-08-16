% Script to analyze protein feature variation as a function of size and
% duration of transcription bursts
clear
close all
addpath('utilities')
% define core ID variables
project = 'Dl-Ven_snaBAC-mCh';
% project = 'Dl-Ven_hbP2P-mCh';
DropboxFolder =  'E:\Nick\LivemRNA\Dropbox\';
[~, DataPath, FigRoot] =   header_function(DropboxFolder, project);
FigPath = [FigRoot '\' project '\_paper_figures\input_output02\'];
mkdir(FigPath)
% load data
load([DataPath 'hmm_input_output_results.mat'])
w = 7;
K = 3;
load([DataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'],'hmm_input_output')


% burst rise range
if strcmp(project,'Dl-Ven_hbP2P-mCh')
    burst_range = 5:10;
    burst_sigma = 3;
else
    burst_range = 2:12;
    burst_sigma = 3;
end
min_buffer_len = 5;
max_buffer_len = 30;
% extract relevant arrays
lag_dur_vec = results_struct.lag_dur_vec;
lead_dur_vec = results_struct.lead_dur_vec;
feature_sign_vec = results_struct.feature_sign_vec;
hmm_array = results_struct.hmm_array;
swap_hmm_array = results_struct.swap_hmm_array;
spot_array = results_struct.spot_array_dt;
swap_array = results_struct.swap_array_dt;
virtual_array = results_struct.virtual_array_dt;
window_size = size(swap_array,2);
% initialize data arrays
burst_rise_dur_hmm_mean = NaN(numel(burst_range),window_size);
burst_rise_dur_swap_hmm_mean = NaN(numel(burst_range),window_size);
burst_rise_dur_virt_mean = NaN(numel(burst_range),window_size);
burst_rise_dur_spot_mean = NaN(numel(burst_range),window_size);
burst_rise_dur_swap_mean = NaN(numel(burst_range),window_size);
for i = 1:numel(burst_range)
    burst_vec = burst_range(i)-burst_sigma:burst_range(i)+burst_sigma;
    burst_ft = feature_sign_vec == 1 & ismember(lag_dur_vec,burst_vec) & ...
        lead_dur_vec>= min_buffer_len & lead_dur_vec < max_buffer_len;
    % calculate averages
    burst_rise_dur_hmm_mean(i,:) = nanmean(hmm_array(burst_ft,:));
    burst_rise_dur_swap_hmm_mean(i,:) = nanmean(swap_hmm_array(burst_ft,:));
    burst_rise_dur_spot_mean(i,:) = nanmean(spot_array(burst_ft,:));
    burst_rise_dur_swap_mean(i,:) = nanmean(swap_array(burst_ft,:));
    burst_rise_dur_virt_mean(i,:) = nanmean(virtual_array(burst_ft,:));
end

% % initialize data arrays
% burst_fall_dur_hmm_mean = NaN(numel(burst_range),window_size);
% burst_fall_dur_swap_hmm_mean = NaN(numel(burst_range),window_size);
% burst_fall_dur_spot_mean = NaN(numel(burst_range),window_size);
% burst_fall_dur_swap_mean = NaN(numel(burst_range),window_size);
% burst_fall_dur_virt_mean = NaN(numel(burst_range),window_size);
% for i = 1:numel(burst_range)
%     burst_vec = burst_range(i)-burst_sigma:burst_range(i)+burst_sigma;
%     burst_ft = feature_sign_vec == -1 & ismember(lead_dur_vec,burst_vec) & lag_dur_vec >= min_buffer_len ;
%     % calculate averages
%     burst_fall_dur_hmm_mean(i,:) = nanmean(hmm_array(burst_ft,:));
%     burst_fall_dur_swap_hmm_mean(i,:) = nanmean(swap_hmm_array(burst_ft,:));
%     burst_fall_dur_spot_mean(i,:) = nanmean(spot_array(burst_ft,:));
%     burst_fall_dur_swap_mean(i,:) = nanmean(swap_array(burst_ft,:));
%     burst_fall_dur_virt_mean(i,:) = nanmean(virtual_array(burst_ft,:));
% end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%% RISE HEATMAPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
burst_rise_dur_hm = figure;
burst_rise_dur_hm.Name = 'target spot burst rise hmm';
pt_hm_cm = flipud(brewermap([],'RdYlBu'));
colormap(pt_hm_cm)
% colormap(jet(128)/1.1);
pcolor(flipud(burst_rise_dur_spot_mean))
xlabel('distance from burst start (minutes)')
set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)
ylabel('burst duration (min)')
set(gca,'yticklabels',fliplr(round(burst_range/3,1)))
h = colorbar;
caxis([-15 20])
ylabel(h, 'Dorsal levels (au)','FontSize',14)
set(gca,'FontSize', 14);
saveas(burst_rise_dur_hm, [FigPath 'burst_rise_hm_target.tif'])
saveas(burst_rise_dur_hm, [FigPath 'burst_rise_hm_target.pdf'])

swap_rise_dur_hm = figure;
swap_rise_dur_hm.Name = 'swap spot burst rise hmm';
colormap(pt_hm_cm)
pcolor(flipud(burst_rise_dur_swap_mean))
xlabel('distance from burst start (minutes)')
set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)
ylabel('burst duration (min)')
set(gca,'yticklabels',fliplr(round(burst_range/3,1)))
h = colorbar;
ylabel(h, 'Dorsal levels (au)','FontSize',14)
set(gca,'FontSize', 12);
caxis([-15 20])
saveas(swap_rise_dur_hm, [FigPath 'burst_rise_hm_swap.tif'])
saveas(swap_rise_dur_hm, [FigPath 'burst_rise_hm_swap.pdf'])

virt_rise_dur_hm = figure;
virt_rise_dur_hm.Name = 'virtual spot burst rise hmm';
colormap(pt_hm_cm)
% colormap(jet(128)/1.1);
pcolor(flipud(burst_rise_dur_virt_mean))
xlabel('distance from burst start (minutes)')
set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)
ylabel('burst duration (min)')
set(gca,'yticklabels',fliplr(round(burst_range/3,1)))
h = colorbar;
caxis([-15 20])
ylabel(h, 'Dorsal levels (au)','FontSize',14)
set(gca,'FontSize', 12);
saveas(virt_rise_dur_hm, [FigPath 'burst_rise_hm_virtual.tif'])
saveas(virt_rise_dur_hm, [FigPath 'burst_rise_hm_virtual.pdf'])
%%
hmm_rise_dur_hm = figure;
hmm_rise_dur_hm.Name = 'target spot burst rise hmm';
tr_hm_cm = flipud(flipud(brewermap([],'Purples')));
colormap(tr_hm_cm)
pcolor(flipud(burst_rise_dur_hmm_mean))
xlabel('distance from burst start (minutes)')
set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)
ylabel('burst duration (min)')
set(gca,'yticklabels',fliplr(round(burst_range/3,1)))
h = colorbar;
caxis([0 1.5])
ylabel(h, 'sna activity (au)','FontSize',14)
set(gca,'FontSize', 14);
saveas(hmm_rise_dur_hm, [FigPath 'burst_dur_rise_hm_hmm.tif'])
saveas(hmm_rise_dur_hm, [FigPath 'burst_dur_rise_hm_hmm.pdf'])

swap_hmm_rise_dur_hm = figure;
swap_hmm_rise_dur_hm.Name = 'swap spot burst rise hmm';
colormap(tr_hm_cm)
pcolor(flipud(burst_rise_dur_swap_hmm_mean))
xlabel('offset (min)')
set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)
ylabel('burst duration (min)')
set(gca,'yticklabels',fliplr(round(burst_range/3,1)))
h = colorbar;
caxis([.3 1.5])
ylabel(h, 'sna activity (au)','FontSize',14)
set(gca,'FontSize', 12);
saveas(swap_hmm_rise_dur_hm, [FigPath 'swap_burst_dur_rise_hm_hmm.tif'])
saveas(swap_hmm_rise_dur_hm, [FigPath 'swap_burst_dur_rise_hm_hmm.pdf'])

%%
%%%%%%%%%%%%%%%%%%%%%%%%%% RISE WATERFALLS %%%%%%%%%%%%%%%%%%%%%%%%%%
inc = floor(128/numel(burst_range));
pt_cm_rise = brewermap(128,'Reds');
tr_cm_rise = brewermap(128,'Purples');
burst_rise_dur_wt = figure;
index_vec = 1:numel(burst_range);
hold on
for i = 1:numel(burst_range)
    temp = burst_rise_dur_spot_mean;
    temp(index_vec~=i,:) = NaN;
    temp = temp(:,7:end-6);
    w = waterfall(temp-nanmin(temp(:)),repmat((1:numel(burst_range))',1,window_size-12));
    w.FaceColor = pt_cm_rise(1+(i-1)*inc,:);
    w.FaceAlpha = .6;
    w.EdgeColor = 'black';
end
xlabel('offset (min)')
set(gca,'xtick',1:3:window_size-12,'xticklabels',-3:3)    
zlabel('Dorsal levels (au)')    
view(-15,20)
set(gca,'Fontsize',14)
% xlim([-3.5 3.5])
grid on
saveas(burst_rise_dur_wt, [FigPath 'burst_waterfall_target.tif'])
saveas(burst_rise_dur_wt, [FigPath 'burst_waterfall_target.pdf'])



hmm_rise_dur_wt = figure;
index_vec = 1:numel(burst_range);
hold on
for i = 1:numel(burst_range)
    temp = burst_rise_dur_hmm_mean;
    temp(index_vec~=i,:) = NaN;
    temp = temp(:,7:end-6);
    w = waterfall(temp-nanmin(temp(:)),repmat((1:numel(burst_range))',1,window_size-12));
    w.FaceColor = tr_cm_rise(1+(i-1)*inc,:);
    w.FaceAlpha = .8;
    w.EdgeColor = 'black';
end
xlabel('offset (min)')
set(gca,'xtick',1:3:window_size-12,'xticklabels',-3:3)    
zlabel('sna activity (au)')    
view(-15,20)
set(gca,'Fontsize',14)
grid on
saveas(hmm_rise_dur_wt, [FigPath 'burst_dur_rise_waterfall_hmm.tif'])
saveas(hmm_rise_dur_wt, [FigPath 'burst_dur_rise_waterfall_hmm.pdf'])

