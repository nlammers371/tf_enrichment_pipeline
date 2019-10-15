% Script to analyze protein feature variation as a function of size and
% duration of transcription bursts
clear
close all
addpath('utilities')
% define core ID variables
project = 'Dl-Ven_snaBAC-mCh';
dropboxFolder =  'E:\Nick\LivemRNA\Dropbox (Personal)\';
dataPath = [dropboxFolder 'ProcessedEnrichmentData\' project '\'];
figPath = [dropboxFolder 'LocalEnrichmentFigures\_paper_figures\input_output02\'];
mkdir(figPath)
% load data
load([dataPath 'hmm_input_output_results.mat'])
w = 7;
K = 3;
load([dataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'],'hmm_input_output')


% burst rise range
burst_range = 2:12;
burst_sigma = 3;
min_buffer_len = 5;
max_buffer_len = 30;

% extract relevant arrays
lag_dur_vec = results_struct.lag_dur_vec;
lead_dur_vec = results_struct.lead_dur_vec;
feature_sign_vec = results_struct.feature_sign_vec;
hmm_array = results_struct.hmm_array; % transcriptional activity at target locus
spot_array = results_struct.spot_array_dt; % protein snips at target locus
window_size = size(spot_array,2);

% initialize data arrays
burst_rise_dur_hmm_mean = NaN(numel(burst_range),window_size);
burst_rise_dur_spot_mean = NaN(numel(burst_range),window_size);

for i = 1:numel(burst_range)
    burst_vec = burst_range(i)-burst_sigma:burst_range(i)+burst_sigma;
    burst_ft = feature_sign_vec == 1 & ismember(lag_dur_vec,burst_vec) & ...
        lead_dur_vec>= min_buffer_len & lead_dur_vec < max_buffer_len;
    % calculate averages
    burst_rise_dur_hmm_mean(i,:) = nanmean(hmm_array(burst_ft,:));  
    burst_rise_dur_spot_mean(i,:) = nanmean(spot_array(burst_ft,:));  
end

%%%%%%%%%%%%%%%%%%%%%%%%%% RISE HEATMAPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% protein channel
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
saveas(burst_rise_dur_hm, [figPath 'burst_rise_hm_protein.tif'])
saveas(burst_rise_dur_hm, [figPath 'burst_rise_hm_protein.pdf'])


% transcription channel
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
saveas(hmm_rise_dur_hm, [figPath 'burst_rise_hm_hmm.tif'])
saveas(hmm_rise_dur_hm, [figPath 'burst_rise_hm_hmm.pdf'])


%%%%%%%%%%%%%%%%%%%%%%%%%% RISE WATERFALLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% protein waterfall
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
saveas(burst_rise_dur_wt, [figPath 'burst_waterfall_target.tif'])
saveas(burst_rise_dur_wt, [figPath 'burst_waterfall_target.pdf'])

% transcription waterfall
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
saveas(hmm_rise_dur_wt, [figPath 'burst_dur_rise_waterfall_hmm.tif'])
saveas(hmm_rise_dur_wt, [figPath 'burst_dur_rise_waterfall_hmm.pdf'])
