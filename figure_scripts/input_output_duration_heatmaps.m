% Script to analyze protein feature variation as a function of size and
% duration of transcription bursts
clear
close all
addpath('utilities')
% define core ID variables
project = 'Dl-Ven_snaBAC-mCh_v3';
% DropboxFolder =  'E:\Meghan\Dropbox\';
DropboxFolder = 'E:\Nick\LivemRNA\Dropbox (Personal)\';
[~, DataPath, FigRoot] =   header_function(DropboxFolder, project);
FigPath = [FigRoot '\_paper_figures\burst_analyses\'];
mkdir(FigPath)

% load data
load([DataPath 'hmm_input_output_results.mat'])
% w = 7;
% K = 3;
% load([DataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'],'hmm_input_output')

% extract relevant arrays
lag_dur_vec = results_struct.lag_dur_vec;
lag_size_vec = results_struct.lag_size_vec;
lead_dur_vec = results_struct.lead_dur_vec;
feature_sign_vec = results_struct.feature_sign_vec;
hmm_array = results_struct.hmm_array; % transcriptional activity at target locus
spot_array = results_struct.spot_array_dt; % protein snips at target locus
window_size = size(spot_array,2);


% burst rise range
burst_range = 2:12;
burst_sigma = 2;
min_buffer_len = 5;
max_burst_dur = 30;


% initialize data arrays
n_boots = 100;
burst_dur_hmm_array = NaN(numel(burst_range),window_size,n_boots);
burst_dur_spot_array = NaN(numel(burst_range),window_size,n_boots);

for i = 1:numel(burst_range)
    burst_vec = [burst_range(i)-burst_sigma burst_range(i)+burst_sigma];
    burst_ft = feature_sign_vec == 1 & lag_dur_vec >= burst_vec(1) & lag_dur_vec <= burst_vec(2) & ...
        lead_dur_vec>= min_buffer_len & lag_dur_vec < max_burst_dur;
    burst_indices = find(burst_ft);
    for n = 1:n_boots
        boot_burst_indices = randsample(burst_indices,numel(burst_indices),true);        
        % calculate averages
        burst_dur_hmm_array(i,:,n) = nanmean(hmm_array(boot_burst_indices,:));  
        burst_dur_spot_array(i,:,n) = nanmean(spot_array(boot_burst_indices,:));  
    end
end
burst_dur_spot_mean = nanmean(burst_dur_spot_array,3);
burst_dur_hmm_mean = nanmean(burst_dur_hmm_array,3);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%% RISE HEATMAPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the x-axis (time from start of burst) limits
time_vec = linspace(-5,5,window_size);
time_lb = -2;   % [min], start x axis here
xlim_lb = find(time_vec == time_lb);
time_ub = 4;
xlim_ub = find(time_vec == time_ub);   % [min], end x axis here

% protein channel
burst_rise_dur_hm = figure;
burst_rise_dur_hm.Name = 'target spot burst rise hmm';
pt_hm_cm = flipud(brewermap([],'RdBu'));
colormap(pt_hm_cm)
pcolor(flipud(burst_dur_spot_mean(:,xlim_lb:xlim_ub)))
xlabel('time from burst start (minutes)')
set(gca,'xtick',1:3:(xlim_ub - xlim_lb + 1),'xticklabels',[time_lb:time_ub])
ylabel('{\itsna} burst duration (min)')
set(gca,'ytick',3:3:(burst_range(end) - burst_range(1) +1),'yticklabels',fliplr([1 2 3]))    %***HARD-CODED***
c = colorbar;
caxis([-.3 .3])
c.Ticks = round(linspace(-.35,.35,11),2);
ylabel(c, 'Dorsal enrichment (au)','FontSize',14)
set(gca,'FontSize', 14);
saveas(burst_rise_dur_hm, [FigPath 'burst_dur_hm_protein.tif'])
saveas(burst_rise_dur_hm, [FigPath 'burst_dur_hm_protein.pdf'])


% transcription channel
hmm_rise_dur_hm = figure;
hmm_rise_dur_hm.Name = 'target spot burst rise hmm';
tr_hm_cm = flipud(flipud(brewermap([],'Greys')));
colormap(tr_hm_cm)
pcolor(flipud(burst_dur_hmm_mean(:,xlim_lb:xlim_ub)))
xlabel('time from burst start (minutes)')
set(gca,'xtick',1:3:(xlim_ub - xlim_lb + 1),'xticklabels',[time_lb:time_ub])
ylabel('{\itsna} burst duration (min)')
set(gca,'ytick',3:3:(burst_range(end) - burst_range(1) +1),'yticklabels',fliplr([1 2 3]))	%***HARD-CODED***
c = colorbar;
caxis([0 1.5])
ylabel(c, '{\itsna} transcriptional activity (au)','FontSize',14)
set(gca,'FontSize', 14);
saveas(hmm_rise_dur_hm, [FigPath 'burst_dur_hm_hmm.tif'])
saveas(hmm_rise_dur_hm, [FigPath 'burst_dur_hm_hmm.pdf'])

