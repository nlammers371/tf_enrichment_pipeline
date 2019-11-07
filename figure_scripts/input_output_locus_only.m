% Script to analyze protein feature variation as a function of size and
% duration of transcription bursts
clear
close all
addpath('utilities')
% define core ID variables
project = 'Dl-Ven_snaBAC-mCh';
% DropboxFolder =  'E:\Meghan\Dropbox\';
DropboxFolder = 'E:\Nick\LivemRNA\Dropbox (Personal)\';
[~, DataPath, FigRoot] =   header_function(DropboxFolder, project);
FigPath = [FigRoot '\_paper_figures\input_output04\'];
mkdir(FigPath)

% load data
load([DataPath 'hmm_input_output_results.mat'])
w = 7;
K = 3;
load([DataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'],'hmm_input_output')

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
pcolor(flipud(burst_rise_dur_spot_mean(:,xlim_lb:xlim_ub)))
xlabel('time from burst start (minutes)')
set(gca,'xtick',1:3:(xlim_ub - xlim_lb + 1),'xticklabels',[time_lb:time_ub])
ylabel('{\itsna} transcription burst duration (min)')
set(gca,'ytick',3:3:(burst_range(end) - burst_range(1) +1),'yticklabels',fliplr([1 2 3]))    %***HARD-CODED***
c = colorbar;
caxis([-20 20])
c.Ticks = linspace(-20,20,5);
ylabel(c, 'Dorsal enrichment (au)','FontSize',14)
set(gca,'FontSize', 14);
saveas(burst_rise_dur_hm, [FigPath 'burst_rise_hm_protein.tif'])
saveas(burst_rise_dur_hm, [FigPath 'burst_rise_hm_protein.pdf'])


% transcription channel
hmm_rise_dur_hm = figure;
hmm_rise_dur_hm.Name = 'target spot burst rise hmm';
tr_hm_cm = flipud(flipud(brewermap([],'Greys')));
colormap(tr_hm_cm)
pcolor(flipud(burst_rise_dur_hmm_mean(:,xlim_lb:xlim_ub)))
xlabel('time from burst start (minutes)')
set(gca,'xtick',1:3:(xlim_ub - xlim_lb + 1),'xticklabels',[time_lb:time_ub])
ylabel('{\itsna} transcription burst duration (min)')
set(gca,'ytick',3:3:(burst_range(end) - burst_range(1) +1),'yticklabels',fliplr([1 2 3]))	%***HARD-CODED***
c = colorbar;
caxis([0 1.5])
ylabel(c, '{\itsna} transcriptional activity (au)','FontSize',14)
set(gca,'FontSize', 14);
saveas(hmm_rise_dur_hm, [FigPath 'burst_rise_hm_hmm.tif'])
saveas(hmm_rise_dur_hm, [FigPath 'burst_rise_hm_hmm.pdf'])

%%
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
saveas(burst_rise_dur_wt, [FigPath 'burst_waterfall_target.tif'])
saveas(burst_rise_dur_wt, [FigPath 'burst_waterfall_target.pdf'])

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
saveas(hmm_rise_dur_wt, [FigPath 'burst_dur_rise_waterfall_hmm.tif'])
saveas(hmm_rise_dur_wt, [FigPath 'burst_dur_rise_waterfall_hmm.pdf'])


%%
%%%%%%%%%%%%%%%%%%%%%%%%%% AREA PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_lb = -1;   % [min], start x axis here
xlim_lb = find(time_vec == time_lb);
time_ub = 4;    % [min], end x axis here
xlim_ub = find(time_vec == time_ub);   

timeFromBurst = xlim_lb:xlim_ub;

% Making a square representation of the average transcription signal for
% visualization/presentation purposes only
zeroIndex = find(time_vec == 0);
durationCohorts = [3 6 9];
durationTimes = [find(time_vec == 1), find(time_vec == 2), find(time_vec == 3)];
cohortLabels = ["short bursts (1 min)", "medium bursts (2 min)", "long bursts (3 min)"];
burst_rise_dur_hmm_square = zeros(length(durationCohorts),window_size);
for i = 1:numel(durationCohorts)
    burst_rise_dur_hmm_square(durationCohorts(i),zeroIndex:durationTimes(i)) = burst_rise_dur_hmm_mean(durationCohorts(i),durationTimes(i));
end

for i = 1:numel(durationCohorts)
%     burstDur_hmmSquare = burst_rise_dur_hmm_square(durationCohorts(i),xlim_lb:xlim_ub);
    burstDur_hmm = burst_rise_dur_hmm_mean(durationCohorts(i),xlim_lb:xlim_ub);
    burstDur_hmm = burstDur_hmm - nanmin(burstDur_hmm);
%     burstDur_hmmSquare = burstDur_hmmSquare - nanmin(burstDur_hmmSquare);
    burstDur_protein = burst_rise_dur_spot_mean(durationCohorts(i),xlim_lb:xlim_ub);
    burstDur_protein_min = burstDur_protein - nanmin(nanmin(burst_rise_dur_spot_mean));
%     burstDur_protein(durationCohorts(i),1) = 0;
%     burstDur_protein(durationCohorts(i),end) = 0;

    burstFig = figure;
    hold on
    yyaxis right
    proteinAreaPlot = area(timeFromBurst, burstDur_protein_min);
    ylabel('Dorsal enrichment (au)')
    ylim([0 40])
    yticks(linspace(0,40,8))
    yticklabels(string(linspace(-15,20,8)))
    set(proteinAreaPlot,'FaceColor',[213,108,85]/255,'FaceAlpha',0.4)
    yyaxis left
    StandardFigurePBoC(proteinAreaPlot, gca)
    hmmAreaPlot = area(timeFromBurst, burstDur_hmm);
%     hmmAreaPlot = area(timeFromBurst, burstDur_hmmSquare);
    ylabel('{\itsna} transcription (au)')
    ylim([0 1])
    yticks(linspace(0,1,3))
    yticklabels(string(linspace(0,1,3)))
    set(hmmAreaPlot,'FaceColor',[115,142,193]/255)
    xlabel('time from start of burst (min)')
    xlim([xlim_lb,xlim_ub])
    xticks(xlim_lb:3:xlim_ub)
    xticklabels(string(time_lb:time_ub))
    title(cohortLabels(i))
    hold off
    StandardFigurePBoC(hmmAreaPlot, gca)
    
    saveas(burstFig, [FigPath 'burstDurRise_hmmSquareProtein_FromMin' num2str(i) '.pdf'])
    saveas(burstFig, [FigPath 'burstDurRise_hmmSquareProtein_FromMin' num2str(i) '.tif'])
    
    
    % Alternatively, Using the actual zero value as the "zero" point
    burstDur_protein_zero = burstDur_protein.*(burstDur_protein >= 0);
    
    burstFigFromZero = figure;
    hold on
    yyaxis right
    proteinAreaPlot = area(timeFromBurst, burstDur_protein_zero);
    ylabel('Dorsal enrichment (au)')
    ylim([0 20])
    yticks(linspace(0,20,3))
    yticklabels(string(linspace(0,20,3)))
    set(proteinAreaPlot,'FaceColor',[213,108,85]/255,'FaceAlpha',0.4)
    yyaxis left
    StandardFigurePBoC(proteinAreaPlot, gca)
    hmmAreaPlot = area(timeFromBurst, burstDur_hmm);
%     hmmAreaPlot = area(timeFromBurst, burstDur_hmmSquare);
    ylabel('{\itsna} transcription (au)')
    ylim([0 1])
    yticks(linspace(0,1,3))
    yticklabels(string(linspace(0,1,3)))
    set(hmmAreaPlot,'FaceColor',[115,142,193]/255)
    xlabel('time from start of burst (min)')
    xlim([xlim_lb,xlim_ub])
    xticks(xlim_lb:3:xlim_ub)
    xticklabels(string(time_lb:time_ub))
    title(cohortLabels(i))
    hold off
    StandardFigurePBoC(hmmAreaPlot, gca)
    
    saveas(burstFigFromZero, [FigPath 'burstDurRise_hmmSquareProtein_FromZero' num2str(i) '.pdf'])
    saveas(burstFigFromZero, [FigPath 'burstDurRise_hmmSquareProtein_FromeZero' num2str(i) '.tif'])
end

