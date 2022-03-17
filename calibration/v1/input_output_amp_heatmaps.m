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


%%% burst amp range
amp_val_vec = lag_size_vec(feature_sign_vec==1);
n_bins = 11;
amp_range = linspace(prctile(amp_val_vec,5),prctile(amp_val_vec,95),n_bins);
amp_sigma = 2*median(diff(amp_range));
min_buffer_len = 5;
max_burst_dur = 30;


% initialize data arrays
n_boots = 100;
burst_size_hmm_mean = NaN(numel(amp_range),window_size,n_boots);
burst_size_spot_mean = NaN(numel(amp_range),window_size,n_boots);

for i = 1:numel(amp_range)
    amp_vec = [amp_range(i)-amp_sigma amp_range(i)+amp_sigma];
    burst_ft = feature_sign_vec == 1 & lag_size_vec >= amp_vec(1) & lag_size_vec < amp_vec(2) & ...
        lead_dur_vec>= min_buffer_len & lag_dur_vec <= max_burst_dur;
    burst_indices = find(burst_ft);
    for n = 1:n_boots
        boot_burst_indices = randsample(burst_indices,numel(burst_indices),true);        
        % calculate averages
        burst_size_hmm_mean(i,:,n) = nanmean(hmm_array(boot_burst_indices,:));  
        burst_size_spot_mean(i,:,n) = nanmean(spot_array(boot_burst_indices,:));  
    end
end
burst_size_spot_mean = nanmean(burst_size_spot_mean,3);
burst_size_hmm_mean = nanmean(burst_size_hmm_mean,3);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%% RISE HEATMAPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
% Set the x-axis (time from start of burst) limits
time_vec = linspace(-5,5,window_size);
time_lb = -2;   % [min], start x axis here
xlim_lb = find(time_vec == time_lb);
time_ub = 4;
xlim_ub = find(time_vec == time_ub);   % [min], end x axis here

% protein channel
burst_amp_dur_hm = figure;
burst_amp_dur_hm.Name = 'target spot burst rise hmm';
pt_hm_cm = flipud(brewermap([],'RdBu'));
colormap(pt_hm_cm)
pcolor(flipud(burst_size_spot_mean(:,xlim_lb:xlim_ub)))
axis equal tight
xlabel('time from burst start (minutes)')
set(gca,'xtick',1:3:(xlim_ub - xlim_lb + 1),'xticklabels',[time_lb:time_ub])
ylabel('{\itsna} burst amplitude (au/min)')
set(gca,'ytick',3:3:numel(amp_range),'yticklabels',round(3*fliplr(amp_range(3:3:numel(amp_range))),1))    %***HARD-CODED***
c = colorbar;
caxis([-.28 .28])
c.Ticks = round(linspace(-.28,.28,5),2);
ylabel(c, 'Dorsal enrichment (au)','FontSize',14)
set(gca,'FontSize', 14);
saveas(burst_amp_dur_hm, [FigPath 'burst_amp_hm_protein.tif'])
saveas(burst_amp_dur_hm, [FigPath 'burst_amp_hm_protein.pdf'])


% transcription channel
hmm_rise_dur_hm = figure;
hmm_rise_dur_hm.Name = 'target spot burst rise hmm';
tr_hm_cm = brewermap([],'Greys');
colormap(tr_hm_cm)
pcolor(flipud(burst_size_hmm_mean(:,xlim_lb:xlim_ub)))
axis equal tight
xlabel('time from burst start (minutes)')
set(gca,'xtick',1:3:(xlim_ub - xlim_lb + 1),'xticklabels',[time_lb:time_ub])
ylabel('{\itsna} burst amplitude (au/min)')
set(gca,'ytick',3:3:numel(amp_range),'yticklabels',round(3*fliplr(amp_range(3:3:numel(amp_range))),1))    %***HARD-CODED***
c = colorbar;
caxis([0 1.5])
ylabel(c, '{\itsna} transcriptional activity (au)','FontSize',14)
set(gca,'FontSize', 14);
saveas(hmm_rise_dur_hm, [FigPath 'burst_rise_hm_hmm.tif'])
saveas(hmm_rise_dur_hm, [FigPath 'burst_rise_hm_hmm.pdf'])

%%
%%%%%%%%%%%%%%%%%%%%%%%%%% RISE WATERFALLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % protein waterfall
% inc = floor(128/numel(amp_range));
% pt_cm_rise = brewermap(128,'Reds');
% tr_cm_rise = brewermap(128,'Purples');
% burst_rise_dur_wt = figure;
% index_vec = 1:numel(amp_range);
% hold on
% for i = 1:numel(amp_range)
%     temp = burst_size_spot_mean;
%     temp(index_vec~=i,:) = NaN;
%     temp = temp(:,7:end-6);
%     w = waterfall(temp-nanmin(temp(:)),repmat((1:numel(amp_range))',1,window_size-12));
%     w.FaceColor = pt_cm_rise(1+(i-1)*inc,:);
%     w.FaceAlpha = .6;
%     w.EdgeColor = 'black';
% end
% xlabel('offset (min)')
% set(gca,'xtick',1:3:window_size-12,'xticklabels',-3:3)    
% zlabel('Dorsal levels (au)')    
% view(-15,20)
% set(gca,'Fontsize',14)
% % xlim([-3.5 3.5])
% grid on
% saveas(burst_rise_dur_wt, [FigPath 'burst_waterfall_target.tif'])
% saveas(burst_rise_dur_wt, [FigPath 'burst_waterfall_target.pdf'])
% 
% % transcription waterfall
% hmm_rise_dur_wt = figure;
% index_vec = 1:numel(amp_range);
% hold on
% for i = 1:numel(amp_range)
%     temp = burst_size_hmm_mean;
%     temp(index_vec~=i,:) = NaN;
%     temp = temp(:,7:end-6);
%     w = waterfall(temp-nanmin(temp(:)),repmat((1:numel(amp_range))',1,window_size-12));
%     w.FaceColor = tr_cm_rise(1+(i-1)*inc,:);
%     w.FaceAlpha = .8;
%     w.EdgeColor = 'black';
% end
% xlabel('offset (min)')
% set(gca,'xtick',1:3:window_size-12,'xticklabels',-3:3)    
% zlabel('sna activity (au)')    
% view(-15,20)
% set(gca,'Fontsize',14)
% grid on
% saveas(hmm_rise_dur_wt, [FigPath 'burst_dur_rise_waterfall_hmm.tif'])
% saveas(hmm_rise_dur_wt, [FigPath 'burst_dur_rise_waterfall_hmm.pdf'])
% 
% 
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%% AREA PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time_lb = -1;   % [min], start x axis here
% xlim_lb = find(time_vec == time_lb);
% time_ub = 4;    % [min], end x axis here
% xlim_ub = find(time_vec == time_ub);   
% 
% timeFromBurst = xlim_lb:xlim_ub;
% 
% % Making a square representation of the average transcription signal for
% % visualization/presentation purposes only
% zeroIndex = find(time_vec == 0);
% durationCohorts = [3 6 9];
% durationTimes = [find(time_vec == 1), find(time_vec == 2), find(time_vec == 3)];
% cohortLabels = ["short bursts (1 min)", "medium bursts (2 min)", "long bursts (3 min)"];
% burst_rise_dur_hmm_square = zeros(length(durationCohorts),window_size);
% for i = 1:numel(durationCohorts)
%     burst_rise_dur_hmm_square(durationCohorts(i),zeroIndex:durationTimes(i)) = burst_size_hmm_mean(durationCohorts(i),durationTimes(i));
% end
% 
% for i = 1:numel(durationCohorts)
% %     burstDur_hmmSquare = burst_rise_dur_hmm_square(durationCohorts(i),xlim_lb:xlim_ub);
%     burstDur_hmm = burst_size_hmm_mean(durationCohorts(i),xlim_lb:xlim_ub);
%     burstDur_hmm = burstDur_hmm - nanmin(burstDur_hmm);
% %     burstDur_hmmSquare = burstDur_hmmSquare - nanmin(burstDur_hmmSquare);
%     burstDur_protein = burst_size_spot_mean(durationCohorts(i),xlim_lb:xlim_ub);
%     burstDur_protein_min = burstDur_protein - nanmin(nanmin(burst_size_spot_mean));
% %     burstDur_protein(durationCohorts(i),1) = 0;
% %     burstDur_protein(durationCohorts(i),end) = 0;
% 
%     burstFig = figure;
%     hold on
%     yyaxis right
%     proteinAreaPlot = area(timeFromBurst, burstDur_protein_min);
%     ylabel('Dorsal enrichment (au)')
%     ylim([0 40])
%     yticks(linspace(0,40,8))
%     yticklabels(string(linspace(-15,20,8)))
%     set(proteinAreaPlot,'FaceColor',[213,108,85]/255,'FaceAlpha',0.4)
%     yyaxis left
%     StandardFigurePBoC(proteinAreaPlot, gca)
%     hmmAreaPlot = area(timeFromBurst, burstDur_hmm);
% %     hmmAreaPlot = area(timeFromBurst, burstDur_hmmSquare);
%     ylabel('{\itsna} transcription (au)')
%     ylim([0 1])
%     yticks(linspace(0,1,3))
%     yticklabels(string(linspace(0,1,3)))
%     set(hmmAreaPlot,'FaceColor',[115,142,193]/255)
%     xlabel('time from start of burst (min)')
%     xlim([xlim_lb,xlim_ub])
%     xticks(xlim_lb:3:xlim_ub)
%     xticklabels(string(time_lb:time_ub))
%     title(cohortLabels(i))
%     hold off
%     StandardFigurePBoC(hmmAreaPlot, gca)
%     
%     saveas(burstFig, [FigPath 'burstDurRise_hmmSquareProtein_FromMin' num2str(i) '.pdf'])
%     saveas(burstFig, [FigPath 'burstDurRise_hmmSquareProtein_FromMin' num2str(i) '.tif'])
%     
%     
%     % Alternatively, Using the actual zero value as the "zero" point
%     burstDur_protein_zero = burstDur_protein.*(burstDur_protein >= 0);
%     
%     burstFigFromZero = figure;
%     hold on
%     yyaxis right
%     proteinAreaPlot = area(timeFromBurst, burstDur_protein_zero);
%     ylabel('Dorsal enrichment (au)')
%     ylim([0 20])
%     yticks(linspace(0,20,3))
%     yticklabels(string(linspace(0,20,3)))
%     set(proteinAreaPlot,'FaceColor',[213,108,85]/255,'FaceAlpha',0.4)
%     yyaxis left
%     StandardFigurePBoC(proteinAreaPlot, gca)
%     hmmAreaPlot = area(timeFromBurst, burstDur_hmm);
% %     hmmAreaPlot = area(timeFromBurst, burstDur_hmmSquare);
%     ylabel('{\itsna} transcription (au)')
%     ylim([0 1])
%     yticks(linspace(0,1,3))
%     yticklabels(string(linspace(0,1,3)))
%     set(hmmAreaPlot,'FaceColor',[115,142,193]/255)
%     xlabel('time from start of burst (min)')
%     xlim([xlim_lb,xlim_ub])
%     xticks(xlim_lb:3:xlim_ub)
%     xticklabels(string(time_lb:time_ub))
%     title(cohortLabels(i))
%     hold off
%     StandardFigurePBoC(hmmAreaPlot, gca)
%     
%     saveas(burstFigFromZero, [FigPath 'burstDurRise_hmmSquareProtein_FromZero' num2str(i) '.pdf'])
%     saveas(burstFigFromZero, [FigPath 'burstDurRise_hmmSquareProtein_FromeZero' num2str(i) '.tif'])
% end
% 
