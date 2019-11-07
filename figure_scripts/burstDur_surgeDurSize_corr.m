%Script to quantify the relationship between burst duration and surge
%duration & magnitude
clear
close all
addpath('utilities')
% define core ID variables
project = 'Dl-Ven_snaBAC-mCh';
% DropboxFolder =  'E:\Meghan\Dropbox\';
DropboxFolder = 'E:\Nick\LivemRNA\Dropbox (Personal)\';
[~, DataPath, FigRoot] =   header_function(DropboxFolder, project);
FigPath = [FigRoot '\_paper_figures\burstDur_surgeSizeDur_corr\'];
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

time_vec_data = linspace(-5,5,window_size);

% Define the x axis (time) vector for the interpolations
time_points_factor_increase = 6;   %increase time resolution by this factor
num_time_points_interp = (window_size - 1)*time_points_factor_increase + 1;
time_vec_interp = linspace(-5,5,num_time_points_interp);

% Trim down time window as desired
time_lb = -1;   % [min], start x axis here
xlim_lb_data = find(time_vec_data == time_lb);
xlim_lb_interp = find(time_vec_interp == time_lb);
time_ub = 4;    % [min], end x axis here
xlim_ub_data = find(time_vec_data == time_ub); 
xlim_ub_interp = find(time_vec_interp == time_ub);

time_vec_data_trim = time_vec_data(xlim_lb_data:xlim_ub_data);
time_vec_interp_trim = time_vec_interp(xlim_lb_interp:xlim_ub_interp);
burst_rise_dur_spot_mean_trim = burst_rise_dur_spot_mean(:,xlim_lb_data:xlim_ub_data);

durationCohorts = [2 3 4 5 6 7 8 9 10 11];
numDurCohorts = numel(durationCohorts);
burstDuration = [0.66 1 1.33 1.66 2 2.33 2.66 3 3.33 3.66];

%% Interpolate protein channel data
interpMethod = {'linear', 'pchip', 'makima'}; %different interpolation methods to try
burstDur_protein_interp = nan(size(burst_rise_dur_spot_mean_trim,1),numel(time_vec_interp_trim),numel(interpMethod));

for i = 1:numel(interpMethod)
    for j = 1:size(burst_rise_dur_spot_mean_trim,1)
        interpResult = interp1(time_vec_data_trim,burst_rise_dur_spot_mean_trim(j,:),time_vec_interp_trim,interpMethod{i});
        burstDur_protein_interp(j,:,i) = interpResult;
    end 
end

currDurCohort = 9;

plot(time_vec_data_trim,burst_rise_dur_spot_mean_trim(currDurCohort,:),'o',...
     time_vec_interp_trim,burstDur_protein_interp(currDurCohort,:,1),':.', ...
     time_vec_interp_trim,burstDur_protein_interp(currDurCohort,:,2),'-.',...
     time_vec_interp_trim,burstDur_protein_interp(currDurCohort,:,3),'-')
legend('Sample Points','linear','pchip','makima','spline')



%% Quantifying Surge Duration
close all

totalEnrich_min = nan(1,numDurCohorts);
partialEnrich_min = nan(numDurCohorts,numel(timeFromBurst));
percentEnrich_min = nan(numDurCohorts,numel(timeFromBurst));
fwhmMinZero = nan(1,numDurCohorts);
deliveryTime_min = nan(1,numDurCohorts);

totalEnrich_zero = nan(1,numDurCohorts);
partialEnrich_zero = nan(numDurCohorts,numel(timeFromBurst));
percentEnrich_zero = nan(numDurCohorts,numel(timeFromBurst));
deliveryTime_zero = nan(1,numDurCohorts);

for i = 1:numel(durationCohorts)
%     burstDur_hmmSquare = burst_rise_dur_hmm_square(durationCohorts(i),xlim_lb:xlim_ub);
%     burstDur_hmm = burst_rise_dur_hmm_mean(durationCohorts(i),xlim_lb:xlim_ub);
%     burstDur_hmm = burstDur_hmm - nanmin(burstDur_hmm);
%     burstDur_hmmSquare = burstDur_hmmSquare - nanmin(burstDur_hmmSquare);
    burstDur_protein = burst_rise_dur_spot_mean(durationCohorts(i),xlim_lb_data:xlim_ub_data);
    
    % Using the minimum recorded enrichment across all cohorst as the 
    % "zero" point
    burstDur_protein_min = burstDur_protein - nanmin(nanmin(burst_rise_dur_spot_mean));
    %Calculate period of time over which x% of all the protein is 
    %"delivered" to the locus, centered around the max enrichment
    totalEnrich_min(i) = trapz(timeFromBurst,burstDur_protein_min);
    maxEnrichIndex_min = find(burstDur_protein_min == max(burstDur_protein_min));
    deliveryPercent = 0.68;
    for j = 1:round(numel(burstDur_protein_min)/2)-1
        startIndex = maxEnrichIndex_min - j;
        endIndex = maxEnrichIndex_min + j;
        if startIndex <= 0 || endIndex > numel(burstDur_protein_min)
            error('Percentage is too large given the location of the maximum. Use a lower %.')
        end
        partialArea = trapz(timeFromBurst(startIndex:endIndex),burstDur_protein_min(startIndex:endIndex));
        percentArea = partialArea/totalEnrich_min(i);
        %Exit for loop once you've hit the desired percent.
        if percentArea >= deliveryPercent
            deliveryTime_min(i) = (timeFromBurst(endIndex) - timeFromBurst(startIndex)) / 3;
            break
        end      
    end
    
    % Using the actual zero value as the "zero" point
    burstDur_protein_zero = burstDur_protein.*(burstDur_protein >= 0);
    %Calculate period of time over which x% of all the protein is 
    %"delivered" to the locus, centered around the max enrichment
    totalEnrich_zero(i) = trapz(timeFromBurst,burstDur_protein_zero);
    maxEnrichIndex_zero = find(burstDur_protein_zero == max(burstDur_protein_zero));
    deliveryPercent = 0.68;
    for j = 1:round(numel(burstDur_protein_min)/2)-1
        startIndex = maxEnrichIndex_zero - j;
        endIndex = maxEnrichIndex_zero + j;
        if startIndex <= 0 || endIndex > numel(burstDur_protein_zero)
            error('Percentage is too large given the location of the maximum. Use a lower %.')
        end
        partialArea = trapz(timeFromBurst(startIndex:endIndex),burstDur_protein_zero(startIndex:endIndex));
        percentArea = partialArea/totalEnrich_zero(i);
        %Exit for loop once you've hit the desired percent.
        if percentArea >= deliveryPercent
            deliveryTime_zero(i) = (timeFromBurst(endIndex) - timeFromBurst(startIndex)) / 3;
            break
        end      
    end
    
%     for j = 2:numel(timeFromBurst)
%         partialArea = trapz(timeFromBurst(1:j),burstDur_protein_zeroed(1:j));
%         percentArea = partialArea/totalEnrichMinZero(i);
%         partialEnrichMinZero(i,j-1) = partialArea;
%         percentEnrichMinZero(i,j-1) = percentArea;
%     end
%     %Calculate period of time during which 68% to 95% of all the protein is
%     %"delivered" to the locus
%     startPercent = 0.68;
%     endPercent = 0.95;
%     dtStartIndex = find(percentEnrichMinZero(i,:) >= startPercent,1,'first');
%     dtEndIndex = find(percentEnrichMinZero(i,:) >= endPercent,1,'first');
%     deliveryTimeMinZero(i) = (timeFromBurst(dtEndIndex) - timeFromBurst(dtStartIndex)) / 3;
    
%     %Calculate full width half max (fwhm)
%     halfMax = max(burstDur_protein_zeroed)/2;
%     hmStartIndex = find(burstDur_protein_zeroed >= halfMax, 1, 'first');
%     hmEndIndex = find(burstDur_protein_zeroed >= halfMax, 1, 'last');
%     fwhmMinZero(i) = (timeFromBurst(hmEndIndex) - timeFromBurst(hmStartIndex)) / 3;
   
    
end

surgeDur_min_fig = figure;
% plot(burstDuration, fwhmMinZero);
% hold on
plot_min = plot(burstDuration, deliveryTime_min,'-o','MarkerFaceColor',[213,108,85]/256,'MarkerEdgeColor','black');
% hold off
box on
ylabel('Dl surge duration (min)')
xlabel('{\itsna} burst duration (min)')
ylim([0 4])
yticks(linspace(0,4,5))
yticklabels(string(linspace(0,4,5)))
xlim([0.66 3.66])
xticks(linspace(0,4,5))
xticklabels(string(linspace(0,4,5)))
set(gca,'FontSize',14)
% legend('Delivery Time')
StandardFigurePBoC(plot_min, gca)
saveas(surgeDur_min_fig, [FigPath 'burstDur_vs_surgeDur_fromMin' num2str(i) '.pdf'])
saveas(surgeDur_min_fig, [FigPath 'burstDur_vs_surgeDur_fromMin' num2str(i) '.tif'])

surgeDur_zero_fig = figure;
% plot(burstDuration, fwhmMinZero);
% hold on
plot_zero = plot(burstDuration, deliveryTime_zero,'-ko','MarkerSize',40,'MarkerFaceColor',[213,108,85]/256,'MarkerEdgeColor','black');
% hold off
box on
ylabel('Dl surge duration (min)')
xlabel('{\itsna} burst duration (min)')
ylim([0 4])
yticks(linspace(0,4,5))
yticklabels(string(linspace(0,4,5)))
xlim([0.66 3.66])
xticks(linspace(0,4,5))
xticklabels(string(linspace(0,4,5)))
set(gca,'FontSize',14)
% legend('Delivery Time')
StandardFigurePBoC(plot_zero, gca)
saveas(surgeDur_zero_fig, [FigPath 'burstDur_vs_surgeDur_fromZero' num2str(i) '.pdf'])
saveas(surgeDur_zero_fig, [FigPath 'burstDur_vs_surgeDur_fromZero' num2str(i) '.tif'])