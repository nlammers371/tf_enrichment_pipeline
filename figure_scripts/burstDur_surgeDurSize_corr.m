%Script to quantify the relationship between burst duration and surge
%duration & magnitude
clear
close all
addpath('utilities')
% define core ID variables
project = 'Dl-Ven_snaBAC-mCh';
DropboxFolder =  'E:\Meghan\Dropbox\';
[~, ~, FigRoot] =   header_function(DropboxFolder, project);
DropboxFolder = 'E:\Nick\LivemRNA\Dropbox (Personal)\';
[~, DataPath, ~] =   header_function(DropboxFolder, project);
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

% Define time vector
time_vec_data = linspace(-5,5,window_size);

% Define trimmed time window, if neeaded
time_lb = -1;   % [min], start x axis here
xlim_lb_data = find(time_vec_data == time_lb);
xlim_lb_interp = find(time_vec_interp == time_lb);
time_ub = 4;    % [min], end x axis here
xlim_ub_data = find(time_vec_data == time_ub); 
xlim_ub_interp = find(time_vec_interp == time_ub);

%% Bootstrap & interpolate the data
%Bootstrap initialization
n_boots = 100;
burstDur_hmm_boots = NaN(numel(burst_range),window_size,n_boots);
burstDur_protein_boots = NaN(numel(burst_range),window_size,n_boots);

% Interpolation initialization
%define the x axis (time) vector for the interpolations
time_points_factor_increase = 6;   %increase time resolution by this factor
num_time_points_interp = (window_size - 1)*time_points_factor_increase + 1;
time_vec_interp = linspace(-5,5,num_time_points_interp);
burstDur_protein_boots_interp = NaN(numel(burst_range),num_time_points_interp,n_boots);

%Surge duration & magnitude quantification initialization
protein_cum_dist_100p = nan(numDurCohorts,n_boots);
protein_cum_dist_68p = nan(numDurCohorts,n_boots);
protein_surge_dur = nan(numDurCohorts,n_boots);
protein_max_index = nan(numDurCohorts,n_boots);
gaussianStd = [0.6827, 0.8664, 0.9545, 0.9973]; %1, 1.5, 2, and 3 std from mean
deliveryPercent = gaussianStd(1);

interpMethod = 'pchip';

for i = 1:numel(burst_range)
   burst_vec = burst_range(i)-burst_sigma:burst_range(i)+burst_sigma;
   burst_filter = feature_sign_vec == 1 & ismember(lag_dur_vec,burst_vec) & ...
       lead_dur_vec>= min_buffer_len & lead_dur_vec < max_buffer_len;
   burst_indices = find(burst_filter);
   
   for n = 1:n_boots
        boot_burst_indices = randsample(burst_indices,numel(burst_indices),true);

        % Calculate averages ...
        % for sample data
        burstDur_hmm_boots(i,:,n) = nanmean(hmm_array(boot_burst_indices,:));
        burstDur_protein_boots(i,:,n) = nanmean(spot_array(boot_burst_indices,:));
        % for interpolated data
        protein_interp = interp1(time_vec_data,burstDur_protein_boots(i,:,n),time_vec_interp,'pchip'); % different methods you can use: 'linear', 'pchip', 'makima', 'spline', etc.
        burstDur_protein_boots_interp(i,:,n) = protein_interp;

        %Calculate surge duration & magnitude
        burstDur_protein_zero = burstDur_protein_boots_interp(i,:,n).*(burstDur_protein_boots_interp(i,:,n) >= 0);  %Taking everything above zero (the average trendline at the locus)
        % Might want to trim already here (or not?)
        protein_cum_dist_100p(i,n) = trapz(time_vec_interp,burstDur_protein_zero);
        protein_max_index(i,n) = find(burstDur_protein_zero == max(burstDur_protein_zero));
        for j = 1:round(numel(burstDur_protein_zero)/2)-1
            startIndex = protein_max_index(i,n) - j;
            endIndex = protein_max_index(i,n) + j;
            if startIndex <= 0 || endIndex > numel(burstDur_protein_zero)
                error('Percentage is too large given the location of the maximum. Use a lower %.')
            end
            partialArea = trapz(time_vec_interp_trim(startIndex:endIndex),burstDur_protein_zero(startIndex:endIndex));
            percentArea = partialArea/totalEnrich_zero(i);
            %Exit for loop once you've hit the desired percent.
            if percentArea >= deliveryPercent
                deliveryTime_zero(i) = (time_vec_interp_trim(endIndex) - time_vec_interp_trim(startIndex));
                break
            end      
        end
   end
end
burstDur_hmm_mean = nanmean(burstDur_hmm_boots,3);
burstDur_hmm_ste = nanstd(burstDur_hmm_boots,0,3);
burstDur_protein_mean = nanmean(burstDur_protein_boots,3);
burstDur_protein_ste = nanstd(burstDur_protein_boots,0,3);
burstDur_protein_interp_mean = nanmean(burstDur_protein_boots_interp,3);
burstDur_protein_interp_ste = nanstd(burstDur_protein_boots_interp,0,3);

currDurCohort = 2;
% plot(time_vec_data,burstDur_protein_mean(currDurCohort,:),'o',...
%      time_vec_interp,burstDur_protein_interp_mean(currDurCohort,:),':.')
hold on
e_interp = errorbar(time_vec_data,burstDur_protein_mean(currDurCohort,:),burstDur_protein_ste(currDurCohort,:),'Color','black','LineWidth',1.5);
e.CapSize = 0;
% legend('Sample Points','pchip')
hold off


% Trim down time window as needed
time_vec_data_trim = time_vec_data(xlim_lb_data:xlim_ub_data);
time_vec_interp_trim = time_vec_interp(xlim_lb_interp:xlim_ub_interp);
burstDur_protein_mean_trim = burstDur_protein_mean(:,xlim_lb_data:xlim_ub_data);
burstDur_protein_ste_trim = burstDur_protein_ste(:,xlim_lb_data:xlim_ub_data);
burstDur_protein_interp_mean_trim = burstDur_protein_interp_mean(:,xlim_lb_interp:xlim_ub_interp);
burstDur_protein_interp_ste_trim = burstDur_protein_interp_ste(:,xlim_lb_interp:xlim_ub_interp);

durationCohorts = [2 3 4 5 6 7 8 9 10 11];
numDurCohorts = numel(durationCohorts);
burstDuration = [0.66 1 1.33 1.66 2 2.33 2.66 3 3.33 3.66];

%% Quantifying Surge Duration

totalEnrich_min = nan(1,numDurCohorts);
partialEnrich_min = nan(numDurCohorts,numel(time_vec_interp_trim));
deliveryTime_min = nan(1,numDurCohorts);

totalEnrich_zero = nan(1,numDurCohorts);
partialEnrich_zero = nan(numDurCohorts,numel(time_vec_interp_trim));
deliveryTime_zero = nan(1,numDurCohorts);

interpMethod = 'pchip';

%Calculate period of time over which x% of all the protein is "delivered" 
%to the locus, centered around the max enrichment

for i = 1:numel(durationCohorts)
    burstDur_protein = burstDur_protein_interp_mean_trim(durationCohorts(i),:);
    
    % Using the minimum recorded enrichment across all cohorst as the 
    % "zero" point
    burstDur_protein_min = burstDur_protein - nanmin(nanmin(burstDur_protein_interp_mean_trim(:,:)));
    totalEnrich_min(i) = trapz(time_vec_interp_trim,burstDur_protein_min);
    maxEnrichIndex_min = find(burstDur_protein_min == max(burstDur_protein_min));
    for j = 1:round(numel(burstDur_protein_min)/2)-1
        startIndex = maxEnrichIndex_min - j;
        endIndex = maxEnrichIndex_min + j;
        if startIndex <= 0 || endIndex > numel(burstDur_protein_min)
            error('Percentage is too large given the location of the maximum. Use a lower %.')
        end
        
        partialArea = trapz(time_vec_interp_trim(startIndex:endIndex),burstDur_protein_min(startIndex:endIndex));
        percentArea = partialArea/totalEnrich_min(i);
        %Exit for loop once you've hit the desired percent.
        if percentArea >= deliveryPercent
            deliveryTime_min(i) = (time_vec_interp_trim(endIndex) - time_vec_interp_trim(startIndex));
            break
        end 

     
    end
    
    % Using the actual zero value as the "zero" point
    
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
saveas(surgeDur_min_fig, [FigPath 'burstDur_vs_surgeDur_fromMin_' interpMethod num2str(i) '.pdf'])
saveas(surgeDur_min_fig, [FigPath 'burstDur_vs_surgeDur_fromMin_' interpMethod num2str(i) '.tif'])

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
saveas(surgeDur_zero_fig, [FigPath 'burstDur_vs_surgeDur_fromZero_' interpMethod num2str(i) '.pdf'])
saveas(surgeDur_zero_fig, [FigPath 'burstDur_vs_surgeDur_fromZero_' interpMethod num2str(i) '.tif'])

%%
%Calculate amount of time surge spends above zero (the average trend of Dl
%at the locus)

maxIndex = nan(1,numDurCohorts);
surgeStartIndex = nan(1,numDurCohorts);
surgeEndIndex = nan(1,numDurCohorts);
timeAboveZero = nan(1,numDurCohorts);

for i = 1:numel(durationCohorts)
    burstDur_protein_mean = burstDur_protein_interp_mean_trim(durationCohorts(i),:);
    burstDur_protein_ste = burstDur_protein_interp_ste_trim(durationCohorts(i),:);
    
    [maxProtein,maxIndex(i)] = max(burstDur_protein_mean);
    surgeStartIndex(i) = find(burstDur_protein_mean > 0, 1, 'first');   % First time surge goes above trend
    
    burstDur_protein_temp = burstDur_protein_mean(maxIndex(i):end);
    burstDur_protein_ste_temp = burstDur_protein_ste(maxIndex(i):end);
    surgeEndIndex(i) = maxIndex(i) + find((burstDur_protein_temp - burstDur_protein_ste_temp) <= 0, 1, 'first') - 1;   %First time surge dips within the ste of the trend
    
    timeAboveZero(i) = (time_vec_interp_trim(surgeEndIndex(i)) - time_vec_interp_trim(surgeStartIndex(i)));
end

timeAboveZero_fig = figure;
plot_zero = plot(burstDuration, timeAboveZero,'-ko','MarkerSize',40,'MarkerFaceColor',[213,108,85]/256,'MarkerEdgeColor','black');
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
saveas(timeAboveZero_fig, [FigPath 'burstDur_vs_surgeDur_timeAboveZero_' interpMethod num2str(i) '.pdf'])
saveas(timeAboveZero_fig, [FigPath 'burstDur_vs_surgeDur_timeAboveZero_' interpMethod num2str(i) '.tif'])