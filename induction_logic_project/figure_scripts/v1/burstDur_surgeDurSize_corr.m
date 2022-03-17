%Script to quantify the relationship between burst duration and surge
%duration & magnitude
clear
close all
addpath('utilities')
% define core ID variables
project = 'Dl-Ven_snaBAC-mCh_v3';
DropboxFolder =  'S:\Nick\Dropbox\';
[~, ~, FigRoot] =   header_function(DropboxFolder, project);
DropboxFolder = 'S:\Nick\Dropbox\';
[~, DataPath, ~] =   header_function(DropboxFolder, project);
FigPath = [FigRoot '\_paper_figures\burstFeature_vs_surgeFeature\'];
mkdir(FigPath)

% load data
load([DataPath 'hmm_input_output_results.mat'])
% w = 7;
% K = 3;
% load([DataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'],'hmm_input_output')

% Extract relevant arrays from results structure
lag_dur_vec = results_struct.lag_dur_vec'; %burst duration
lead_dur_vec = results_struct.lead_dur_vec'; %burst duration
lag_amp_vec = results_struct.lag_size_vec'; %burst amplitude
lead_amp_vec = results_struct.lead_size_vec'; %burst amplitude
feature_sign_vec = results_struct.feature_sign_vec';
hmm_array = results_struct.hmm_array; %transcriptional activity at target locus
locus_protein_array = results_struct.spot_array_dt; %protein snips at target locus

%Filter out bursts where all protein values are NaN
protein_filter = all(isnan(locus_protein_array),2); %finds All-NaN rows
locus_protein_array = locus_protein_array(~protein_filter,:);
lag_dur_vec = lag_dur_vec(~protein_filter);
lead_dur_vec = lead_dur_vec(~protein_filter);
lag_amp_vec = lag_amp_vec(~protein_filter);
lead_amp_vec = lead_amp_vec(~protein_filter);
feature_sign_vec = feature_sign_vec(~protein_filter);
hmm_array = hmm_array(~protein_filter);

% Burst duration range & sigma
burst_dur_vec = lag_dur_vec(feature_sign_vec == 1);
burst_dur_range = 2:12;
burst_dur_range_inMin = linspace(0.6667,4,numel(burst_dur_range));
min_buffer_length = 5;
max_buffer_length = 30;
burst_dur_sigma = 3;

% Burst amplitude range & sigma
burst_amp_vec = lag_amp_vec(feature_sign_vec == 1);
n_bins = 15;
burst_amp_range = linspace(prctile(burst_amp_vec,5),prctile(burst_amp_vec,95),n_bins);
burst_amp_sigma = median(diff(burst_amp_range));

% Burst filters
burst_filter = (feature_sign_vec == 1) & ...
                    lead_dur_vec >= min_buffer_length & ...
                    lead_dur_vec < max_buffer_length;

% Define time vectors
window_size = size(locus_protein_array,2);  %total number of time points included
abs_window_size = 15; %number of time points included before and after burst start
time_zero_index = abs_window_size + 1; % t = 0
time_vec_data = linspace(-5,5,window_size); %[min]
time_step = 20/60;  %[sec]
% Define the x axis (time) vector for the interpolations
time_points_factor_increase = 6;   %increase time resolution by this factor
num_time_points_interp = (window_size - 1)*time_points_factor_increase + 1;
time_vec_interp = linspace(-5,5,num_time_points_interp);
% Define trimmed time window, if needed
time_lb = -1;   % [min], start x axis here
xlim_lb_data = find(time_vec_data == time_lb);
xlim_lb_interp = find(time_vec_interp == time_lb);
time_ub = 4;    % [min], end x axis here
xlim_ub_data = find(time_vec_data == time_ub); 
xlim_ub_interp = find(time_vec_interp == time_ub);
xaxis_burst_start = find(time_vec_interp == 0);
% Trim down time window as needed
time_vec_data_trim = time_vec_data(xlim_lb_data:xlim_ub_data);
time_vec_interp_trim = time_vec_interp(xlim_lb_interp:xlim_ub_interp);

% Plotting variables
PBoC_red = [213,108,85]/256;
PBoC_blue = [115 143 193]/256;
PBoC_green = [122 169 116]/256;
PBoC_yellow = [234 194 100]/256;
PBoC_purple = [171 133 172]/256;

%% Quantifying Surge Duration

length_burst_dur_range = numel(burst_dur_range);

%different methods you can use for interp1: 'linear', 'pchip', 'makima', 
%'spline', etc.
protein_interp = interp1(time_vec_data,locus_protein_array',time_vec_interp,'pchip');
protein_interp = protein_interp';

% burstDur_protein_boots_interp = NaN(numel(burst_dur_range),num_time_points_interp,n_boots);
% 
% %Surge duration & magnitude quantification initialization
% protein_cum_dist_100p = nan(length_burst_dur_range,n_boots);
% protein_cum_dist_68p = nan(length_burst_dur_range,n_boots);
% protein_surge_dur = nan(length_burst_dur_range,n_boots);
% % protein_max_index = nan(length_burst_range,n_boots);
% protein_surge_amp_2min = nan(length_burst_dur_range,n_boots);
% gaussianStd = [0.6827, 0.8664, 0.9545, 0.9973]; %1, 1.5, 2, and 3 std from mean
% percent_delivered = 0.5; %gaussianStd(1);

%Calculate surge duration & magnitude
%         burstDur_protein_trim = burstDur_protein_boots_interp(i,xlim_lb_interp:xlim_ub_interp,n);
%         burstDur_protein_zero = burstDur_protein_trim.*(burstDur_protein_trim >= 0);  %Taking everything above zero (the average trendline at the locus)
%         % Might want to trim already here (or not?)
%         protein_cum_dist_100p(i,n) = trapz(time_vec_interp_trim,burstDur_protein_zero);
%         protein_max_index(i,n) = find(burstDur_protein_zero == max(burstDur_protein_zero));
%         for j = 1:round(numel(burstDur_protein_zero)/2)-1
%             startIndex = protein_max_index(i,n) - j;
%             endIndex = protein_max_index(i,n) + j;
%             if startIndex <= 0 || endIndex > numel(burstDur_protein_zero)
%                 warning('Percentage is too large given the location of the maximum. Use a lower %.')
%                 protein_surge_dur(i,n) = NaN;
%                 break
%             end
%             protein_cum_dist_68p(i,n) = trapz(time_vec_interp_trim(startIndex:endIndex),burstDur_protein_zero(startIndex:endIndex));
%             frac_cum_dist = protein_cum_dist_68p(i,n) / protein_cum_dist_100p(i,n);
%             %Exit for loop once you've hit the desired percent.
%             if frac_cum_dist >= percent_delivered
%                 protein_surge_dur(i,n) = (time_vec_interp_trim(endIndex) - time_vec_interp_trim(startIndex));
%                 break
%             end      
%         end


%% Explore different ways of defining surge magnitude

% Make locus protein zero point the value at start of burst (t = 0)
locus_protein_time_zero = locus_protein_array(:,time_zero_index);
locus_protein_array_zeroed = locus_protein_array - locus_protein_time_zero;

% Surge magnitude = total protein delivered in fixed 2 min window after start of burst
est_burst_dur = 2; %[min]
est_surge_window = est_burst_dur / time_step;
surge_mag_vec_2min_fixed = nansum(locus_protein_array(:,time_zero_index:time_zero_index + est_surge_window),2);
surge_mag_vec_zeroed_2min_fixed = nansum(locus_protein_array_zeroed(:,time_zero_index:time_zero_index + est_surge_window),2);

% Surge magnitude = total protein delivered during actual duration of burst
burst_window = lag_dur_vec;
surge_mag_vec_actual_burst_dur = nan(numel(lag_dur_vec),1);
for i = 1:numel(lag_dur_vec)
    if (time_zero_index + burst_window(i)) <= window_size
        surge_mag_vec_actual_burst_dur(i) = nansum(locus_protein_array_zeroed(i,time_zero_index:time_zero_index + burst_window(i)),2) / burst_window(i);
%         surge_mag_vec_actual_burst_dur(i) = nansum(locus_protein_array(i,time_zero_index:time_zero_index + burst_window(i)),2) / burst_window(i);
    else
      surge_mag_vec_actual_burst_dur(i) = NaN;
    end
end

% Surge magnitude = Rate of protein delivery in 1st [20, 40, 60] seconds of
% burst
rate_window = [1 2 3]; %20, 40, 60 seconds
locus_protein_array_20s = locus_protein_array_zeroed(:,time_zero_index:time_zero_index + rate_window(1));
% locus_protein_array_20s = locus_protein_array(:,time_zero_index:time_zero_index + rate_window(1));
surge_mag_vec_rate_20s = (locus_protein_array_20s(:,end) - locus_protein_array_20s(:,1)) ./ rate_window(1)*time_step;

locus_protein_array_40s = locus_protein_array_zeroed(:,time_zero_index:time_zero_index + rate_window(2));
% locus_protein_array_40s = locus_protein_array(:,time_zero_index:time_zero_index + rate_window(2));
surge_mag_vec_rate_40s = (locus_protein_array_40s(:,end) - locus_protein_array_40s(:,1)) ./ rate_window(2)*time_step;

locus_protein_array_60s = locus_protein_array_zeroed(:,time_zero_index:time_zero_index + rate_window(3));
% locus_protein_array_60s = locus_protein_array(:,time_zero_index:time_zero_index + rate_window(3));
surge_mag_vec_rate_60s = (locus_protein_array_60s(:,end) - locus_protein_array_60s(:,1)) ./ rate_window(3)*time_step;

%% Burst amplitude vs surge magnitude

% Bootstrap data
length_burst_amp_range = numel(burst_amp_range);
n_boots = 100;
% burst_amp_hmm_boots = NaN(length_burst_amp_range,window_size,n_boots);
% burst_amp_protein_boots = NaN(length_burst_amp_range,window_size,n_boots);
% burst_amp_hmm_boots = NaN(n_boots,length_burst_amp_range);
burst_amp_surgeMag2min_boots = NaN(n_boots,length_burst_amp_range);
burst_amp_surgeMag2minZeroed_boots = NaN(n_boots,length_burst_amp_range);
burst_amp_surgeMagActBurstDur_boots = NaN(n_boots,length_burst_amp_range);
burst_amp_surgeMagRate20s_boots = NaN(n_boots,length_burst_amp_range);
burst_amp_surgeMagRate40s_boots = NaN(n_boots,length_burst_amp_range);
burst_amp_surgeMagRate60s_boots = NaN(n_boots,length_burst_amp_range);

for a = 1:numel(burst_amp_range)
    % Filter for bursts in the current or directly neighboring amplitude 
    % bins (shouldn't change the weighted average much, doing it just to be
    % cautious)
    amp_sigma_limit = 2;
    burst_amp_bin_limits = [burst_amp_range(a) - amp_sigma_limit*burst_amp_sigma, ...
                            burst_amp_range(a) + amp_sigma_limit*burst_amp_sigma];
    burst_amp_range_filter = ~isnan(lag_amp_vec) & ...
                             lag_amp_vec >= burst_amp_bin_limits(1) & ...
                             lag_amp_vec < burst_amp_bin_limits(2);
    burst_amp_indices = find(burst_filter & burst_amp_range_filter);
    amp_curr = burst_amp_range(a);
    for n = 1:n_boots
        % Get random bootstrap sample
        boots_indices = randsample(burst_amp_indices,numel(burst_amp_indices),true);
        boots_burst_amp = lag_amp_vec(boots_indices);
        % Calc. gaussian weighting
        gauss_wt_vec = exp(-.5*((boots_burst_amp - amp_curr)/burst_amp_sigma).^2); 
        
        % Fixed 2 min surge, raw
        boots_surge_mag_2min = surge_mag_vec_2min_fixed(boots_indices);
        % Take the Gaussian weighted average
        burst_amp_surgeMag2min_boots(n,a) = nansum(gauss_wt_vec .* ...
                                            boots_surge_mag_2min) ./ ...
                                            nansum(gauss_wt_vec);  
        % Fixed 2 min surge, zeroed to start
        boots_surge_mag_2min_zeroed = surge_mag_vec_zeroed_2min_fixed(boots_indices);
        burst_amp_surgeMag2minZeroed_boots(n,a) = nansum(gauss_wt_vec .* ...
                                                boots_surge_mag_2min_zeroed) ./ ...
                                                nansum(gauss_wt_vec);
                                            
        % Actual burst dur
        boots_surge_mag_actual_burst_dur = surge_mag_vec_actual_burst_dur(boots_indices);
        burst_amp_surgeMagActBurstDur_boots(n,a) = nansum(gauss_wt_vec .* ...
                                                boots_surge_mag_actual_burst_dur) ./ ...
                                                nansum(gauss_wt_vec);
        
        % Surge magnitude = Rate of protein delivery in 1st [20, 40, 60] seconds of
        % burst
        boots_surge_mag_vec_rate_20s = surge_mag_vec_rate_20s(boots_indices);
        burst_amp_surgeMagRate20s_boots(n,a) = nansum(gauss_wt_vec .* ...
                                                boots_surge_mag_vec_rate_20s) ./ ...
                                                nansum(gauss_wt_vec);
        boots_surge_mag_vec_rate_40s = surge_mag_vec_rate_40s(boots_indices);
        burst_amp_surgeMagRate40s_boots(n,a) = nansum(gauss_wt_vec .* ...
                                                boots_surge_mag_vec_rate_40s) ./ ...
                                                nansum(gauss_wt_vec);                                    
        boots_surge_mag_vec_rate_60s = surge_mag_vec_rate_60s(boots_indices);
        burst_amp_surgeMagRate60s_boots(n,a) = nansum(gauss_wt_vec .* ...
                                                boots_surge_mag_vec_rate_60s) ./ ...
                                                nansum(gauss_wt_vec);
        
        % same as ^, but average boots first, then calculate slope
%         locus_protein_array_zeroed(:,time_zero_index:time_zero_index + rate_window(1))
    end
end

%Fixed 2 min,, ra
amp_surgeMag2min_mean = nanmean(burst_amp_surgeMag2min_boots);
amp_surgeMag2min_ste = nanstd(burst_amp_surgeMag2min_boots);
%Fixed 2 min, zeroed to start
amp_surgeMag2minZeroed_mean = nanmean(burst_amp_surgeMag2minZeroed_boots);
amp_surgeMag2minZeroed_ste = nanstd(burst_amp_surgeMag2minZeroed_boots);
% Actual burst duration
amp_surgeMagActBurstDur_mean = nanmean(burst_amp_surgeMagActBurstDur_boots);
amp_surgeMagActBurstDur_ste = nanstd(burst_amp_surgeMagActBurstDur_boots);
% Rate at start of burst
amp_surgeMagRate20s_mean = nanmean(burst_amp_surgeMagRate20s_boots);
amp_surgeMagRate20s_ste= nanstd(burst_amp_surgeMagRate20s_boots);
amp_surgeMagRate40s_mean = nanmean(burst_amp_surgeMagRate40s_boots);
amp_surgeMagRate40s_ste= nanstd(burst_amp_surgeMagRate40s_boots);
amp_surgeMagRate60s_mean = nanmean(burst_amp_surgeMagRate60s_boots);
amp_surgeMagRate60s_ste= nanstd(burst_amp_surgeMagRate60s_boots);

% Plot results
%Fixed 2 min, raw
burstAmp_vs_surgeMag_fig = figure(3);
burst_amp_range_auPerMin = burst_amp_range * (1/time_step);
hold on
e_2min = errorbar(burst_amp_range_auPerMin,amp_surgeMag2min_mean,amp_surgeMag2min_ste,'Color','black','LineWidth',1.5);
e_2min.CapSize = 0;
s_2min = scatter(burst_amp_range_auPerMin,amp_surgeMag2min_mean,75,'filled','MarkerFaceColor',PBoC_red);
p = plot(0,0);
box on
xlim([0.95*burst_amp_range_auPerMin(1)  burst_amp_range_auPerMin(end)+burst_amp_range_auPerMin(1)-0.95*burst_amp_range_auPerMin(1)]);
ylim([-0.5 1.5])    %hard-coded to match amp and dur graphs
yticks([-0.5 0.5 1.5])
xlabel('burst amplitude (au/min)')
ylabel('Dorsal surge magnitude (au)')
% legend([s_2min s_2minZeroed],'2min surge dur','2min surge dur, zeroed','Location','northwest')
set(gca,'FontSize',14)
StandardFigure(p,gca)
set(gca,'Color',[228 220 209]/255) 

saveas(burstAmp_vs_surgeMag_fig, [FigPath 'burstAmp_vs_surgeMag_2min' '.pdf'])
saveas(burstAmp_vs_surgeMag_fig, [FigPath 'burstAmp_vs_surgeMag_2min' '.tif'])

%Fixed 2 min window, zeroed to start
burstAmp_vs_surgeMag_zeroed_fig = figure(4);
burst_amp_range_auPerMin = burst_amp_range * (1/time_step);
hold on
e_2minZeroed = errorbar(burst_amp_range_auPerMin,amp_surgeMag2minZeroed_mean,amp_surgeMag2minZeroed_ste,'Color','black','LineWidth',1.5);
e_2minZeroed.CapSize = 0;
s_2minZeroed = scatter(burst_amp_range_auPerMin,amp_surgeMag2minZeroed_mean,75,'filled','MarkerFaceColor',PBoC_blue);
p = plot(0,0);
box on
xlim([0.95*burst_amp_range_auPerMin(1)  burst_amp_range_auPerMin(end)+burst_amp_range_auPerMin(1)-0.95*burst_amp_range_auPerMin(1)]);
ylim([0 2 ])    %hard-coded to match amp and dur graphs
yticks([0 1 2])
xlabel('burst amplitude (au/min)')
ylabel('Dorsal surge magnitude (au)')
% legend([s_2min s_2minZeroed],'2min surge dur','2min surge dur, zeroed','Location','northwest')
set(gca,'FontSize',14)
StandardFigure(p,gca)
set(gca,'Color',[228 220 209]/255) 

saveas(burstAmp_vs_surgeMag_zeroed_fig, [FigPath 'burstAmp_vs_surgeMag_2minZeroed' '.pdf'])
saveas(burstAmp_vs_surgeMag_zeroed_fig, [FigPath 'burstAmp_vs_surgeMag_2minZeroed' '.tif'])


burstAmp_vs_surgeMagActBurstDur_fig = figure(5);
hold on
%Actual burst duration
e_ActBurstDur = errorbar(burst_amp_range_auPerMin,amp_surgeMagActBurstDur_mean,amp_surgeMagActBurstDur_ste,'Color','black','LineWidth',1.5);
e_ActBurstDur.CapSize = 0;
s_ActBurstDur = scatter(burst_amp_range_auPerMin,amp_surgeMagActBurstDur_mean,75,'filled','MarkerFaceColor',PBoC_purple);
p = plot(0,0);
box on
xlim([0.95*burst_amp_range_auPerMin(1)  burst_amp_range_auPerMin(end)+burst_amp_range_auPerMin(1)-0.95*burst_amp_range_auPerMin(1)]);
ylim([0 0.2])     %hard-coded amp and dur graphs match
yticks([0 0.1 0.2])
xlabel('burst amplitude (au/min)')
ylabel('Dorsal surge magnitude (au/min)')
legend(s_ActBurstDur, 'Actual burst dur','Location','northwest')
set(gca,'FontSize',14)
StandardFigure(p,gca)
set(gca,'Color',[228 220 209]/255) 

saveas(burstAmp_vs_surgeMagActBurstDur_fig, [FigPath 'burstAmp_vs_surgeMagActBurstDur' '.pdf'])
saveas(burstAmp_vs_surgeMagActBurstDur_fig, [FigPath 'burstAmp_vs_surgeMagActBurstDur' '.tif'])


% Plot surge mag = rate at start results
burstAmp_vs_surgeRate_fig = figure(6);
burst_amp_range_auPerMin = burst_amp_range * (1/time_step);
hold on
%Rate over first 20s
e_20s = errorbar(burst_amp_range_auPerMin,amp_surgeMagRate20s_mean,amp_surgeMagRate20s_ste,'Color','black','LineWidth',1.5);
e_20s.CapSize = 0;
sca_20s = scatter(burst_amp_range_auPerMin,amp_surgeMagRate20s_mean,75,'filled','MarkerFaceColor',PBoC_red);
%Rate over first 40s
e_40s = errorbar(burst_amp_range_auPerMin,amp_surgeMagRate40s_mean,amp_surgeMagRate20s_ste,'Color','black','LineWidth',1.5);
e_40s.CapSize = 0;
sca_40s = scatter(burst_amp_range_auPerMin,amp_surgeMagRate40s_mean,75,'filled','MarkerFaceColor',PBoC_blue);
%Rate over first 60s
e_60s = errorbar(burst_amp_range_auPerMin,amp_surgeMagRate60s_mean,amp_surgeMagRate20s_ste,'Color','black','LineWidth',1.5);
e_60s.CapSize = 0;
sca_60s = scatter(burst_amp_range_auPerMin,amp_surgeMagRate60s_mean,75,'filled','MarkerFaceColor',PBoC_green);
p = plot(0,0);
box on
xlim([0.95*burst_amp_range_auPerMin(1)  burst_amp_range_auPerMin(end)+burst_amp_range_auPerMin(1)-0.95*burst_amp_range_auPerMin(1)]);
ylim([-0.01 0.06])      %hard-coded amp and dur graphs match
yticks([0 0.03 0.06])
xlabel('burst amplitude (au/min)')
ylabel('Dorsal surge rate (au/min)')
legend([sca_20s sca_40s sca_60s],'20s', '40s','60s','Location','northwest')
set(gca,'FontSize',14)
StandardFigure(p,gca)
set(gca,'Color',[228 220 209]/255) 

saveas(burstAmp_vs_surgeRate_fig, [FigPath 'burstAmp_vs_surgeRate' '.pdf'])
saveas(burstAmp_vs_surgeRate_fig, [FigPath 'burstAmp_vs_surgeRate' '.tif'])

%% Burst duration vs surge magnitude

length_burst_dur_range = numel(burst_dur_range);
n_boots = 100;
burst_dur_hmm_boots = NaN(n_boots,length_burst_dur_range);
burst_dur_surgeMag2min_boots = NaN(n_boots,length_burst_dur_range);
burst_dur_surgeMag2minZeroed_boots = NaN(n_boots,length_burst_dur_range);
burst_dur_surgeMagActBurstDur_boots = NaN(n_boots,length_burst_dur_range);
burst_dur_surgeMagRate20s_boots = NaN(n_boots,length_burst_dur_range);
burst_dur_surgeMagRate40s_boots = NaN(n_boots,length_burst_dur_range);
burst_dur_surgeMagRate60s_boots = NaN(n_boots,length_burst_dur_range);

for d = 1:numel(burst_dur_range)
    % Filter for bursts in the current or directly neighboring amplitude 
    % bins (shouldn't change the weighted average much, doing it just to be
    % cautious)
    dur_sigma_limit = 2;
    burst_dur_bin_limits = [burst_dur_range(d) - dur_sigma_limit*burst_dur_sigma, ...
                            burst_dur_range(d) + dur_sigma_limit*burst_dur_sigma];
    burst_dur_range_filter = ~isnan(lag_dur_vec) & ...
                             lag_dur_vec >= burst_dur_bin_limits(1) & ...
                             lag_dur_vec < burst_dur_bin_limits(2);
    burst_dur_indices = find(burst_filter & burst_dur_range_filter);
    dur_curr = burst_dur_range(d);
    for n = 1:n_boots
        % Pick random bootstrap sample
        boots_indices = randsample(burst_dur_indices,numel(burst_dur_indices),true);
        boots_burst_dur = lag_dur_vec(boots_indices);
        % Calc. gaussian weighting
        gauss_wt_vec = exp(-.5*((boots_burst_dur - dur_curr)/burst_dur_sigma).^2);
        
        % Fixed 2 min, raw
        boots_surge_mag_2min = surge_mag_vec_2min_fixed(boots_indices);
        burst_dur_surgeMag2min_boots(n,d) = nansum(gauss_wt_vec .* boots_surge_mag_2min) ./ nansum(gauss_wt_vec);  
        
        % Fixed 2 min, zeroed to start
        boots_surge_mag_2min_zeroed = surge_mag_vec_zeroed_2min_fixed(boots_indices);
        burst_dur_surgeMag2minZeroed_boots(n,d) = nansum(gauss_wt_vec .* boots_surge_mag_2min_zeroed) ./ nansum(gauss_wt_vec);
        
        % Actual burst dur
        boots_surge_mag_actual_burst_dur = surge_mag_vec_actual_burst_dur(boots_indices);
        burst_dur_surgeMagActBurstDur_boots(n,d) = nansum(gauss_wt_vec .* boots_surge_mag_actual_burst_dur) ./ nansum(gauss_wt_vec);
        
        % Surge magnitude = Rate of protein delivery in 1st [20, 40, 60] seconds of
        % burst
        boots_surge_mag_vec_rate_20s = surge_mag_vec_rate_20s(boots_indices);
        burst_dur_surgeMagRate20s_boots(n,d) = nansum(gauss_wt_vec .* ...
                                                boots_surge_mag_vec_rate_20s) ./ ...
                                                nansum(gauss_wt_vec);
        boots_surge_mag_vec_rate_40s = surge_mag_vec_rate_40s(boots_indices);
        burst_dur_surgeMagRate40s_boots(n,d) = nansum(gauss_wt_vec .* ...
                                                boots_surge_mag_vec_rate_40s) ./ ...
                                                nansum(gauss_wt_vec);                                    
        boots_surge_mag_vec_rate_60s = surge_mag_vec_rate_60s(boots_indices);
        burst_dur_surgeMagRate60s_boots(n,d) = nansum(gauss_wt_vec .* ...
                                                boots_surge_mag_vec_rate_60s) ./ ...
                                                nansum(gauss_wt_vec);
    end
end

%Fixed 2 min, raw
dur_surgeMag2min_mean = nanmean(burst_dur_surgeMag2min_boots);
dur_surgeMag2min_ste = nanstd(burst_dur_surgeMag2min_boots);
% Fixed 2 min, zeroed
dur_surgeMag2minZeroed_mean = nanmean(burst_dur_surgeMag2minZeroed_boots);
dur_surgeMag2minZeroed_ste = nanstd(burst_dur_surgeMag2minZeroed_boots);
% Actual burst duration
dur_surgeMagActBurstDur_mean = nanmean(burst_dur_surgeMagActBurstDur_boots);
dur_surgeMagActBurstDur_ste = nanstd(burst_dur_surgeMagActBurstDur_boots);
% Rate at start of burst
dur_surgeMagRate20s_mean = nanmean(burst_dur_surgeMagRate20s_boots);
dur_surgeMagRate20s_ste= nanstd(burst_dur_surgeMagRate20s_boots);
dur_surgeMagRate40s_mean = nanmean(burst_dur_surgeMagRate40s_boots);
dur_surgeMagRate40s_ste= nanstd(burst_dur_surgeMagRate40s_boots);
dur_surgeMagRate60s_mean = nanmean(burst_dur_surgeMagRate60s_boots);
dur_surgeMagRate60s_ste= nanstd(burst_dur_surgeMagRate60s_boots);

% Fixed 2 min, raw
burstDur_vs_surgeMag_fig = figure(7);
hold on
e_2min = errorbar(burst_dur_range_inMin,dur_surgeMag2min_mean,dur_surgeMag2min_ste,'Color','black','LineWidth',1.5);
e_2min.CapSize = 0;
s_2min = scatter(burst_dur_range_inMin,dur_surgeMag2min_mean,75,'filled','MarkerFaceColor',PBoC_red);
p = plot(0,0);
box on
xlim([.95*burst_dur_range_inMin(1)  burst_dur_range_inMin(end)+burst_dur_range_inMin(1)-.95*burst_dur_range_inMin(1)]);
ylim([-0.5 1.5])    %hard-coded to match amp and dur graphs
yticks([-0.5 0.5 1.5])
xlabel('burst duration (min)')
ylabel('Dorsal surge magnitude (au)')
% legend([s_2min s_2minZeroed],'2min surge dur','2min surge dur, zeroed','Location','northwest')
set(gca,'FontSize',14)
StandardFigure(p,gca)
set(gca,'Color',[228 220 209]/255)

saveas(burstDur_vs_surgeMag_fig, [FigPath 'burstDur_vs_surgeDur_2min' '.pdf'])
saveas(burstDur_vs_surgeMag_fig, [FigPath 'burstDur_vs_surgeDur_2min' '.tif'])

% Fixed 2min, zeroed to start
burstDur_vs_surgeMag_zeroed_fig = figure(8);
hold on
e_2minZeroed = errorbar(burst_dur_range_inMin,dur_surgeMag2minZeroed_mean,dur_surgeMag2minZeroed_ste,'Color','black','LineWidth',1.5);
e_2minZeroed.CapSize = 0;
s_2minZeroed = scatter(burst_dur_range_inMin,dur_surgeMag2minZeroed_mean,75,'filled','MarkerFaceColor',PBoC_blue);
p = plot(0,0);
box on
xlim([.95*burst_dur_range_inMin(1)  burst_dur_range_inMin(end)+burst_dur_range_inMin(1)-.95*burst_dur_range_inMin(1)]);
ylim([0 2 ])    %hard-coded to match amp and dur graphs
yticks([0 1 2])
xlabel('burst duration (min)')
ylabel('Dorsal surge magnitude (au)')
% legend([s_2min s_2minZeroed],'2min surge dur','2min surge dur, zeroed','Location','northwest')
set(gca,'FontSize',14)
StandardFigure(p,gca)
set(gca,'Color',[228 220 209]/255)

saveas(burstDur_vs_surgeMag_zeroed_fig, [FigPath 'burstDur_vs_surgeDur_2minZeroed' '.pdf'])
saveas(burstDur_vs_surgeMag_zeroed_fig, [FigPath 'burstDur_vs_surgeDur_2minZeroed' '.tif'])


burstDur_vs_surgeMagActBurstDur_fig = figure(9);
hold on
% Actual burst duration
e_ActBurstDur = errorbar(burst_dur_range_inMin,dur_surgeMagActBurstDur_mean,dur_surgeMagActBurstDur_ste,'Color','black','LineWidth',1.5);
e_ActBurstDur.CapSize = 0;
s_ActBurstDur = scatter(burst_dur_range_inMin,dur_surgeMagActBurstDur_mean,75,'filled','MarkerFaceColor',PBoC_purple);
p = plot(0,0);
box on
xlim([.95*burst_dur_range_inMin(1)  burst_dur_range_inMin(end)+burst_dur_range_inMin(1)-.95*burst_dur_range_inMin(1)]);
ylim([0 0.2])     %hard-coded so amp and dur graphs match
yticks([0 0.1 0.2])
xlabel('burst duration (min)')
ylabel('Dorsal surge magnitude (au/min)')
legend(s_ActBurstDur,'Total, actual burst dur','Location','northwest')
set(gca,'FontSize',14)
StandardFigure(p,gca)
set(gca,'Color',[228 220 209]/255)

saveas(burstDur_vs_surgeMagActBurstDur_fig, [FigPath 'burstDur_vs_surgeMagActBurstDur' '.pdf'])
saveas(burstDur_vs_surgeMagActBurstDur_fig, [FigPath 'burstDur_vs_surgeMagActBurstDur' '.tif'])


burstDur_vs_surgeRate_fig = figure(10);
hold on
% 20s
e_20s = errorbar(burst_dur_range_inMin,dur_surgeMagRate20s_mean,dur_surgeMagRate20s_ste,'Color','black','LineWidth',1.5);
e_20s.CapSize = 0;
s_20s = scatter(burst_dur_range_inMin,dur_surgeMagRate20s_mean,75,'MarkerFaceColor',PBoC_red,'MarkerEdgeColor','black');
% 40s
e_40s = errorbar(burst_dur_range_inMin,dur_surgeMagRate40s_mean,dur_surgeMagRate40s_ste,'Color','black','LineWidth',1.5);
e_40s.CapSize = 0;
s_40s = scatter(burst_dur_range_inMin,dur_surgeMagRate40s_mean,75,'MarkerFaceColor',PBoC_blue,'MarkerEdgeColor','black');
% 40s
e_60s = errorbar(burst_dur_range_inMin,dur_surgeMagRate60s_mean,dur_surgeMagRate60s_ste,'Color','black','LineWidth',1.5);
e_60s.CapSize = 0;
s_60s = scatter(burst_dur_range_inMin,dur_surgeMagRate60s_mean,75,'MarkerFaceColor',PBoC_green,'MarkerEdgeColor','black');
p = plot(0,0);
box on
xlim([.95*burst_dur_range_inMin(1)  burst_dur_range_inMin(end)+burst_dur_range_inMin(1)-.95*burst_dur_range_inMin(1)]);
ylim([-0.01 0.06])      %hard-coded amp and dur graphs match
yticks([0 0.03 0.06])
xlabel('burst duration (min)')
ylabel('Dorsal surge rate (au/min)')
legend([s_20s s_40s s_60s],'first 20 sec', 'first 40 sec','first 60 sec','Location','northwest')
set(gca,'FontSize',14)
StandardFigure(p,gca)
set(gca,'Color',[228 220 209]/255)

saveas(burstDur_vs_surgeRate_fig, [FigPath 'burstDur_vs_surgeRate' '.pdf'])
saveas(burstDur_vs_surgeRate_fig, [FigPath 'burstDur_vs_surgeRate' '.tif'])

% for i = 1:numel(burst_dur_range)
%    burst_vec = burst_dur_range(i)-burst_dur_sigma:burst_dur_range(i)+burst_dur_sigma;
%    burst_filter = feature_sign_vec == 1 & ismember(lag_dur_vec,burst_vec) & ...
%        lead_dur_vec>= min_buffer_len & lead_dur_vec < max_buffer_len;
%    burst_indices = find(burst_filter);
%    for n = 1:n_boots
%         boot_burst_indices = randsample(burst_dur_indices_vec,numel(burst_dur_indices_vec),true);
%
%         % Calculate averages ...
%         % for sample data
%         burstDur_hmm_boots(i,:,n) = nanmean(hmm_array(boot_burst_indices,:));
%         burstDur_protein_boots(i,:,n) = nanmean(locus_protein_array(boot_burst_indices,:));
%         % for interpolated data
%         protein_interp = interp1(time_vec_data,burstDur_protein_boots(i,:,n),time_vec_interp,'pchip'); % different methods you can use: 'linear', 'pchip', 'makima', 'spline', etc.
%         burstDur_protein_boots_interp(i,:,n) = protein_interp;
% 
%         %Calculate surge duration & magnitude
%         burstDur_protein_trim = burstDur_protein_boots_interp(i,xlim_lb_interp:xlim_ub_interp,n);
%         burstDur_protein_zero = burstDur_protein_trim.*(burstDur_protein_trim >= 0);  %Taking everything above zero (the average trendline at the locus)
%         % Might want to trim already here (or not?)
%         protein_cum_dist_100p(i,n) = trapz(time_vec_interp_trim,burstDur_protein_zero);
%         protein_max_index(i,n) = find(burstDur_protein_zero == max(burstDur_protein_zero));
%         for j = 1:round(numel(burstDur_protein_zero)/2)-1
%             startIndex = protein_max_index(i,n) - j;
%             endIndex = protein_max_index(i,n) + j;
%             if startIndex <= 0 || endIndex > numel(burstDur_protein_zero)
%                 warning('Percentage is too large given the location of the maximum. Use a lower %.')
%                 protein_surge_dur(i,n) = NaN;
%                 break
%             end
%             protein_cum_dist_68p(i,n) = trapz(time_vec_interp_trim(startIndex:endIndex),burstDur_protein_zero(startIndex:endIndex));
%             frac_cum_dist = protein_cum_dist_68p(i,n) / protein_cum_dist_100p(i,n);
%             %Exit for loop once you've hit the desired percent.
%             if frac_cum_dist >= percent_delivered
%                 protein_surge_dur(i,n) = (time_vec_interp_trim(endIndex) - time_vec_interp_trim(startIndex));
%                 break
%             end      
%         end


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