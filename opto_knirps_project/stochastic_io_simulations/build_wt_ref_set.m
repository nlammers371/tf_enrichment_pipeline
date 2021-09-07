% Script to build io silencing dataset to compare with stochastic
% simulations
function io_ref_wt = build_wt_ref_set(projectNameWT)
% clear
% close all
% clc

addpath(genpath('./lib'))
addpath(genpath('../../utilities'))

% Load data
% projectNameWT = 'optokni_eve4+6_WT_FUN'; 

liveProject = LiveEnrichmentProject(projectNameWT);
resultsRoot = [liveProject.dataPath filesep];

% load data
load([resultsRoot 'spot_struct.mat'],'spot_struct')

% Filter for only nuclei in correct region with "significant" activity
min_dp = 10;
time_bounds = [8 25];%[7 30]; % nuclei must have been around for full extent of this interval 
ap_bounds = [-0.09 0.07];
knirps_offset = 375000 / 1e5;
TresInterp = 20;

% define time axis 
time_index = 0:TresInterp:40*60;

% initialize data structure to store results
io_ref_wt = struct;

% Initialize vectors to capture key average transcriptional behaviors
io_ref_wt.off_time_vec = NaN(size(spot_struct));
io_ref_wt.on_time_vec = NaN(size(spot_struct));
io_ref_wt.off_spot_fluo = NaN(size(spot_struct));
io_ref_wt.off_knirps_vec = NaN(size(spot_struct)); % knirps level at time of shutoff
io_ref_wt.mean_ap = NaN(size(spot_struct)); % average AP position (over active liftime)
io_ref_wt.mean_fluo = NaN(size(spot_struct)); % average spot intensity (over active liftime)
io_ref_wt.mean_fluo_n = NaN(size(spot_struct));

% initialize arrays to store full trace info
io_ref_wt.knirps_array = NaN(length(time_index),length(spot_struct));
io_ref_wt.fluo_array = NaN(length(time_index),length(spot_struct));
io_ref_wt.fluo_full_array = zeros(length(time_index),length(spot_struct));
io_ref_wt.off_flag_array = NaN(length(time_index),length(spot_struct));

% track which entries to keep
keep_flags = false(1,length(spot_struct));

% iterate through all struct entries
for i = 1:length(spot_struct)
  
    % extract core vectors 
    fluo_vec = spot_struct(i).fluo;
    time_vec = spot_struct(i).time;
    knirps_vec = spot_struct(i).rawNCProtein / 1e5;
    ap_vec = spot_struct(i).APPosNucleus;
    
    % find first active time
    first_time = time_vec(find(~isnan(fluo_vec),1))/60;
    if isempty(first_time)
        first_time = Inf;
    end
    % filter for nuclei that are present and that had activity
    if time_vec(end)/60 >= time_bounds(2) && time_vec(1) <=time_bounds(1)/60 && first_time < 15 &&...
          sum(~isnan(fluo_vec))>=min_dp && mean(~isnan(knirps_vec)) > 0.9
      
        % get off and on indices        
        post_on_vec = zeros(size(ap_vec));
        post_off_vec = zeros(size(ap_vec));
        
        % get start and stop indices
        start_i = find(~isnan(fluo_vec),1);
        stop_i = find(~isnan(fluo_vec),1,'last');
        
        % save lifetime values
        io_ref_wt.mean_ap(i) = nanmean(ap_vec(start_i:stop_i));
        io_ref_wt.mean_knirps(i) = nanmean(knirps_vec(start_i:stop_i));

        % shut-off info
        io_ref_wt.off_time_vec(i) = min([time_vec(stop_i), 40*60]);                
        post_off_vec(stop_i+1:end) = 1;

        if start_i > 1
            io_ref_wt.on_time_vec(i) = time_vec(start_i);            
        end
        post_on_vec(start_i+1:end) = 1;
        [~, on_i] = min(abs(time_index-time_vec(start_i)));
        [~, off_i] = min(abs(time_index-io_ref_wt.off_time_vec(i)));
        io_ref_wt.off_flag_array(on_i:end,i) = 1;
        io_ref_wt.off_flag_array(off_i+1:end,i) = 0;
        % make regression vectors              
        
        % make vectors where missing data points are replaced with zeros
        % only between first and last detections
        fluo_zeros = fluo_vec;        
        fluo_zeros(post_on_vec&~post_off_vec&isnan(fluo_vec)) = 0;
        time_vec_active = time_vec(post_on_vec&~post_off_vec);
        
        % interpolate
        time_interp = time_index(time_index>=time_vec_active(1) & time_index<=time_vec_active(end));
        fluo_zeros_interp = interp1(time_vec_active,fluo_zeros(post_on_vec&~post_off_vec),time_interp,'nearest');
        
        % store
        io_ref_wt.fluo_full_array(ismember(time_index,time_interp),i) = fluo_zeros_interp;
        io_ref_wt.fluo_array(ismember(time_index,time_interp),i) = fluo_zeros_interp;
        
        % store average
        io_ref_wt.mean_fluo(i) = nanmean(fluo_zeros_interp);
        io_ref_wt.mean_fluo_n(i) = sum(~isnan(fluo_zeros_interp));
        
        % interpolate the knirps info
        kni_nan_ft = ~isnan(knirps_vec);
        time_interp = time_index(time_index>=min(time_vec(kni_nan_ft)) & time_index<=max(time_vec(kni_nan_ft)));
        knirps_interp = interp1(time_vec(kni_nan_ft),knirps_vec(kni_nan_ft),time_interp,'linear');
        
        % record
        kni_temp = io_ref_wt.knirps_array(:,i);
        kni_temp(ismember(time_interp,time_index)) = knirps_interp;
        kni_temp_filled = kni_temp;
        si = find(~isnan(kni_temp_filled),1);
        fi = find(~isnan(kni_temp_filled),1,'last');
        kni_temp_filled(1:si-1) = kni_temp_filled(si);
        kni_temp_filled(fi+1:end) = kni_temp_filled(fi); 
        io_ref_wt.knirps_array(:,i) = kni_temp_filled - knirps_offset;

        if ~any(isnan(kni_temp_filled)) && io_ref_wt.mean_ap(i)>=ap_bounds(1) && io_ref_wt.mean_ap(i)<=ap_bounds(2)
            keep_flags(i) = 1;
        end
    end
end

% Calculate AP trends for comparison with simulations
% first let's remove unwanted entries
io_ref_wt.knirps_array = io_ref_wt.knirps_array(:,keep_flags);
io_ref_wt.knirps_array(io_ref_wt.knirps_array<0) = 1e-6; % setting this to small nonzero value saves pain later
io_ref_wt.fluo_array = io_ref_wt.fluo_array(:,keep_flags);
io_ref_wt.fluo_full_array = io_ref_wt.fluo_full_array(:,keep_flags);
io_ref_wt.off_flag_array = io_ref_wt.off_flag_array(:,keep_flags);
io_ref_wt.mean_fluo = io_ref_wt.mean_fluo(keep_flags);
io_ref_wt.mean_fluo_n = io_ref_wt.mean_fluo_n(keep_flags);
io_ref_wt.on_time_vec = io_ref_wt.on_time_vec(keep_flags);
io_ref_wt.off_time_vec = io_ref_wt.off_time_vec(keep_flags);
io_ref_wt.mean_knirps = io_ref_wt.mean_knirps(keep_flags);
io_ref_wt.mean_ap = io_ref_wt.mean_ap(keep_flags);

% keep_flags = true(1,sum(keep_flags));

% now calculate mean trends vs AP
ap_bins = linspace(ap_bounds(1),ap_bounds(2),13); % used in mean figure
ap_groups_mean = discretize(io_ref_wt.mean_ap,ap_bins); 

nBoots = 100;
off_time_vec_mean = NaN(nBoots,length(ap_bins)-1);
fluo_vec_mean = NaN(nBoots,length(ap_bins)-1);
knirps_vec_mean = NaN(nBoots,length(ap_bins)-1);

for a = 1:length(ap_bins)-1          
    ap_indices = find(ap_groups_mean==a);
    for n = 1:nBoots
        boot_indices = randsample(ap_indices,length(ap_indices),true);
        off_time_vec_mean(n,a) = mean(io_ref_wt.off_time_vec(boot_indices));
        fluo_z_ft = io_ref_wt.fluo_array(:,boot_indices);
        fluo_vec_mean(n,a) = nanmean(fluo_z_ft(:));%.*io_ref_wt.mean_fluo_n(boot_indices))/sum(io_ref_wt.mean_fluo_n(boot_indices));
        knirps_vec_mean(n,a) = mean(io_ref_wt.mean_knirps(boot_indices));                         
    end
end

% store
io_ref_wt.off_time_vec_mean = mean(off_time_vec_mean);
io_ref_wt.off_time_vec_ste = std(off_time_vec_mean);
io_ref_wt.fluo_vec_mean = mean(fluo_vec_mean);
io_ref_wt.fluo_vec_ste = std(fluo_vec_mean);
io_ref_wt.knirps_vec_mean = knirps_vec_mean;
io_ref_wt.ap_axis_mean = ap_bins(1:end-1) + diff(ap_bins)/2;
io_ref_wt.time_axis = time_index;

% Generate fraction off vs knirps curve
knirps_vec_long = io_ref_wt.knirps_array(:);

% off_flags_long = io_ref_wt.off_flag_array(:);
knirps_bins = logspace(log10(prctile(knirps_vec_long,1)),log10(prctile(knirps_vec_long,99)),26);
knirps_axis = knirps_bins(1:end-1) + diff(knirps_bins)/2;
% knirps_groups = discretize(knirps_vec_long, knirps_bins); 

% fraction_on_vec = NaN(1,length(knirps_bins)-1);
% 
% for k = 1:length(knirps_bins)-1
%     fraction_on_vec(k) = mean(~off_flags_long(knirps_groups==k));
% end
% 
% io_ref_wt.fraction_on_vec = fraction_on_vec;
io_ref_wt.knirps_axis_fon = knirps_axis;
io_ref_wt.knirps_bins_fon = knirps_bins;

%% Estimate our detection threshold
% get empirical distribution of min fluo values
min_fluo_vec = NaN(1,length(spot_struct));
for i = 1:length(spot_struct)
    min_fluo_vec(i) = nanmin(spot_struct(i).fluo);
end
bins = linspace(0,5e5);
counts = histcounts(min_fluo_vec,bins);
% counts = counts/sum(counts);
bins_fit = bins(1:end-1) + diff(bins)/2;

% fit a simple Gaussian function
gauss_fun = @(x) x(3)*exp(-0.5*((x(1)-bins_fit)./x(2)).^2);
ob_fun = @(x) counts-gauss_fun(x);
fit_parameters = lsqnonlin(ob_fun,[1e4 1e4 100],[0 0 0],[Inf Inf Inf]);

% record
io_ref_wt.min_fluo_bins = bins_fit;
io_ref_wt.min_fluo_counts = counts;
io_ref_wt.fit_parameters = fit_parameters;
io_ref_wt.gauss_fit = gauss_fun(fit_parameters);
io_ref_wt.min_spot_values = min_fluo_vec;
io_ref_wt.F_min_fit = fit_parameters(1);
io_ref_wt.F_min_std_fit = fit_parameters(2);

%% Estimate p-still-on titration curve
ap_limit = [-.02 .02];
ap_filter = io_ref_wt.mean_ap>=ap_limit(1) & io_ref_wt.mean_ap<=ap_limit(2);

knirps_array = io_ref_wt.knirps_array(:,ap_filter);
knirps_vec_long = knirps_array(:);

off_flag_array = io_ref_wt.off_flag_array(:,ap_filter);
off_vec_long = off_flag_array(:);

knirps_index = linspace(nanmin(knirps_vec_long),nanmax(knirps_vec_long),51);

knirps_axis = knirps_index(1:end-1) + diff(knirps_index)/2;
knirps_groups = discretize(knirps_vec_long,knirps_index);
fraction_still_on_array = NaN(length(knirps_axis),nBoots);
for k = 1:length(knirps_axis)
    k_options = find(knirps_groups==k);
    for n = 1:nBoots
        boot_indices = randsample(k_options,length(k_options),true);
        fraction_still_on_array(k,n) = nanmean(off_vec_long(boot_indices));
    end
end    

io_ref_wt.knirps_axis_still_on = knirps_axis;
io_ref_wt.knirps_bins_still_on = knirps_index;
io_ref_wt.fraction_still_on_mean = nanmean(fraction_still_on_array,2);
io_ref_wt.fraction_still_on_ste = nanstd(fraction_still_on_array,[],2);
io_ref_wt.ap_limits_still_on = ap_limit;

%% estimate mean fluorescence over time amongst active traces
boot_options = find(ap_filter);
mean_fluo_time = NaN(length(time_index),nBoots);
for n = 1:nBoots
    boot_indices = randsample(boot_options,length(boot_options),true);
    mean_fluo_time(:,n) = nanmean(io_ref_wt.fluo_array(:,boot_indices),2);
end
% filter 
time_filter = time_index/60>=time_bounds(1)&time_index/60<=time_bounds(2);
io_ref_wt.time_axis_mf = time_index(time_filter)/60;
% calculate mean and standard error
t_counts = sum(~isnan(io_ref_wt.fluo_array(:,ap_filter)),2);
mean_fluo_time(t_counts<10,:) = NaN;
io_ref_wt.fluo_time_mean = nanmean(mean_fluo_time(time_filter,:),2);
io_ref_wt.fluo_time_ste = nanstd(mean_fluo_time(time_filter,:),[],2);

% now let's do something a little mor sophisticated: try to to estimate
% probabilities of missing a detection as a function of spot fluorescence
io_ref_wt = estimateDetectionThreshold(io_ref_wt,spot_struct);
% %% save
% save([resultsRoot 'io_ref_wt.mat'],'io_ref_wt')
