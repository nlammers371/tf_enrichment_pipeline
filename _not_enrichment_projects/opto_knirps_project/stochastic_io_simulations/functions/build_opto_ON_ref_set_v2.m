% Script to build io silencing dataset to compare with stochastic
% simulations
function io_ref_ON = build_opto_ON_ref_set_v2(projectNameCONST)

addpath(genpath('./lib'))
addpath(genpath('../../utilities'))

% Load data
% projectNameWT = 'optokni_eve4+6_WT_FUN'; 

% liveProject = LiveEnrichmentProject(projectNameWT);
% resultsRoot = [liveProject.dataPath filesep];
dataRoot = 'S:\Nick\Dropbox\ProcessedEnrichmentData\';
resultsRoot = [dataRoot projectNameCONST filesep];
% load data
load([resultsRoot 'spot_struct.mat'],'spot_struct')

% indicate frames where opto light was turned ON
setID_vec = [1 2 3];
on_frame_vec = [65 52 43];

% Filter for only nuclei in correct region with "significant" activity
min_dp = 10;
time_bounds = [8 30];%[7 30]; % nuclei must have been around for full extent of this interval 
ap_bounds = [-0.09 0.07];
knirps_offset = 375000 / 1e5;
cal_slope = 1.243;
cal_intercept = 1.079e5 / 1e5; %NL: dividing everything through by 1e5 for simplicity

TresInterp = 20;

% define time axis 
time_index = 0:TresInterp:40*60;

% initialize data structure to store results
io_ref_ON = struct;

% Initialize vectors to capture key average transcriptional behaviors
io_ref_ON.off_time_vec = NaN(size(spot_struct));
io_ref_ON.on_time_vec = NaN(size(spot_struct));
io_ref_ON.off_spot_fluo = NaN(size(spot_struct));
io_ref_ON.off_knirps_vec = NaN(size(spot_struct)); % knirps level at time of shutoff
io_ref_ON.mean_ap = NaN(size(spot_struct)); % average AP position (over active liftime)
io_ref_ON.mean_fluo = NaN(size(spot_struct)); % average spot intensity (over active liftime)
io_ref_ON.mean_fluo_n = NaN(size(spot_struct));

% initialize arrays to store full trace info
io_ref_ON.knirps_array = NaN(length(time_index),length(spot_struct));
io_ref_ON.fluo_array = NaN(length(time_index),length(spot_struct));
io_ref_ON.fluo_full_array = zeros(length(time_index),length(spot_struct));
io_ref_ON.fluo_offset_array = NaN(length(time_index),length(spot_struct));
io_ref_ON.off_flag_array = NaN(length(time_index),length(spot_struct));

% track which entries to keep
keep_flags = false(1,length(spot_struct));

% iterate through all struct entries
for i = 1:length(spot_struct)
  
    % extract core vectors 
    fluo_vec = spot_struct(i).fluo;
    offset_vec = spot_struct(i).fluoOffset;
    time_vec = spot_struct(i).time;
    knirps_vec = spot_struct(i).rawNCProtein / 1e5;
    ap_vec = spot_struct(i).APPosNucleus;
    setID = spot_struct(i).setID;
    frame_vec = spot_struct(i).frames;
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
        io_ref_ON.mean_ap(i) = nanmean(ap_vec(start_i:stop_i));
        io_ref_ON.mean_knirps(i) = nanmean(knirps_vec(start_i:stop_i));

        % shut-off info
        io_ref_ON.off_time_vec(i) = min([time_vec(stop_i), 40*60]);                
        post_off_vec(stop_i+1:end) = 1;

        if start_i > 1
            io_ref_ON.on_time_vec(i) = time_vec(start_i);            
        end
        post_on_vec(start_i+1:end) = 1;
        [~, on_i] = min(abs(time_index-time_vec(start_i)));
        [~, off_i] = min(abs(time_index-io_ref_ON.off_time_vec(i)));
        io_ref_ON.off_flag_array(on_i:end,i) = 1;
        io_ref_ON.off_flag_array(off_i+1:end,i) = 0;
        % make regression vectors              
        
        % make vectors where missing data points are replaced with zeros
        % only between first and last detections
        fluo_zeros = fluo_vec;        
        fluo_zeros(post_on_vec&~post_off_vec&isnan(fluo_vec)) = 0;
        time_vec_active = time_vec(post_on_vec&~post_off_vec);
        
        % interpolate
        time_interp = time_index(time_index>=time_vec_active(1) & time_index<=time_vec_active(end));
        fluo_zeros_interp = interp1(time_vec_active,fluo_zeros(post_on_vec&~post_off_vec),time_interp,'nearest');
        offset_interp = interp1(time_vec,offset_vec,time_index,'nearest');
        
        % store
        io_ref_ON.fluo_offset_array(:,i) = offset_interp;
        io_ref_ON.fluo_full_array(ismember(time_index,time_interp),i) = fluo_zeros_interp;
        io_ref_ON.fluo_array(ismember(time_index,time_interp),i) = fluo_zeros_interp;
        
        % store average
        io_ref_ON.mean_fluo(i) = nanmean(fluo_zeros_interp);
        io_ref_ON.mean_fluo_n(i) = sum(~isnan(fluo_zeros_interp));
        
        % interpolate the knirps info
        kni_nan_ft = ~isnan(knirps_vec);
        time_interp = time_index(time_index>=min(time_vec(kni_nan_ft)) & time_index<=max(time_vec(kni_nan_ft)));
        knirps_interp = interp1(time_vec(kni_nan_ft),knirps_vec(kni_nan_ft),time_interp,'linear');
        
        % record
        kni_temp = io_ref_ON.knirps_array(:,i);
        kni_temp(ismember(time_interp,time_index)) = knirps_interp;
        kni_temp_filled = kni_temp;
        si = find(~isnan(kni_temp_filled),1);
        fi = find(~isnan(kni_temp_filled),1,'last');
        kni_temp_filled(1:si-1) = kni_temp_filled(si);
        kni_temp_filled(fi+1:end) = kni_temp_filled(fi); 
        % apply laser correction
        % Account for added laser intensity
        pert_ind = find(frame_vec==on_frame_vec(setID_vec==setID),1);
        kni_temp_norm = kni_temp_filled;
        kni_temp_norm(pert_ind:end) = ...
                (kni_temp_norm(pert_ind:end)-cal_intercept)/cal_slope;
        kni_temp_norm = kni_temp_norm - knirps_offset;
        io_ref_ON.knirps_array(:,i) = kni_temp_norm;
        
        if ~any(isnan(kni_temp_filled)) && io_ref_ON.mean_ap(i)>=ap_bounds(1) && io_ref_ON.mean_ap(i)<=ap_bounds(2)
            keep_flags(i) = 1;
        end
    end
end

% Calculate AP trends for comparison with simulations
% first let's remove unwanted entries
io_ref_ON.knirps_array = io_ref_ON.knirps_array(:,keep_flags);
io_ref_ON.knirps_array(io_ref_ON.knirps_array<0) = 1e-6; % setting this to small nonzero value saves pain later
io_ref_ON.fluo_array = io_ref_ON.fluo_array(:,keep_flags);
io_ref_ON.fluo_full_array = io_ref_ON.fluo_full_array(:,keep_flags);
io_ref_ON.off_flag_array = io_ref_ON.off_flag_array(:,keep_flags);
io_ref_ON.mean_fluo = io_ref_ON.mean_fluo(keep_flags);
io_ref_ON.mean_fluo_n = io_ref_ON.mean_fluo_n(keep_flags);
io_ref_ON.on_time_vec = io_ref_ON.on_time_vec(keep_flags);
io_ref_ON.off_time_vec = io_ref_ON.off_time_vec(keep_flags);
io_ref_ON.mean_knirps = io_ref_ON.mean_knirps(keep_flags);
io_ref_ON.mean_ap = io_ref_ON.mean_ap(keep_flags);

% now calculate mean trends vs AP
ap_bins = linspace(ap_bounds(1),ap_bounds(2),13); % used in mean figure
ap_groups_mean = discretize(io_ref_ON.mean_ap,ap_bins); 

nBoots = 100;
off_time_vec_mean = NaN(nBoots,length(ap_bins)-1);
fluo_vec_mean = NaN(nBoots,length(ap_bins)-1);
knirps_vec_mean = NaN(nBoots,length(ap_bins)-1);

for a = 1:length(ap_bins)-1          
    ap_indices = find(ap_groups_mean==a);
    for n = 1:nBoots
        boot_indices = randsample(ap_indices,length(ap_indices),true);
        off_time_vec_mean(n,a) = mean(io_ref_ON.off_time_vec(boot_indices));
        fluo_z_ft = io_ref_ON.fluo_array(:,boot_indices);
        fluo_vec_mean(n,a) = nanmean(fluo_z_ft(:));%.*io_ref_wt.mean_fluo_n(boot_indices))/sum(io_ref_wt.mean_fluo_n(boot_indices));
        knirps_vec_mean(n,a) = mean(io_ref_ON.mean_knirps(boot_indices));                         
    end
end

% store
io_ref_ON.off_time_vec_mean = mean(off_time_vec_mean);
io_ref_ON.off_time_vec_ste = std(off_time_vec_mean);
io_ref_ON.fluo_vec_mean = mean(fluo_vec_mean);
io_ref_ON.fluo_vec_ste = std(fluo_vec_mean);
io_ref_ON.knirps_vec_mean = knirps_vec_mean;
io_ref_ON.ap_axis_mean = ap_bins(1:end-1) + diff(ap_bins)/2;
io_ref_ON.time_axis = time_index;

% Generate fraction off vs knirps curve
knirps_vec_long = io_ref_ON.knirps_array(:);

% off_flags_long = io_ref_wt.off_flag_array(:);
knirps_bins = logspace(log10(prctile(knirps_vec_long,1)),log10(prctile(knirps_vec_long,99)),26);
knirps_axis = knirps_bins(1:end-1) + diff(knirps_bins)/2;
% knirps_groups = discretize(knirps_vec_long, knirps_bins); 

% 
% io_ref_wt.fraction_on_vec = fraction_on_vec;
io_ref_ON.knirps_axis_fon = knirps_axis;
io_ref_ON.knirps_bins_fon = knirps_bins;

%% Estimate p-still-on titration curve
ap_limit = [-.02 .02];
ap_filter = io_ref_ON.mean_ap>=ap_limit(1) & io_ref_ON.mean_ap<=ap_limit(2);

knirps_array = io_ref_ON.knirps_array(:,ap_filter);
knirps_vec_long = knirps_array(:);

off_flag_array = io_ref_ON.off_flag_array(:,ap_filter);
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

io_ref_ON.knirps_axis_still_on = knirps_axis;
io_ref_ON.knirps_bins_still_on = knirps_index;
io_ref_ON.fraction_still_on_mean = nanmean(fraction_still_on_array,2);
io_ref_ON.fraction_still_on_ste = nanstd(fraction_still_on_array,[],2)+1e-2;
io_ref_ON.ap_limits_still_on = ap_limit;

%% estimate mean fluorescence over time amongst active traces
boot_options = find(ap_filter);
mean_fluo_time = NaN(length(time_index),nBoots);
mean_knirps_time = NaN(length(time_index),nBoots);
mean_offset_time = NaN(length(time_index),nBoots);
for n = 1:nBoots
    boot_indices = randsample(boot_options,length(boot_options),true);
    fluo_temp = io_ref_ON.fluo_array(:,boot_indices);
    knirps_temp = io_ref_ON.knirps_array(:,boot_indices);
    offset_temp = io_ref_ON.fluo_offset_array(:,boot_indices);
    knirps_temp(isnan(fluo_temp)) = NaN;
    
    mean_fluo_time(:,n) = nanmean(fluo_temp,2);
    mean_knirps_time(:,n) = nanmean(knirps_temp,2);
    mean_offset_time(:,n) = nanmean(offset_temp,2);
end

% filter 
time_filter = time_index/60>=time_bounds(1)&time_index/60<=time_bounds(2);
io_ref_ON.time_axis_mf = time_index(time_filter)/60;

% calculate mean and standard error
t_counts = sum(~isnan(io_ref_ON.fluo_array(:,ap_filter)),2);
mean_fluo_time(t_counts<10,:) = NaN;
io_ref_ON.fluo_time_mean = nanmean(mean_fluo_time(time_filter,:),2);
io_ref_ON.knirps_time_mean = nanmean(mean_knirps_time(time_filter,:),2);
io_ref_ON.offset_time_mean = nanmean(mean_offset_time(time_filter,:),2);
io_ref_ON.fluo_time_ste = nanstd(mean_fluo_time(time_filter,:),[],2);

