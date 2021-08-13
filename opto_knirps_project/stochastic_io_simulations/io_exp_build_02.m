% Figure to build io silencing dataset to compare with stochastic
% simulations
clear
close all
% clc

addpath(genpath('./lib'))
addpath(genpath('../../utilities'))

% Load data
projectName = 'optokni_eve4+6_WT_FUN'; 

liveProject = LiveEnrichmentProject(projectName);
resultsRoot = [liveProject.dataPath filesep];

% load data
load([resultsRoot 'spot_struct.mat'])

% Filter for only nuclei in correct region with "significant" activity
min_dp = 10;
time_bounds = [7 30]; % nuclei must have been around for full extent of this interval 
ap_bounds = [-0.07 0.05];
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

% initialize arrays to store full trace info
io_ref_wt.knirps_array = NaN(length(time_index),length(spot_struct));
io_ref_wt.fluo_array = NaN(length(time_index),length(spot_struct));
io_ref_wt.fluo_full_array = zeros(length(time_index),length(spot_struct));
io_ref_wt.off_flag_array = zeros(length(time_index),length(spot_struct));

% initialize longform vectors for regression
ap_vec_long = [];
fluo_zeros_long = [];
post_turn_on_flags = [];
post_turn_off_flags = [];

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
          sum(~isnan(fluo_vec))>=min_dp && mean(~isnan(knirps_vec)) >0.9
      
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
        io_ref_wt.off_time_vec(i) = time_vec(stop_i);                
        post_off_vec(stop_i+1:end) = 1;

        if start_i > 1
            io_ref_wt.on_time_vec(i) = time_vec(start_i);            
        end
        post_on_vec(start_i+1:end) = 1;
                
        [~, off_i] = min(abs(time_index-io_ref_wt.off_time_vec(i)));
        io_ref_wt.off_flag_array(off_i+1:end,i) = 1;
        % make regression vectors       
%         post_turn_on_flags = [post_turn_on_flags post_on_vec];
%         post_turn_off_flags = [post_turn_off_flags post_off_vec];        
        
        % make vectors where missing data points are replaced with zeros
        % only between first and last detections
        fluo_zeros = fluo_vec;        
        fluo_zeros(post_on_vec&~post_off_vec&isnan(fluo_vec)) = 0;
        time_vec_active = time_vec(post_on_vec&~post_off_vec);
        
        % interpolate
        time_interp = time_index(time_index>=time_vec_active(1) & time_index<=time_vec_active(end));
        fluo_zeros_interp = interp1(time_vec_active,fluo_zeros(post_on_vec&~post_off_vec),time_interp,'nearest');
        
        % store
        io_ref_wt.fluo_full_array(ismember(time_interp,time_index),i) = fluo_zeros_interp;
        io_ref_wt.fluo_array(ismember(time_interp,time_index),i) = fluo_zeros_interp;
        
        % store average
        io_ref_wt.mean_fluo(i) = nanmean(fluo_zeros_interp);
        
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
        % longform stuff
%         knirps_vec_long_raw = [knirps_vec_long_raw knirps_vec];
%         ap_vec_long = [ap_vec_long ap_vec];
%         time_vec_long = [time_vec_long time_vec];
        if ~any(isnan(kni_temp_filled)) && io_ref_wt.mean_ap(i)>=ap_bounds(1) && io_ref_wt.mean_ap(i)<=ap_bounds(2)
            keep_flags(i) = 1;
        end
    end
end

% Calculate AP trends for comparison with simulations
% first let's remove unwanted entries
io_ref_wt.knirps_array = io_ref_wt.knirps_array(:,keep_flags);
io_ref_wt.knirps_array(io_ref_wt.knirps_array<0) = 0;
io_ref_wt.fluo_array = io_ref_wt.fluo_array(:,keep_flags);
io_ref_wt.fluo_full_array = io_ref_wt.fluo_full_array(:,keep_flags);
io_ref_wt.off_flag_array = io_ref_wt.off_flag_array(:,keep_flags);
io_ref_wt.mean_fluo = io_ref_wt.mean_fluo(keep_flags);
io_ref_wt.on_time_vec = io_ref_wt.on_time_vec(keep_flags);
io_ref_wt.off_time_vec = io_ref_wt.off_time_vec(keep_flags);
io_ref_wt.mean_knirps = io_ref_wt.mean_knirps(keep_flags);
io_ref_wt.mean_ap = io_ref_wt.mean_ap(keep_flags);

keep_flags = true(1,sum(keep_flags));

% now calculate mean trends
ap_bins = linspace(ap_bounds(1),ap_bounds(2),11); % used in mean figure
ap_groups_mean = discretize(io_ref_wt.mean_ap,ap_bins); 

off_time_vec_mean = NaN(1,length(ap_bins)-1);
fluo_vec_mean = NaN(1,length(ap_bins)-1);
knirps_vec_mean = NaN(1,length(ap_bins)-1);


for a = 1:length(ap_bins)-1          
    ap_filter = ap_groups_mean==a;
    off_time_vec_mean(a) = mean(io_ref_wt.off_time_vec(ap_filter));
    fluo_vec_mean(a) = mean(io_ref_wt.mean_fluo(ap_filter));
    knirps_vec_mean(a) = mean(io_ref_wt.mean_knirps(ap_filter));                         
end

% store
io_ref_wt.off_time_vec_mean = off_time_vec_mean;
io_ref_wt.fluo_vec_mean = fluo_vec_mean;
io_ref_wt.knirps_vec_mean = knirps_vec_mean;
io_ref_wt.ap_axis_mean = ap_bins(1:end-1) + diff(ap_bins)/2;
io_ref_wt.time_axis = time_index;

% Generate fraction off vs knirps curve
knirps_vec_long = io_ref_wt.knirps_array(:);
off_flags_long = io_ref_wt.off_flag_array(:);
knirps_bins = logspace(log10(prctile(knirps_vec_long,1)),log10(prctile(knirps_vec_long,99)),26);
knirps_axis = knirps_bins(1:end-1) + diff(knirps_bins)/2;
knirps_groups = discretize(knirps_vec_long, knirps_bins); 

fraction_on_vec = NaN(1,length(knirps_bins)-1);

for k = 1:length(knirps_bins)-1
    fraction_on_vec(k) = mean(~off_flags_long(knirps_groups==k));
end

io_ref_wt.fraction_on_vec = fraction_on_vec;
io_ref_wt.knirps_axis_fon = knirps_axis;
io_ref_wt.knirps_bins_fon = knirps_bins;

% save
save([resultsRoot 'io_ref_wt.mat'],'io_ref_wt')
