% Figure to build io silencing dataset to compare with stochastic
% simulations
clear
close all
clc

addpath(genpath('./lib'))
addpath(genpath('../../utilities'))

%% Load data

projectName = 'optokni_eve4+6_ON'; 

liveProject = LiveEnrichmentProject(projectName);
resultsRoot = [liveProject.dataPath filesep];

% load data
load([resultsRoot 'spot_struct.mat'])
% load([resultsRoot 'spot_struct_protein.mat'])

%% Filter for only nuclei in correct region with "significant" activity
min_dp = 10;
time_bounds = [5 30]; % nuclei must have been around for full extent of this interval 
ap_bounds = [-0.025 0.025];
sets_to_use = [5];
cal_factor  = 1.35;
switch_time = 940/60;

% create filtered dataset
n_frames_active_vec = [spot_struct.N];
n_filter = n_frames_active_vec >= min_dp;

set_vec = [spot_struct.setID];
set_filter = ismember(set_vec,sets_to_use);

particle_id_vec = [spot_struct.particleID];

ap_pos_mean_vec = NaN(size(n_filter));
for i = 1:length(spot_struct)
    time_vec = spot_struct(i).time/60;
    knirps_vec = spot_struct(i).rawNCProtein;
    if nanmin(time_vec)<=time_bounds(1) && nanmax(time_vec)>=time_bounds(2) && mean(isnan(knirps_vec)) < 0.25
        ap_pos_mean_vec(i) = nanmean(spot_struct(i).APPosNucleus);
    end
end

ap_filter = ap_pos_mean_vec<=ap_bounds(2) & ap_pos_mean_vec>=ap_bounds(1);
spot_struct_trunc = spot_struct(n_filter & ap_filter & set_filter);
particle_id_filt = particle_id_vec(n_filter & ap_filter & set_filter);
% now generate arrays containing fraction of active nuclei and knirps as a
% function of time
time_index = unique([spot_struct_trunc.timeInterp])/60;
time_ref = time_index(time_index>=time_bounds(1)&time_index<=time_bounds(2));
index_vec = 1:length(time_ref);

io_ref_struct = struct;
io_ref_struct.on_off_array = false(length(time_ref),length(spot_struct_trunc));
io_ref_struct.fluo_array = zeros(length(time_ref),length(spot_struct_trunc));
io_ref_struct.knirps_array = NaN(length(time_ref),length(spot_struct_trunc));

for i = 1:length(spot_struct_trunc)
    % extract vectors
    time_vec_interp = spot_struct_trunc(i).timeInterp/60;
    time_vec_raw = spot_struct_trunc(i).time/60;
    fluo_vec_interp = spot_struct_trunc(i).fluoInterp;
    knirps_vec = spot_struct_trunc(i).rawNCProtein;
    knirps_vec_interp = interp1(time_vec_raw,knirps_vec,time_ref);
    % perform calculations
    frac_on_vec = false(length(time_ref),1);
    match_to = ismember(time_ref,time_vec_interp);
    match_from = ismember(time_vec_interp,time_ref);
    
    io_ref_struct.fluo_array(match_to,i) = fluo_vec_interp(match_from);
    io_ref_struct.on_off_array(:,i) = io_ref_struct.fluo_array(:,i) > 0;
    io_ref_struct.knirps_array(:,i) = knirps_vec_interp;
    
end    

% remove NaNs
nn_filter = ~any(isnan(io_ref_struct.knirps_array));

io_ref_struct.particle_id_vec = particle_id_filt(nn_filter);
io_ref_struct.knirps_array = io_ref_struct.knirps_array(:,nn_filter);
io_ref_struct.fluo_array = io_ref_struct.fluo_array(:,nn_filter);
io_ref_struct.on_off_array = io_ref_struct.on_off_array(:,nn_filter);

% adjust knirps array (kind of ad-hoc atm)
mean_knirps_vec = nanmean(io_ref_struct.knirps_array,2);
switch_i = find(time_ref==switch_time);
k_before = mean_knirps_vec(switch_i);
% k_after = mean_knirps_vec(switch_i+1);
% cv_factor = k_after/k_before;
io_ref_struct.knirps_array_norm = io_ref_struct.knirps_array;
io_ref_struct.knirps_array_norm(switch_i+1:end,:) = io_ref_struct.knirps_array_norm(switch_i+1:end,:)/cal_factor;

% add some basic metadata
io_ref_struct.projectName = projectName;
io_ref_struct.deltaT = spot_struct_trunc(1).tresInterp;
io_ref_struct.min_dp = min_dp;
io_ref_struct.ap_bounds = ap_bounds;
io_ref_struct.time_bounds = time_bounds;
io_ref_struct.time_vec = time_ref;
io_ref_struct.set_ids_used = sets_to_use;

% save
save([resultsRoot 'io_ref_struct.mat'],'io_ref_struct')

