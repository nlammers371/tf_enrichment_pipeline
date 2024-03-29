function io_ref_ra = build_ra_ref_set_v2(projectNameRA, dataRoot)
% Script to build io silencing dataset to compare with stochastic
% simulations
% clear
% close all
% clc

addpath(genpath('./lib'))
addpath(genpath('../../utilities'))

% Load data
% projectNameRA = 'optokni_eve4+6_ON'; 

% liveProject = LiveEnrichmentProject(projectNameRA);
% dataRoot = 'S:\Nick\Dropbox (Personal)\ProcessedEnrichmentData\';
resultsRoot = [dataRoot projectNameRA filesep];

% load data
load([resultsRoot 'spot_struct.mat'],'spot_struct')

% Filter for only nuclei in correct region with "significant" activity
knirps_offset = 3.75e5 / 1e5;
min_dp = 10;
window_size = 40; % size of lead and lag interval, in frames
time_bounds = [8 25]; % nuclei must have been around for full extent of this interval 
ap_bounds = [-0.01 0.01];
sets_to_use = [2 3 5 6];
on_frames = [36 43 45 43];
cal_slope = 1.243;
cal_intercept = 1.079e5 / 1e5; %NL: dividng everything through by 1e5 for simplicity
off_frame_ref = -1; % let's say it needs to have been off for at least two minutes prior

% iterate through sets one-by-one and extract relevant io data
master_struct = struct;

% generate key reference vectors
set_vec = [spot_struct.setID];
n_frames_active_vec = [spot_struct.N];
n_filter = n_frames_active_vec >= min_dp;
time_ref = (-window_size:window_size)*spot_struct(1).tresInterp;
index_vec = -window_size:window_size;
full_denom_count = 0;    
for s = 1:length(sets_to_use)
  
    % extract key indentifying parameters
    setID = sets_to_use(s);
    switch_frame = on_frames(s);

    % create filtered dataset            
    set_filter = set_vec==setID;
    particle_id_vec = [spot_struct.particleID];
    ap_pos_mean_vec = NaN(size(n_filter));
    first_fluo_time_vec = NaN(size(n_filter));
    first_nuc_time_vec = NaN(size(n_filter));
    last_nuc_time_vec = NaN(size(n_filter));
    
    % geneate set-specific time ref vectors    
    time_index_interp = unique([spot_struct(set_filter).timeInterp])/60;
    time_index_interp = time_index_interp(~isnan(time_index_interp));
    time_index_raw = unique([spot_struct(set_filter).time])/60;
    time_index_raw = time_index_raw(~isnan(time_index_raw));
    frame_index_raw = unique([spot_struct(set_filter).frames]);
    switch_time = time_index_raw(frame_index_raw==switch_frame); % actual time of change
    
    % get set-specific index
    switch_i = find(time_index_interp<switch_time,1,'last'); % index for interpolated vectors
    
    % iterate through indices that correspond to desired set
    set_indices = find(set_filter);
    for i = set_indices
        time_vec = spot_struct(i).time/60;        
        fluo_vec_interp = spot_struct(i).fluoInterp;
        time_vec_interp = spot_struct(i).timeInterp/60;
        obs_times = time_vec_interp(fluo_vec_interp>0);        
        if isempty(obs_times)
            obs_times = 0;
        end
        first_fluo_time_vec(i) = nanmin(obs_times);
        first_nuc_time_vec(i) = nanmin(time_vec);
        last_nuc_time_vec(i) = nanmax(time_vec);
        [~,closest_index] = min(abs(time_vec-switch_time));
        ap_pos_mean_vec(i) = spot_struct(i).APPosNucleus(closest_index);          
        
    end

    % apply position and obs number filter
    time_filter = first_fluo_time_vec < switch_time & last_nuc_time_vec >= ...
                        time_bounds(2) & first_nuc_time_vec <= time_bounds(1);
    ap_filter = ap_pos_mean_vec<=ap_bounds(2) & ap_pos_mean_vec>=ap_bounds(1);
    spot_struct_trunc = spot_struct(n_filter & ap_filter & time_filter);
    particle_id_filt = particle_id_vec(n_filter & ap_filter & time_filter);
        
    % now generate arrays containing fraction of active nuclei and knirps as a
    % function of time      
    io_ref_struct_temp = struct;
%     io_ref_struct_temp.on_off_array = NaN(length(index_vec),length(spot_struct_trunc));
    io_ref_struct_temp.fluo_array = NaN(length(index_vec),length(spot_struct_trunc));        
    io_ref_struct_temp.knirps_array = NaN(length(index_vec),length(spot_struct_trunc));
    io_ref_struct_temp.knirps_array_filled = NaN(length(index_vec),length(spot_struct_trunc));    
    
    % track wt for reactivation
    io_ref_struct_temp.reactivation_time_vec = NaN(1,length(spot_struct_trunc));

    % loop through qualifying particles
    for i = 1:length(spot_struct_trunc)

        % extract vectors
        time_vec_interp = spot_struct_trunc(i).timeInterp/60;        
        frame_vec = find(ismember(time_index_interp,time_vec_interp));
        frame_vec_rel = frame_vec - switch_i;
        time_vec_raw = spot_struct_trunc(i).time/60;
        knirps_vec = spot_struct_trunc(i).rawNCProtein / 1e5;
        
        if max(diff(time_vec_raw)) < 2 && mean(~isnan(knirps_vec)) > 0.9

            % actual data            
            fluo_vec_interp = spot_struct_trunc(i).fluoInterp;
            fluo_vec_orig = spot_struct_trunc(i).fluo;   
            nan_frames = isnan(knirps_vec);
            knirps_vec_interp = interp1(time_vec_raw(~nan_frames),knirps_vec(~nan_frames),time_index_interp);   
            time_filter = time_ref/60<=max(time_vec_raw-switch_time) & time_ref/60>=min(time_vec_raw-switch_time);
            
            % Check to see if the trace properly "turns off" prior to
            % perturbation            
            off_flag = all(fluo_vec_interp(ismember(frame_vec_rel,off_frame_ref))==0);            
            full_denom_count = full_denom_count + 1;
            
            % slightly less restrictive
            if off_flag                                                                                                 
                % estimate reactivation time
                ra_i = find(time_vec_raw > switch_time & ~isnan(fluo_vec_orig),1);
                if isempty(ra_i)
                    io_ref_struct_temp.reactivation_time_vec(i) = 1e6;                
                else
                    io_ref_struct_temp.reactivation_time_vec(i) = time_vec_raw(ra_i)-switch_time;                
                end
            end
      
            % generate indexing vectors            
            match_to = ismember(index_vec,frame_vec_rel);
            match_from = ismember(frame_vec_rel,index_vec);

            % save
            io_ref_struct_temp.fluo_array(time_filter,i) = 0;
            io_ref_struct_temp.fluo_array(match_to,i) = fluo_vec_interp(match_from);

            io_ref_struct_temp.knirps_array(match_to,i) = knirps_vec_interp(match_from);            

            % extend knirps profile to fill early and late periods
            io_ref_struct_temp.knirps_array_filled(:,i) = io_ref_struct_temp.knirps_array(:,i);            
            si = find(~isnan(io_ref_struct_temp.knirps_array(:,i)),1);
            fi = find(~isnan(io_ref_struct_temp.knirps_array(:,i)),1,'last');
            io_ref_struct_temp.knirps_array_filled(1:si-1,i) = io_ref_struct_temp.knirps_array_filled(si,i);
            io_ref_struct_temp.knirps_array_filled(fi+1:end,i) = io_ref_struct_temp.knirps_array_filled(fi,i);                                    
        end
    end    
    
    % cv_factor = k_after/k_before;
    io_ref_struct_temp.knirps_array_norm = io_ref_struct_temp.knirps_array_filled;
    io_ref_struct_temp.knirps_array_norm(window_size+1:end,:) = ...
      (io_ref_struct_temp.knirps_array_norm(window_size+1:end,:)-cal_intercept)/cal_slope;
    io_ref_struct_temp.knirps_array_norm = io_ref_struct_temp.knirps_array_norm - knirps_offset;

    % add some basic metadata
    io_ref_struct_temp.projectName = projectNameRA;
    io_ref_struct_temp.deltaT = spot_struct_trunc(1).tresInterp;
    io_ref_struct_temp.min_dp = min_dp;
    io_ref_struct_temp.ap_bounds = ap_bounds;
    io_ref_struct_temp.time_bounds = time_bounds;
    io_ref_struct_temp.time_vec = time_ref';
    io_ref_struct_temp.setID = setID;
    io_ref_struct_temp.switch_time = switch_time;
    io_ref_struct_temp.particle_id_vec = particle_id_filt;
    fnames = fieldnames(io_ref_struct_temp);
    for f = 1:length(fnames)
        master_struct(s).(fnames{f}) = io_ref_struct_temp.(fnames{f});
    end
    clear io_ref_struct_temp
end

% merge sets
fnames = fieldnames(master_struct);
io_ref_ra = struct;
for f = 1:length(fnames)
    io_ref_ra.(fnames{f}) = [master_struct.(fnames{f})];
end
io_ref_ra.set_index_full = floor(io_ref_ra.particle_id_vec);
io_ref_ra.time_vec = io_ref_ra.time_vec(:,1);
io_ref_ra.ap_bounds = io_ref_ra.ap_bounds(1:2);

% get rid of observations with NaN knirps values
nan_flags = any(isnan(io_ref_ra.knirps_array_norm));

io_ref_ra.knirps_array = io_ref_ra.knirps_array_norm(:,~nan_flags);
io_ref_ra = rmfield(io_ref_ra,'knirps_array_norm');
io_ref_ra.knirps_array(io_ref_ra.knirps_array<0) = 0;
io_ref_ra.fluo_array = io_ref_ra.fluo_array(:,~nan_flags);
io_ref_ra.knirps_array_filled = io_ref_ra.knirps_array_filled(:,~nan_flags);
io_ref_ra.reactivation_time_vec = io_ref_ra.reactivation_time_vec(~nan_flags);
io_ref_ra.particle_id_vec = io_ref_ra.particle_id_vec(~nan_flags);
io_ref_ra.set_index_full = io_ref_ra.set_index_full(~nan_flags);
io_ref_ra.set_index = unique(io_ref_ra.set_index_full);
io_ref_ra.off_frame_ref = off_frame_ref;

% use bootstrap resampling to estimate standard error
nBoots = 100;
max_ra_time =15*60;%ceil(max(ra_times)*60/Tres)*Tres;
ra_time_vec = 0:spot_struct(1).tresInterp:max_ra_time;
ra_time_array = NaN(nBoots,length(ra_time_vec));
ra_time_array_full = NaN(nBoots,length(ra_time_vec));
% sample_indices = 1:length(io_ref_ra.reactivation_time_vec);

for n = 1:nBoots
    % randomly select sets to include
    boot_sets = randsample(io_ref_ra.set_index,length(io_ref_ra.set_index),true);
    boot_indices = [];
    for s = 1:length(boot_sets)
        boot_indices = [boot_indices find(io_ref_ra.set_index_full==boot_sets(s))];
    end
%     boot_indices = randsample(find(io_ref_ra.set_index_full)sample_indices,length(sample_indices),true);
    reactivation_time_vec_boot = io_ref_ra.reactivation_time_vec(boot_indices);
    ra_times_boot = reactivation_time_vec_boot(~isnan(reactivation_time_vec_boot));
    [ra_times_sorted,~] = sort(ra_times_boot); 
    ra_count_raw = (0:length(ra_times_boot))/length(ra_times_boot);
    dummy_time_vec = [0 sort(ra_times_sorted*60+rand(size(ra_times_sorted))*1)];
    if max(dummy_time_vec) < ra_time_vec(end)
        dummy_time_vec(end+1) = ra_time_vec(end);
        ra_count_raw(end+1) = ra_count_raw(end);
    end
    ra_count_interp = interp1(dummy_time_vec,ra_count_raw,ra_time_vec,'previous');
    ra_time_array(n,:) = ra_count_interp;
    
    % reincorporate baseline
    ra_count_interp_full = ra_count_interp + sum(isnan(reactivation_time_vec_boot))/sum(~isnan(reactivation_time_vec_boot));
    ra_count_interp_full = ra_count_interp_full*sum(~isnan(reactivation_time_vec_boot))/length(reactivation_time_vec_boot);
    ra_time_array_full(n,:) = ra_count_interp_full;
end

% store results
io_ref_ra.reactivation_time_cdf = mean(ra_time_array);
io_ref_ra.reactivation_time_cdf_ste = std(ra_time_array) + 1e-2;
io_ref_ra.reactivation_time_cdf_full = mean(ra_time_array_full);
io_ref_ra.reactivation_time_cdf_ste_full = std(ra_time_array_full) + 1e-2;
io_ref_ra.reactivation_time_axis = ra_time_vec;

% % save
% save([resultsRoot 'io_ref_ra.mat'],'io_ref_ra')
