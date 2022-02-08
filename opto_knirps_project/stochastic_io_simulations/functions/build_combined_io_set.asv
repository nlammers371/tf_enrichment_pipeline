function io_struct_io = build_combined_io_set(dataRoot)

addpath(genpath('./lib'))
addpath(genpath('../../utilities'))

% Load data
% projectNameWT = 'optokni_eve4+6_WT_FUN'; 
% clear
close all

addpath(genpath('./lib'))
addpath(genpath('../../utilities'))

% indicate projects to use
projectNameCell = {'optokni_eve4+6_ON_LOW_FULL','optokni_eve4+6_WT','optokni_eve4+6_ON_CONST'};

% specify correction parameters for blue light
knirps_offset = 375698.13 / 1e5;
cal_slope = 1.243;
cal_intercept = 1.079e5 / 1e5; %NL: dividing everything through by 1e5 for simplicity

% define key cleaning parameters
time_range = [10 30]*60; % all nuclei considered must be observed for full period
use_range = [15 30]*60; % all nuclei considered must be observed for full period
ap_range = [-0.01 0.01];
min_dp = 10; % must have at least 10 active frames observed

% initialize data structure
io_struct_io = struct;

% load data
spot_struct_full = [];
for p = 1:length(projectNameCell)
    projectName = projectNameCell{p};    
    load([dataRoot projectName filesep 'spot_struct.mat'],'spot_struct')
    for s = 1:length(spot_struct)
        spot_struct(s).projectName = projectName;
        spot_struct(s).projectID = p;
    end
    spot_struct_full = [spot_struct_full spot_struct];
end

spot_struct = spot_struct_full;

% generate a new set of unique identifiers
set_vec = [spot_struct.setID];
project_vec = [spot_struct.projectID];
project_set_array = unique([[spot_struct.projectID]' [spot_struct.setID]'],'rows');

for i = 1:size(project_set_array,1)
    projectID = project_set_array(i,1);
    setID = project_set_array(i,2);
    temp_ids = find(project_vec==projectID&set_vec==setID);
    for t = temp_ids
       spot_struct(t).masterID = i;
       particleID = spot_struct(t).particleID;
       particleIDNew = particleID - floor(particleID) + i;
       spot_struct(t).particleID = particleIDNew;
       spot_struct(t).particleIDOrig = particleID;
    end
end    

% identify and correct for effects of blue light laser
% spot_struct_temp = spot_struct;
master_id_vec = [spot_struct.masterID];
master_id_index = unique(master_id_vec);
time_index_interp = 0:spot_struct(1).tresInterp:(40*60);

%% NL: obtained these frames via manual inspection of protein trends
blue_light_frame_vec = [NaN(1,8) 41 33 36];

% initialize arrays to store results
io_struct_io.master_id_vec = master_id_vec;
io_struct_io.particle_id_vec = [spot_struct.particleID];
io_struct_io.time_axis = time_index_interp;

io_struct_io.mean_ap_vec = NaN(1,length(master_id_index));
io_struct_io.knirps_array_raw = NaN(length(time_index_interp),length(master_id_vec));
io_struct_io.knirps_array = NaN(length(time_index_interp),length(master_id_vec));
io_struct_io.fluo_array = NaN(length(time_index_interp),length(master_id_vec));
io_struct_io.mean_knirps_array = NaN(length(time_index_interp),length(master_id_index));
io_struct_io.mean_fluo_array = NaN(length(time_index_interp),length(master_id_index));
io_struct_io.project_id_index = NaN(1,length(master_id_index));

keep_flags = false(1,length(master_id_vec));

for m = 1:length(master_id_index)
    master_ids = find(master_id_vec==master_id_index(m));
    time_index = unique(round([spot_struct(master_ids).time],0));    
    io_struct_io.project_id_index(m) = spot_struct(master_ids(1)).projectID;
    % find changepoint  
    shift_frame = blue_light_frame_vec(m);
    if ~isnan(shift_frame)
        shift_time = time_index(shift_frame);               
    else
        shift_time = Inf;
    end
    
    % check that corrections look reasonable
    protein_array_temp = NaN(length(time_index),length(master_ids));        
    for i = master_ids
        t_vec = round(spot_struct(i).time,0);
        f_vec_raw = spot_struct(i).fluo;
        t_vec_fluo = spot_struct(i).timeInterp;
        start_i = find(time_index_interp<=t_vec(1)&time_index_interp<=t_vec_fluo(1),1,'last');
        stop_i = find(time_index_interp>=t_vec(end)&time_index_interp>=t_vec_fluo(end),1);
        t_vec_interp = time_index_interp(start_i:stop_i);           

        % extract basic vectors
        pt_vec = spot_struct(i).rawNCProtein/1e5;
        fluo_vec = spot_struct(i).fluoInterp;        
        ap_vec = spot_struct(i).APPosNucleus;
        time_ft = t_vec>=time_range(1)&t_vec<=time_range(2);
        mean_ap = nanmean(ap_vec(time_ft));
        
        % check that this trace qualifies for keeping
        time_keep_flag = t_vec(1)<=time_range(1)&&t_vec(end)>=time_range(2);                
        ap_keep_flag = mean_ap>=ap_range(1)&&mean_ap<=ap_range(2);
        dp_keep_flag = sum(~isnan(f_vec_raw))>=min_dp;
        
        if time_keep_flag && ap_keep_flag && dp_keep_flag &&~any(isnan(pt_vec))
            % update keep flag
            keep_flags(i) = true;
            time_filter = ismember(time_index_interp,t_vec_interp);
            start_i = find(time_filter,1);
            stop_i = find(time_filter,1,'last');
            
            % generate adjusted pt vector
            raw_adjusted = pt_vec - knirps_offset;
            pert_ind = find(t_vec>=shift_time,1);
            if ~isempty(pert_ind) 
                raw_adjusted(pert_ind:end) = ...
                (raw_adjusted(pert_ind:end)-cal_intercept)/cal_slope;
            end
            
            % interpolate protein and asign to array
            pt_interp = interp1(t_vec,raw_adjusted,t_vec_interp,'linear','extrap');
            io_struct_io.knirps_array_raw(time_filter,i) = pt_interp;
            
            % fill obs fowards and backwards
            io_struct_io.knirps_array(:,i) = io_struct_io.knirps_array_raw(:,i);
            io_struct_io.knirps_array(1:start_i-1,i) = io_struct_io.knirps_array(start_i,i);
            io_struct_io.knirps_array(stop_i+1:end,i) = io_struct_io.knirps_array(stop_i,i);
            
            % initialize fluo with zeros
            io_struct_io.fluo_array(:,i) = 0;
            io_struct_io.fluo_array(ismember(time_index_interp, t_vec_fluo),i) = fluo_vec;                        
        end     
    end        

    % check that corrections look reasonable    
    io_struct_io.mean_protein_array(:,m) = nanmean(io_struct_io.knirps_array(:,master_ids),2);    
    io_struct_io.mean_fluo_array(:,m) = nanmean(io_struct_io.fluo_array(:,master_ids),2);    
end    

% Apply space/time filters
io_struct_io.knirps_array_ft = io_struct_io.knirps_array(:,keep_flags);
io_struct_io.fluo_array_ft = io_struct_io.fluo_array(:,keep_flags);
io_struct_io.master_id_vec_ft = master_id_vec(keep_flags);
io_struct_io.particle_id_vec_ft = io_struct_io.particle_id_vec(keep_flags);

% perform bootstrap resampling to obtain estimates of mean and standard
% error
nBoots = 100;
nBins = 20;

% set_index = unique(io_struct_io.master_id_vec_ft);
mean_filter = time_index_interp>=use_range(1)&time_index_interp<=use_range(2);

% generate knirps axis
% knirps_vec_long = io_struct_io.knirps_array(mean_filter,:);
% knirps_vec_long = knirps_vec_long(:);
knirps_axis = linspace(2,12,nBins+1);

fluo_array_temp = NaN(length(knirps_axis)-1,nBoots);
% for each iteration, randomly select embryos to include, and then which
%% nuclei within those projects
for n = 1:nBoots
    % first, randomly select subset to sample
    boot_set_ids = randsample(master_id_index,length(master_id_index),true);
    set_options = [];
    for s = 1:length(boot_set_ids)
        set_options = [set_options find(io_struct_io.master_id_vec_ft==boot_set_ids(s))];
    end
    boot_trace_ids = randsample(set_options,length(io_struct_io.master_id_vec_ft),true);
%     boot_trace_ids = randsample(1:length(io_struct_io.master_id_vec_ft),length(io_struct_io.master_id_vec_ft),true);
    % obtain fluo and knirps arrays
    boot_fluo_vec = io_struct_io.fluo_array_ft(mean_filter,boot_trace_ids);
    boot_fluo_vec = boot_fluo_vec(:);
    boot_knirps_vec = io_struct_io.knirps_array_ft(mean_filter,boot_trace_ids);
    boot_knirps_vec = boot_knirps_vec(:);
    % now iterate through knirps values
    for k = 1:length(knirps_axis)-1
        kni_filter = boot_knirps_vec>=knirps_axis(k) & boot_knirps_vec<knirps_axis(k+1);
        fluo_array_temp(k,n) = nanmean(boot_fluo_vec(kni_filter));
    end
end

io_struct_io.mean_knirps_trend = knirps_axis(1:end-1) + diff(knirps_axis);
io_struct_io.mean_fluo_trend = nanmean(fluo_array_temp,2);
io_struct_io.ste_fluo_trend = nanstd(fluo_array_temp,[],2);
% %%
% io_struct_io.knirps_array = io_struct_io.knirps_array(:,keep_flags==1);
% io_struct_io.cumulative_protein_array = cumsum(io_struct_io.mean_protein_array);
% 
% time_array = repmat(time_index_interp',1,length(master_id_index));
% set_id_array = repmat(io_struct_io.project_id_index,length(time_index_interp),1);
% time_filter = time_array>=15*60&time_array<=30*60;
% 
% % close all
% figure;
% cmap = brewermap(3,'Set2');
% colormap(cmap);
% scatter(io_struct_io.mean_protein_array(time_filter),io_struct_io.mean_fluo_array(time_filter),50,set_id_array(time_filter))
% % h = colorbar('XTick', 1:3);
% xlabel('[Knirps]')
% ylabel('mean activity')
% % %%
% % close all
% % figure;
% % cmap = brewermap(3,'Set2');
% % colormap(cmap);
% % scatter(io_struct_io.mean_protein_array(time_filter),io_struct_io.cumulative_protein_array(time_filter),50,set_id_array(time_filter))
% % colorbar