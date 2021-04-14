% script to build data sets for stripe dynamics figures, as well as
% subsequent analyses
clear
close all

projectNameCell = {'MSE-WT','NSv1','Rand1','Rand4'};
master_struct = struct;
hm_info_struct = struct;

for p = 1:length(projectNameCell)
    projectName = projectNameCell{p};
    liveProject = LiveEnrichmentProject(projectName);
    resultsRoot = [liveProject.dataPath filesep];

    % load data
    load([resultsRoot 'spot_struct.mat'])
    master_struct(p).spot_struct = spot_struct;
    master_struct(p).projectName = projectName;
    clear spot_struct
end

slashesData = strfind(liveProject.dataPath,'\');
DataPath = [liveProject.dataPath(1:slashesData(end-1)) 'eve2_recon_analyses' filesep];
mkdir(DataPath);

hm_info_struct.legend_str = {'WT-MSE2','Neutral Sequence v1', 'Randomized v1','Randomized v4'};
hm_info_struct.legend_str_short = {'MSE2','NS v1', 'Rand v1','Rand v4'};
hm_info_struct.projectNameCell = projectNameCell;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate heatmaps showing the average activity over space and time

disp('performing calculations for stripe dynamics heatmaps...')

% generate space/time indexing vectors
hm_info_struct.max_time = 60;
dT = 1;
hm_info_struct.time_index = 0:dT:hm_info_struct.max_time; % time range
hm_info_struct.ap_index = 20:1:56; % ap range
hm_info_struct.ap_sigma = 0.75;
hm_info_struct.time_sigma = 0.75;
% determine how many embryos there are total
hm_info_struct.embryo_id_vec = [];
hm_info_struct.project_id_vec = [];

for p = 1:length(master_struct)
    % extract and generate indexing vector
    spot_struct = master_struct(p).spot_struct;
    set_vec = [spot_struct.setID];
    set_index = unique(set_vec);
    % update
    hm_info_struct.embryo_id_vec = [hm_info_struct.embryo_id_vec set_index];
    hm_info_struct.project_id_vec = [hm_info_struct.project_id_vec repelem(p,length(set_index))];
end    
    
hm_info_struct.legend_ids = [];
for s = 1:length(master_struct)
    hm_info_struct.legend_ids(s) = find(hm_info_struct.project_id_vec==s,1);
end

% initialize array
hm_info_struct.full_dynamics_array = NaN(length(hm_info_struct.time_index)-1,length(hm_info_struct.ap_index)-1,length(hm_info_struct.project_id_vec));
hm_info_struct.offset_dynamics_array = NaN(length(hm_info_struct.time_index)-1,length(hm_info_struct.ap_index)-1,length(hm_info_struct.project_id_vec));
hm_info_struct.instant_frac_dynamics_array = NaN(length(hm_info_struct.time_index)-1,length(hm_info_struct.ap_index)-1,length(hm_info_struct.project_id_vec));
hm_info_struct.fluo_dynamics_array = NaN(length(hm_info_struct.time_index)-1,length(hm_info_struct.ap_index)-1,length(hm_info_struct.project_id_vec));

% loop through the embryos
iter_i = 1;
for p = 1:length(master_struct)
  
    % extract and generate indexing vector
    spot_struct = master_struct(p).spot_struct;
    set_vec = [spot_struct.setID];
    set_index = unique(set_vec);
    
    % iterate through the sets    
    for s = 1:length(set_index)
      
        % extract key vectors
        time_vec = [spot_struct(set_vec==set_index(s)).time]/60;
        if p == 3 && ismember(s,[1,3]) % NL: temporary fix for early FPs in rand1
            time_vec(time_vec<10) = NaN;
        end
        ap_vec = [spot_struct(set_vec==set_index(s)).APPosNucleus]*100;
        fluo_vec = [spot_struct(set_vec==set_index(s)).fluo];
        off_vec = [spot_struct(set_vec==set_index(s)).fluoOffset];
        fluo_vec_bin = ~isnan(fluo_vec);
        fluo_vec_active = fluo_vec;
        fluo_vec(isnan(fluo_vec)) = 0;
               
        for a = 1:length(hm_info_struct.ap_index)-1
            for t = 1:length(hm_info_struct.time_index)-1
                center_ap = mean(hm_info_struct.ap_index(a:a+1));
                center_time = mean(hm_info_struct.time_index(t:t+1));
                
                % calculate weights
                ap_diffs = (center_ap-ap_vec)./hm_info_struct.ap_sigma /sqrt(2);
                ap_filter = abs(ap_diffs) <= sqrt(2);
                time_diffs = (center_time-time_vec)./hm_info_struct.time_sigma /sqrt(2);
                time_filter = abs(time_diffs) <= sqrt(2);  
                wt_vec = exp(-0.5*(ap_diffs(ap_filter&time_filter).^2 + time_diffs(ap_filter&time_filter).^2));
                
                % assign values
                hm_info_struct.full_dynamics_array(t, a, iter_i) = nansum(wt_vec.*fluo_vec(ap_filter&time_filter))/nansum(wt_vec);
                hm_info_struct.instant_frac_dynamics_array(t, a, iter_i) = nansum(wt_vec.*fluo_vec_bin(ap_filter&time_filter))/nansum(wt_vec);
                hm_info_struct.fluo_dynamics_array(t, a, iter_i) = nansum(wt_vec.*fluo_vec_active(ap_filter&time_filter))/nansum(wt_vec);
                hm_info_struct.offset_dynamics_array(t, a, iter_i) = nanmean(off_vec(ap_filter&time_filter));
            end
        end
        iter_i = iter_i + 1;
    end
end

% save datasets
% save([DataPath 'hm_info_struct.full_dynamics_array.mat'],'hm_info_struct.full_dynamics_array')
% save([DataPath 'hm_info_struct.instant_frac_dynamics_array.mat'],'hm_info_struct.instant_frac_dynamics_array')
% save([DataPath 'hm_info_struct.fluo_dynamics_array.mat'],'hm_info_struct.fluo_dynamics_array')
% save([DataPath 'hm_info_struct.offset_dynamics_array.mat'],'hm_info_struct.offset_dynamics_array')
hm_info_struct.projectNameCell = projectNameCell;
save([DataPath 'hm_info_struct.mat'],'hm_info_struct')
disp('Done.')

% %%%%%%%%%%%%%%% Fit Gaussians to track stripes over time %%%%%%%%%%%%%%%
disp('extracting stripe parameters for each embryo...')

stripe_param_struct = struct;

% define indices to fit
fit_indices = 1:length(hm_info_struct.time_index)-1;
ap_axis = hm_info_struct.ap_index(1:end-1) + diff(hm_info_struct.ap_index)/2;

% define initial conditions and bounds
init_vec = [1e4 38 3];
ub_vec = [5e5 ap_axis(end) 20];
lb_vec = [ap_axis(1) 1 0.5];

% Initialize param array
fit_param_array = NaN(length(fit_indices),3,length(hm_info_struct.project_id_vec));
fit_param_raw_array = NaN(length(fit_indices),4,length(hm_info_struct.project_id_vec));
fit_profile_array = NaN(size(hm_info_struct.full_dynamics_array));

options = optimoptions(@lsqnonlin,'Display','off');

for p = 1:length(hm_info_struct.project_id_vec)    
    for f = 1:length(fit_indices)
        % extract AP profile
        activity_vec = hm_info_struct.full_dynamics_array(fit_indices(f),:,p);
        activity_vec(isnan(activity_vec)) = 0;

        % define function
        gauss_fun = @(x) x(1)*exp(-0.5*((x(2)-ap_axis)./x(3)).^2);
        ob_fun = @(x) gauss_fun(x)-activity_vec;

        % perform fit
        fit_param_array(f,:,p) = lsqnonlin(ob_fun,init_vec,lb_vec,ub_vec,options);
        fit_profile_array(fit_indices(f),:,p) = gauss_fun(fit_param_array(f,:,p));
        
        % use robust metrics as alternative
        wt_center = nansum(activity_vec.*ap_axis)./nansum(activity_vec);
        fit_param_raw_array(f,2,p) = wt_center;
        wstd = sqrt(var(ap_axis,activity_vec));
        fit_param_raw_array(f,3,p) = 2*wstd;
        awt = nanmean(activity_vec(ap_axis>=wt_center-wstd&ap_axis<=wt_center+wstd));
        fit_param_raw_array(f,1,p) = awt;
        
        % check deviation from gaussian
        gauss_pd = gauss_fun([awt wt_center wstd]);
        fit_param_raw_array(f,4,p) = sum((gauss_pd-activity_vec).^2);
    end    
end
    
% add to structure and save
stripe_param_struct.time_index = hm_info_struct.time_index;
time_axis = stripe_param_struct.time_index(1:end-1) + diff(stripe_param_struct.time_index)/2;
stripe_param_struct.fit_indices = fit_indices;
stripe_param_struct.time_to_plot = time_axis(fit_indices);
stripe_param_struct.ap_axis = ap_axis;
stripe_param_struct.fit_param_array = fit_param_array;
stripe_param_struct.fit_param_raw_array = fit_param_raw_array;
stripe_param_struct.fit_profile_array = fit_profile_array;


% First calculate average and standard deviation for each genotype
mean_param_array = NaN(size(fit_param_raw_array,1),size(fit_param_raw_array,2),length(master_struct));
std_param_array = NaN(size(fit_param_raw_array,1),size(fit_param_raw_array,2),length(master_struct));
for m = 1:length(master_struct)
    p_filter = hm_info_struct.project_id_vec == m;
    mean_param_array(:,:,m) = nanmean(fit_param_raw_array(:,:,p_filter),3);
    std_param_array(:,:,m) = nanstd(fit_param_raw_array(:,:,p_filter),[],3);
end    

stripe_param_struct.mean_param_array = mean_param_array;
stripe_param_struct.std_param_array = std_param_array;

% save
save([DataPath 'stripe_param_struct.mat'],'stripe_param_struct')
disp('Done.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Look at Fraction ON, time ON and OFF, and mean spot intensity
disp('Performing calcu;ations for cumulative fraction ON/OFF...')

fraction_on_info = struct;

% generate space/time indexing vectors
fraction_on_info.max_time = hm_info_struct.max_time;
dT = 1;
fraction_on_info.time_index = 0:dT:fraction_on_info.max_time; % time range
fraction_on_info.ap_index = 20:1:56; % ap range
fraction_on_info.mean_ap_window = [20 30];
fraction_on_info.stripe_width = 2;

% initialize array
fraction_on_info.param_cell = {'fraction_on_abs','fraction_on_rel','fraction_off_rel','cum_mRNA_full','cum_mRNA_act'};
fraction_on_info.param_labels = {'cumulative fraction ON','cumulative fraction ON (relative)','cumulative fraction OFF','cumulative mRNA Full','Cumulative mRNA active'};
fraction_on_info.frac_dynamics_array = NaN(length(fraction_on_info.time_index)-1,length(hm_info_struct.project_id_vec),...
                                        length(fraction_on_info.param_cell));

% loop through the embryos
iter_i = 1;
for p = 1:length(master_struct)
  
    % extract and generate indexing vector
    spot_struct = master_struct(p).spot_struct;
    set_vec = [spot_struct.setID];
    set_index = unique(set_vec);
    
    % iterate through the sets    
    for s = 1:length(set_index)
      
        % define average nucleus AP as the mean position from 15-25 minutes
        set_indices = find(set_vec==set_index(s));
        max_set_time = nanmax([spot_struct(set_indices).time])/60;
        max_set_index = find(max_set_time<fraction_on_info.time_index(1:end-1),1,'first');
        
        ap_vec = NaN(1,length(set_indices));                
        ap_steps = zeros(1,length(set_indices));                
        for a = 1:length(ap_vec)
            time_temp = spot_struct(set_indices(a)).time/60;
            time_filter = time_temp >= fraction_on_info.mean_ap_window(1) & time_temp <= fraction_on_info.mean_ap_window(2);
            if ~isempty(time_filter)
                ap_vec(a) = nanmean(spot_struct(set_indices(a)).APPosNucleus(time_filter))*100;
                ap_steps(a) = sum(time_filter);
            end
        end
        ep_filter = hm_info_struct.embryo_id_vec == s & hm_info_struct.project_id_vec == p;
        
        % find the embryo-specific stripe center
        stripe_pos_vec = fit_param_raw_array(:,2,ep_filter);        
        mean_stripe_ap = nanmean(stripe_pos_vec(stripe_param_struct.time_to_plot >= fraction_on_info.mean_ap_window(1) & ...
                                              stripe_param_struct.time_to_plot <= fraction_on_info.mean_ap_window(2)));
                       
        % recenter
        ap_vec_norm = ap_vec - mean_stripe_ap;
        stripe_center_filter = abs(ap_vec_norm) <= fraction_on_info.stripe_width & ap_steps > 10; 
        
        % extract relevant vectors        
        spot_struct_filtered = spot_struct(set_vec==set_index(s));
        spot_struct_filtered = spot_struct_filtered(stripe_center_filter);
        % add flags to spot_struct
        iter_j = 1;
        for sp = set_indices
            spot_struct(sp).stripeCenterFlag = stripe_center_filter(iter_j);
            spot_struct(sp).stripeCenterAP = mean_stripe_ap;
            iter_j = iter_j + 1;
        end
        master_struct(p).spot_struct = spot_struct;
        % make array to track activation status of each nucleus in stripe
        % center
        fluo_on_array = zeros(length(fraction_on_info.time_index)-1,length(spot_struct_filtered));
        fluo_off_array = NaN(length(fraction_on_info.time_index)-1,length(spot_struct_filtered));
        cum_fluo_array = NaN(length(fraction_on_info.time_index)-1,length(spot_struct_filtered));       
        for f = 1:size(fluo_on_array,2)
            time_vec = [spot_struct_filtered(f).time]/60;            
            fluo_vec = [spot_struct_filtered(f).fluo];         
            if p == 3 && ismember(s,[1,3]) % NL: temporary fix for early FPs in rand1
                fluo_vec(time_vec<10) = NaN;
            end
            fluo_vec_bin = ~isnan(fluo_vec);
            fluo_vec_tot = fluo_vec;
            fluo_vec_tot(isnan(fluo_vec_tot)) = 0;
            
            % find first and last times active
            si = find(fluo_vec_bin,1);
            fi = find(fluo_vec_bin,1,'last');
            if ~isempty(si) && min(time_vec) <= 5
                start_i = find(time_vec(si)>hm_info_struct.time_index(1:end),1,'last');
%                 stop_i = find(time_vec(fi)>hm_info_struct.time_index,1,'last');
                fluo_on_array(start_i:end,f) = 1;
            end
                        
            % cumulative mRNA
            cum_fluo_tot = fluo_vec_tot;%cumsum(fluo_vec_tot);
            ttl_interp = interp1(time_vec,cum_fluo_tot,time_axis);
            ttl_interp(time_axis<time_vec(1)) = 0;
            cum_fluo_array(:,f) = ttl_interp;
            cum_fluo_array(max_set_index+1:end,f) = NaN;
            % turn OFF dynamics
            if ~isempty(fi) % this basically filters out nuclei that never turned on
                if fi ~= length(time_vec) || time_vec(fi)==max_set_time % this filters out nuclei that drifted out of frame...
  %                 start_i = find(time_vec(si)>hm_info_struct.time_index,1,'last');
                  stop_i = find(time_vec(fi)<hm_info_struct.time_index(1:end-1),1,'first');
                  fluo_off_array(stop_i+1:max_set_index,f) = 1;
                  fluo_off_array(1:stop_i,f) = 0;
                  
                  % deal with cumulative mRNA
                  cum_fluo_array(stop_i+1:max_set_index,f) = cum_fluo_array(stop_i,f);
                elseif fi ~= length(time_vec)
                  stop_i = find(time_vec(fi)<hm_info_struct.time_index(1:end-1),1,'first');
                  ttl_interp(stop_i+1:max_set_index) = NaN;
                end
            end
                      
            
        end

        % store results                    
        fraction_on_info.frac_dynamics_array(:,iter_i,2) = sum(fluo_on_array,2)/max(sum(fluo_on_array,2));
        fraction_on_info.frac_dynamics_array(:,iter_i,1) = sum(fluo_on_array,2)/size(fluo_on_array,2);
        
        % get denominator for off nuclei
        off_denom = sum(~isnan(fluo_off_array(1,:)),2);
        off_trend = nansum(fluo_off_array,2)/off_denom;
        fraction_on_info.frac_dynamics_array(1:max_set_index,iter_i,3) = off_trend(1:max_set_index);
        
        % deal with mRNA subtleties
        fraction_on_info.frac_dynamics_array(:,iter_i,4) = cumsum(nanmean(cum_fluo_array,2));        
        n_extant = sum(fluo_on_array & fluo_off_array~=1,2);
        fraction_on_info.frac_dynamics_array(:,iter_i,5) = cumsum(nansum(cum_fluo_array,2)./n_extant);
        % get cumulative amongst just active nuclei
%         n_on = sum(fluo_on_array & ~fluo_off_array,2);
%         cumulative_mRNA_array_act(:,iter) = nansum(cum_fluo_array,2);
        iter_i = iter_i + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% take average for each genotype
% First calculate average and standard deviation for each genotype
mean_frac_array = NaN(size(fraction_on_info.frac_dynamics_array,1),length(master_struct),length(fraction_on_info.param_cell));
std_frac_array = NaN(size(fraction_on_info.frac_dynamics_array,1),length(master_struct),length(fraction_on_info.param_cell));
fraction_on_info.loc_str = {'southeast','southeast','northwest'};

for m = 1:length(master_struct)
    p_filter = hm_info_struct.project_id_vec == m;
    mean_frac_array(:,m,:) = nanmean(fraction_on_info.frac_dynamics_array(:,p_filter,:),2);
    std_frac_array(:,m,:) = nanstd(fraction_on_info.frac_dynamics_array(:,p_filter,:),[],2);
end   

fraction_on_info.mean_frac_array = mean_frac_array;
fraction_on_info.std_frac_array = std_frac_array;


% note that we must extrapolate for one WT MSE set to estimate HM point

% Need to revisit this
extrap_id = 4;
fraction_on_info.frac_dynamics_array_ex = fraction_on_info.frac_dynamics_array;
extrap_vec = fraction_on_info.frac_dynamics_array_ex(:,extrap_id,3);
time_ex = time_axis(~isnan(extrap_vec));
extrap_vec_out = interp1(time_ex,extrap_vec(~isnan(extrap_vec)),time_axis,'cubic','extrap');
fraction_on_info.frac_dynamics_array_ex(:,extrap_id,3) = extrap_vec_out;

% now find 50% on and off times
mean_on_50 = NaN(1,length(master_struct));
on_50 = cell(1,length(master_struct));
ste_on_50 = NaN(1,length(master_struct));
mean_off_50 = NaN(1,length(master_struct));
off_50 = cell(1,length(master_struct));
ste_off_50 = NaN(1,length(master_struct));


for m = 1:length(master_struct)
    p_ids = find(hm_info_struct.project_id_vec == m);
    on_times = [];
    off_times = [];
    for p = 1:length(p_ids)
        % on time
        on_i = find(fraction_on_info.frac_dynamics_array_ex(:,p_ids(p),2)>=0.5,1);
        on_times(end+1) = time_axis(on_i);
        % off time
        off_i = find(fraction_on_info.frac_dynamics_array_ex(:,p_ids(p),3)>=0.5,1);
        off_times(end+1) = time_axis(off_i);
    end
    on_50{m} = on_times;
    mean_on_50(m) = mean(on_times);
    ste_on_50(m) = std(on_times);
    
    off_50{m} = off_times;
    mean_off_50(m) = mean(off_times);
    ste_off_50(m) = std(off_times);    
end   

fraction_on_info.on_50 = on_50;
fraction_on_info.mean_on_50 = mean_on_50;
fraction_on_info.ste_on_50 = ste_on_50;

fraction_on_info.off_50 = off_50;
fraction_on_info.mean_off_50 = mean_off_50;
fraction_on_info.ste_off_50 = ste_off_50;

% save
save([DataPath 'fraction_on_info.mat'],'fraction_on_info')
disp('Done.')

% save updated spot structures
clear spot_struct
disp('saving updated spot structures...')
for p = 1:length(projectNameCell)
    projectName = projectNameCell{p};
    liveProject = LiveEnrichmentProject(projectName);
    resultsRoot = [liveProject.dataPath filesep];
    spot_struct = master_struct(p).spot_struct;
    % load data
    save([resultsRoot 'spot_struct.mat'],'spot_struct')    
    clear spot_struct
end
disp('Done.')