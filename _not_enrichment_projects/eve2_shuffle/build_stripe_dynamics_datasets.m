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
        else
            time_vec(time_vec<5) = NaN;
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

% save dataset
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
disp('Performing calculations for cumulative fraction ON/OFF...')

fraction_on_info = struct;

% generate space/time indexing vectors
fraction_on_info.max_time = hm_info_struct.max_time;
dT = 1;
fraction_on_info.time_index = 0:dT:fraction_on_info.max_time; % time range
fraction_on_info.ap_index = 20:1:56; % ap range
fraction_on_info.mean_ap_window = [20 30];
fraction_on_info.stripe_width = 2;

% initialize array
fraction_on_info.param_cell = {'fraction_on_abs','fraction_on_rel','fraction_off_rel','cum_mRNA_full','mean_mRNA_act','mean_frac_on'};
fraction_on_info.param_labels = {'cumulative fraction ON','cumulative fraction ON (relative)','cumulative fraction OFF','cumulative mRNA Full','mean mRNA active','mean fraction active'};
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
        trace_qc_filter = isnan([spot_struct.TraceQCFlag]) | [spot_struct.TraceQCFlag]==1;
        set_indices = find(set_vec==set_index(s)&trace_qc_filter);
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
        spot_struct_filtered = spot_struct(set_vec==set_index(s)&trace_qc_filter);
        spot_struct_filtered = spot_struct_filtered(stripe_center_filter);                
        
        % add flags to spot_struct
        iter_j = 1;
        for sp = set_indices
            spot_struct(sp).stripeCenterFlag = stripe_center_filter(iter_j);
            spot_struct(sp).stripeCenterAP = mean_stripe_ap;
            iter_j = iter_j + 1;
        end
        master_struct(p).spot_struct = spot_struct;
        
        % Make array to track activation status of each nucleus in stripe
        % center
        fluo_on_array = NaN(length(fraction_on_info.time_index)-1,length(spot_struct_filtered));
        fluo_off_array = NaN(length(fraction_on_info.time_index)-1,length(spot_struct_filtered));
        mean_fluo_full_array = NaN(length(fraction_on_info.time_index)-1,length(spot_struct_filtered));     
        mean_fluo_act_array = NaN(length(fraction_on_info.time_index)-1,length(spot_struct_filtered));   
        mean_fluo_array2 = NaN(length(fraction_on_info.time_index)-1,length(spot_struct_filtered));   
        for f = 1:size(fluo_on_array,2)
            time_vec = [spot_struct_filtered(f).time]/60;            
            fluo_vec = [spot_struct_filtered(f).fluo];         
%             if p == 3 && ismember(s,[1,3]) % NL: temporary fix for early FPs in rand1
%                 fluo_vec(time_vec<10) = NaN;
%             else
%                 fluo_vec(time_vec<5) = NaN;
%             end
            % find first and last times active
            fluo_vec_bin = ~isnan(fluo_vec);
            si = find(fluo_vec_bin,1);
            fi = find(fluo_vec_bin,1,'last');
            on_off_flag = min(time_vec) <= 10 & max(time_vec) >= max_set_time-60;
%             on_flag = ~isempty(si) && min(time_vec) <= 10 && si ~= 1;
            if on_off_flag && ~isempty(si) % if start observed and initial time reasonable
                start_i = find(time_vec(si)>hm_info_struct.time_index(1:end),1,'last');
                fluo_on_array(start_i:end,f) = 1;
                fluo_on_array(1:start_i-1,f) = 0; 
            elseif on_off_flag && isempty(si) % if no start observed 
                fluo_on_array(:,f) = 0; 
            end
            % turn OFF dynamics
            if on_off_flag && ~isempty(fi) % this basically filters out nuclei that never turned on
%                 off_flag = fi ~= length(time_vec) || time_vec(fi)==max_set_time; % this filters out nuclei that drifted out of frame...  
%                 if off_flag
                stop_i = find(time_vec(fi)<hm_info_struct.time_index(1:end-1),1,'first');
                fluo_off_array(stop_i+1:max_set_index,f) = 1;
                fluo_off_array(1:stop_i,f) = 0;
%                 end
            end
            
            % now deal with fluoresncence
            
            % initialize all observed time points to 0
            fluo_full = fluo_vec;
            fluo_full(isnan(fluo_full)) = 0;
            time_indices = discretize(time_vec,fraction_on_info.time_index);
            time_u = unique(time_indices);
            mean_fluo_full = grpstats(fluo_full,time_indices,'mean');%accumarray(time_indices',fluo_full',[],@mean);
            mean_fluo_full_array(time_u,f) = mean_fluo_full;                       
            
            if ~isempty(si)% && ~all(isnan(spot_struct_filtered(f).fluoInterp))
%                 fluo_act = fluo_vec(si:fi);
%                 time_act = time_vec(si:fi);
                fluo_act = spot_struct_filtered(f).fluoInterp;
                time_act = spot_struct_filtered(f).timeInterp/60;                
%                 fluo_act(isnan(fluo_act)) = 0;
                        
%                 if mean(fluo_act)*1.5 <= mean(fluo_act2) && ~ismember(f,[16 17 35 44]) && p>3
%                   error('wtf')
%                 end
                if ~all(isnan(time_act))
                  
                    time_indices = discretize(time_act,fraction_on_info.time_index);
                    time_u = unique(time_indices);
                    mean_fluo_act = grpstats(fluo_act,time_indices,'mean');%accumarray(time_indices',fluo_act',[],@mean);
                    mean_fluo_act_array(time_u,f) = mean_fluo_act;
%                     
%                     time_indices2 = discretize(time_act2,fraction_on_info.time_index);
%                     time_u = unique(time_indices2);
%                     mean_fluo_act2 = grpstats(fluo_act2,time_indices2,'mean');%accumarray(time_indices',fluo_act',[],@mean);
%                     mean_fluo_array2(time_u,f) = mean_fluo_act2;
                    
                end
            end
                
        end
        
        % find and interpolate over missing time points 
        pop_ids = all(~isnan(mean_fluo_full_array),2);
        first_i = find(pop_ids,1);
        last_i = find(pop_ids,1,'last');
        ref_index = 1:size(mean_fluo_act_array,1);
        gap_filter = all(isnan(mean_fluo_full_array),2) & ref_index' < last_i & ref_index' > first_i;
        if any(gap_filter)
            mean_fluo_full_array(gap_filter,:) = interp1(time_axis(~gap_filter)',mean_fluo_full_array(~gap_filter,:),time_axis(gap_filter)');
            mean_fluo_act_array(gap_filter,:) = interp1(time_axis(~gap_filter)',mean_fluo_act_array(~gap_filter,:),time_axis(gap_filter)');
        end 
        
        % store results                    
        fraction_on_info.frac_dynamics_array(:,iter_i,2) = nansum(fluo_on_array,2)/max(nansum(fluo_on_array,2));
        n_obs = sum(max(~isnan(fluo_on_array)));
        fraction_on_info.frac_dynamics_array(:,iter_i,1) = nansum(fluo_on_array,2)/n_obs;
        
        % get denominator for off nuclei
        off_denom = sum(~isnan(fluo_off_array(1,:)),2);
        off_trend = nansum(fluo_off_array,2)/off_denom;
        fraction_on_info.frac_dynamics_array(1:max_set_index,iter_i,3) = off_trend(1:max_set_index);
        
        % deal with mRNA 
        mf_full_vec = nanmean(mean_fluo_full_array,2);
        mf_full_vec(isnan(mf_full_vec) & (1:length(mf_full_vec))'<=max_set_index) = 0;
        fraction_on_info.frac_dynamics_array(:,iter_i,4) = cumsum(mf_full_vec);
        
        % cumulative active time steps 
        n_act_vec = nansum(~isnan(mean_fluo_act_array),2);
        n_total_vec = nansum(~isnan(mean_fluo_full_array),2);
        
        frac_active_vec = n_act_vec ./ n_total_vec;
        frac_active_vec(isnan(frac_active_vec)) = 0;
        fraction_on_info.frac_dynamics_array(:,iter_i,6) = cumsum(frac_active_vec);
        
        % now use these two pieces to infer the average fluorescence
        fraction_on_info.frac_dynamics_array(:,iter_i,5) = ...
                fraction_on_info.frac_dynamics_array(:,iter_i,4)./fraction_on_info.frac_dynamics_array(:,iter_i,6);
        
%         mf_act_vec(isnan(mf_act_vec) & (1:length(mf_full_vec))'<=max_set_index) = 0;
%         fraction_on_info.frac_dynamics_array(:,iter_i,4) = cumsum(mf_full_vec);
%         
%         for z = 1:size(mean_fluo_act_array,1)
%             fraction_on_info.frac_dynamics_array(z,iter_i,5) = nanmean(nanmean(mean_fluo_act_array(1:z,:)),2);
%             fraction_on_info.frac_dynamics_array(z,iter_i,6) = nanmean(nanmean(~isnan(mean_fluo_act_array(1:z,:)),2));
%         end

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


%{
%% note that we must extrapolate for one WT MSE set to estimate HM point

% Need to revisit this
extrap_ids = [2,4];
fraction_on_info.frac_dynamics_array_ex = fraction_on_info.frac_dynamics_array;
for e = extrap_ids
    extrap_vec = fraction_on_info.frac_dynamics_array_ex(:,e,3);
    time_ex = time_axis(~isnan(extrap_vec));
    extrap_vec_out = interp1(time_ex,extrap_vec(~isnan(extrap_vec)),time_axis,'linear','extrap');
    fraction_on_info.frac_dynamics_array_ex(:,e,3) = extrap_vec_out;
end
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
%}

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