clear
close all

addpath(genpath('./lib'))
addpath(genpath('../../utilities'))

% Load data
dataRootJake = 'S:\Jake\Dropbox\ProcessedEnrichmentData\';
dataRoot = 'S:\Nick\Dropbox (Personal)\ProcessedEnrichmentData\';

% make new director for combined dataset
savePath = [dataRoot 'combinedOptoSets' filesep];
mkdir(savePath)

% indicate projects to use
projectNameCell = {'optokni_eve4+6_ON_LOW_FULL','optokni_eve4+6_WT','optokni_eve4+6_ON_CONST'};

% specify correction parameters for blue light
knirps_offset = 375000 / 1e5;
cal_slope = 1.243;
cal_intercept = 1.079e5 / 1e5; %NL: dividing everything through by 1e5 for simplicity

% load data
spot_struct_full = [];
for p = 1:length(projectNameCell)
    projectName = projectNameCell{p};    
    load([dataRootJake projectName filesep 'spot_struct.mat'],'spot_struct')
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
minDP = 20;
spot_struct_temp = spot_struct;
master_id_vec = [spot_struct.masterID];
master_id_index = unique(master_id_vec);
time_index_interp = 0:spot_struct(1).tresInterp:(50*60);

% NL: obtained these frames via manual inspection of protein trends
blue_light_frame_vec = [NaN(1,8) 41 33 36];
% blue_light_frame_vec = [52 NaN 52 73 46 NaN 50 36 38 miNaN 42 43 NaN 67 NaN 45 43];% OG

% initialize array to store mean protein trend
mean_protein_array = NaN(length(time_index_interp),length(master_id_index));

for m = 1:length(master_id_index)
    master_ids = find(master_id_vec==master_id_index(m));
    time_index = unique(round([spot_struct(master_ids).time],0));
    projectName = spot_struct(master_ids(1)).projectName;
                      
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
        t_vec_fluo = spot_struct(i).timeInterp;
        start_i = find(time_index_interp<=t_vec(1)&time_index_interp<=t_vec_fluo(1),1,'last');
        stop_i = find(time_index_interp>=t_vec(end)&time_index_interp>=t_vec_fluo(end),1);
        t_vec_interp = time_index_interp(start_i:stop_i);           

        % extract basic vectors
        pt_vec = spot_struct(i).rawNCProtein/1e5;
        fluo_vec = spot_struct(i).fluoInterp;        
        ap_vec = spot_struct(i).APPosNucleus;
        % update time vector            
%         spot_struct(i).timeInterpOrig = spot_struct(i).timeInterp;
        spot_struct(i).timeNew = t_vec_interp;
        
        if length(t_vec) >= minDP && ~any(isnan(pt_vec)) && sum(~isnan(t_vec_interp)) >= minDP/2
            spot_struct(i).useFlag = true;
            spot_struct(i).fluoFlag = false;
            raw_adjusted = pt_vec;
            pert_ind = find(t_vec>=shift_time,1);
%             pert_ind = max([pert_ind,pert_ind-1])
            if ~isempty(pert_ind) 
                raw_adjusted(pert_ind:end) = ...
                (raw_adjusted(pert_ind:end)-cal_intercept)/cal_slope;
            end
            pt_interp = interp1(t_vec,raw_adjusted,t_vec_interp,'linear','extrap');
            ap_interp = interp1(t_vec,ap_vec,t_vec_interp,'linear','extrap');
            spot_struct(i).rawNCProteinNew = pt_interp;
            spot_struct(i).fluoNew = NaN(size(pt_interp));
            spot_struct(i).apPosNucleusNew = ap_interp;
            
            if sum(~isnan(t_vec_fluo)) >= minDP/2
                spot_struct(i).fluoFlag = true;
                fluo_new = zeros(size(t_vec_interp));
                fluo_new(ismember(t_vec_interp,t_vec_fluo)) = fluo_vec;           
                spot_struct(i).fluoNew = fluo_new;
            end
        else
            spot_struct(i).useFlag = false;
            spot_struct(i).fluoFlag = false;
            spot_struct(i).rawNCProteinNew = NaN(size(t_vec_interp));
            spot_struct(i).apPosNucleusNew = NaN(size(t_vec_interp));
            spot_struct(i).fluoNew = NaN(size(t_vec_interp));
        end    
    end        

    % check that corrections look reasonable
    temp_array = NaN(length(time_index_interp),length(master_ids));    
    iter = 1;
    for i = master_ids
        t_vec = spot_struct(i).timeNew;
        if spot_struct(i).useFlag
            temp_array(ismember(time_index_interp,t_vec),iter) = spot_struct(i).rawNCProteinNew;
        end
        iter = iter + 1;
    end
    mean_protein_array(:,m) = nanmean(temp_array,2);    
end    


%% Chop traces up into ~15 minute pieces to try and extend the dynamic range
shift_inc = round(5*60 / spot_struct(1).tresInterp);
window_size = round(15*60 / spot_struct(1).tresInterp);
ap_bounds = [-0.15 0.15];

inf_ind_vec = find([spot_struct.fluoFlag]);
inference_data = struct;

% initialize key values vectors
inference_data.fluo_vec_cell = {};
inference_data.knirps_vec_cell = {};
inference_data.time_vec_cell = {};
inference_data.ap_vec_cell = {};
inference_data.mean_fluo_vec = [];
inference_data.mean_knirps_vec = [];
inference_data.mean_time_vec = [];
inference_data.mean_ap_vec = [];

% initialize key metadata vectors
inference_data.project_id_vec = [];
inference_data.particle_id_vec = [];
inference_data.rep_id_vec = [];
inference_data.zero_flag_vec = [];

iter = 1;
wb = waitbar(0,'building inference set');
for i = inf_ind_vec
   waitbar(iter/length(inf_ind_vec),wb);
   % meta data
   particleID = spot_struct(i).particleID;
   projectID = spot_struct(i).projectID;
   
   % values
   fluo_vec = spot_struct(i).fluoNew;
   first_i = find(fluo_vec>0,1)+6;
   fluo_vec = fluo_vec(first_i:end);
   time_vec = spot_struct(i).timeNew(first_i:end);   
   kni_vec = spot_struct(i).rawNCProteinNew(first_i:end);
   ap_vec = spot_struct(i).apPosNucleusNew(first_i:end);
   
   % calculat how many distinct shifts we have
   n_shifts = 1 + floor((length(time_vec)-window_size)/shift_inc);
   last_obs = find(fluo_vec>0,1,'last');
   for n = 1:n_shifts
      start_i = shift_inc*(n-1)+1;
      last_i = start_i + window_size-1;
      ap_temp = ap_vec(start_i:last_i);
      ap_mean = nanmean(ap_temp);
      
      if ap_mean >= ap_bounds(1) && ap_mean <= ap_bounds(2) && ~isempty(last_obs)
          % record values
          inference_data.fluo_vec_cell{iter} = fluo_vec(start_i:last_i);
          inference_data.knirps_vec_cell{iter} = kni_vec(start_i:last_i);
          inference_data.time_vec_cell{iter} = time_vec(start_i:last_i);
          inference_data.ap_vec_cell{iter} = ap_vec(start_i:last_i);
          inference_data.mean_fluo_vec(iter) = nanmean(fluo_vec(start_i:last_i));
          inference_data.mean_knirps_vec(iter) = nanmean(kni_vec(start_i:last_i));
          inference_data.mean_time_vec(iter) = nanmean(time_vec(start_i:last_i));
          inference_data.mean_ap_vec(iter) = ap_mean;
          % metadata
          inference_data.particle_id_vec(iter) = particleID;
          inference_data.project_id_vec(iter) = projectID;
          inference_data.rep_id_vec(iter) = n;
          inference_data.zero_flag_vec(iter) = last_i>last_obs;
          % increment
          iter = iter + 1;
      end
   end
end  
delete(wb)

save([savePath 'inference_data.mat'],'inference_data')
save([savePath 'spot_struct.mat'],'spot_struct')

% % now, conduct 2 state viterbi fits using generic parameters
% nWorkersMax = 24;
% 
% % get parameters to use
% fitParameters = struct;
% fitParameters = getMarkovSystemInfo(fitParameters);
% fitParameters.A = expm(fitParameters.R2_orig*fitParameters.deltaT);
% 
% A_log = log(fitParameters.A);
% v = double(fitParameters.r2);
% sigma = fitParameters.noise;
% pi0_log = log(fitParameters.pi0');
% nStates = size(A_log,1);
% nSteps = fitParameters.memory;
% alpha = fitParameters.t_MS2;
% 
% % obtain subset of valid traces                
% use_flags = [spot_struct.fluoFlag]==1;
% fit_trace_indices = find(use_flags);                             
% 
% fluo_values = cell(length(fit_trace_indices),1);                
% for i = 1:numel(fit_trace_indices)    
%     fluo_values{i} = spot_struct(fit_trace_indices(i)).fluoNew;    
% end    
% 
% % initialize pool 
% p = gcp('nocreate'); % If no pool, do not create new one.
% if isempty(p)
%     parpool(nWorkersMax);
% elseif p.NumWorkers~= nWorkersMax
%     delete(p);
%     parpool(nWorkersMax);
% end      
% 
% disp('conducting viterbi trace fits...')
% 
% v_fits = struct;
% parfor f = 1:length(fluo_values)
%     viterbi_out = viterbi (fluo_values{f}, v, sigma, pi0_log, A_log, nStates, nSteps, alpha);
%     fnames = fieldnames(viterbi_out);
%     for fn = 1:numel(fnames)
%         v_fits(f).(fnames{fn}) = viterbi_out.(fnames{fn});
%     end        
% end
% 
% % Add to main structure
% for s = 1:length(spot_struct)
%     fit_filter = ismember(fit_trace_indices,s);
%     if any(fit_filter)
%         ind = find(fit_filter);
%         spot_struct(s).fluo_viterbi = v_fits(ind).fluo_viterbi;
%         spot_struct(s).path_viterbi = v_fits(ind).z_viterbi;
%         spot_struct(s).init_viterbi = v_fits(ind).v_viterbi;
%         spot_struct(s).logL = v_fits(ind).logL;
%     else
%         spot_struct(s).fluo_viterbi = NaN(size(spot_struct(s).fluoNew));
%         spot_struct(s).path_viterbi = NaN(size(spot_struct(s).fluoNew));
%         spot_struct(s).init_viterbi = NaN(size(spot_struct(s).fluoNew));
%         spot_struct(s).logL = NaN;
%     end
% end    
% 
% % let's remove entries without viterbi fit for now
% spot_struct = spot_struct(fit_trace_indices);
% 
% % remove unimportant fields to save space
% save([savePath 'spot_struct.mat'],'spot_struct')
