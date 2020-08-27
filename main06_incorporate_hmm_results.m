% DESCRIPTION
% Script to conduct HMM inference
%
% ARGUMENTS
% project: master ID variable 
%
% wInf: memory used for inference
%
% KInf: number of states used for inference
%
% OPTIONS
% dropboxFolder: Path to data folder where you wish to save
%                pipeline-generated data sets and figures. If this
%                var is not specified, output will be saved one level
%                above git repo in the folder structure
%
% controlProject: specifies a project to use as an external control
%
% OUTPUT: hmm_input_output, structure containing vectors of protein and MS2
% intensities, along with corresponding HMM-decoded activity trajectories

function hmm_input_output = main06_incorporate_hmm_results(project,DropboxFolder,varargin)

close all
addpath('./utilities')
%%%%% These options will remain fixed for now
alphaFrac = 1302 / 6444;
if contains(project,'hbP2P')
    alphaFrac = 1275 / 4670;
elseif contains(project,'snaBAC')
    alphaFrac = 1302 / 6444;
end
[~, DataPath, ~] =   header_function(DropboxFolder, project);
w = 7;
K = 3;
nWorkers = 24;
fluo_dim = 2;
protein_dim = fluo_dim;
n_boots_max = 1; % max num protein bootstraps to use for soft fits

%%%%%%%%%%%%%%
for i = 1:numel(varargin)       
    if ischar(varargin{i}) && i < numel(varargin)        
        eval([varargin{i} '=varargin{i+1};']);        
    end
end
% load master nucleus data set
load([DataPath '/nucleus_struct_protein.mat'],'nucleus_struct_protein') % load data
load([DataPath '/nucleus_struct.mat'],'nucleus_struct') % load data
minDP = 1;%nucleus_struct_protein(1).minDP;

Tres = nucleus_struct_protein(1).TresInterp; % Time Resolution
maxDT = 1.2*Tres; % maximum distance from observed data point
alpha = alphaFrac*w;
% generate alpha kernel for estimating predicted HMM fluroescence
alpha_kernel = NaN(1,w);
for i = 1:w
    if i < alpha
        alpha_kernel(i) = ((i-1) / alpha  + .5 * 1/alpha )*Tres;
    elseif i > alpha && (i-1) < alpha
        alpha_kernel(i) = Tres*(1 - .5*(alpha-i+1)*(1-(i-1)/alpha));
    else
        alpha_kernel(i) = Tres;
    end
end
% Set write path (inference results are now written to external directory)
if contains(project,'snaBAC')
    hmm_suffix =  ['hmm_inference_protein/w' num2str(w) '_K' num2str(K) '_f' num2str(fluo_dim) 'D/']; 
else
    hmm_suffix =  ['hmm_inference_mf/w' num2str(w) '_K' num2str(K) '_f' num2str(fluo_dim) 'D/']; 
end
file_list = dir([DataPath hmm_suffix 'hmm_results*.mat']);

%%%%%%%%%%%%%%%%%
% load in inference results
%%%%%%%%%%%%%%%%%
inference_results = struct;
n_boots = min([numel(file_list) n_boots_max]);
iter = 1;
for inf = 1:n_boots
    load([DataPath hmm_suffix file_list(inf).name]);    
    fnames = fieldnames(output);
    if numel(fnames) > 1
        for fn = 1:numel(fnames)
            inference_results(iter).(fnames{fn}) = output.(fnames{fn});
        end
        inference_results(iter).source = [DataPath hmm_suffix file_list(inf).name];
        iter = iter + 1;
    end
end

% check for existence of fit structure
trace_fit_flag = 1;
if exist([DataPath hmm_suffix 'trace_fit_struct.mat']) > 0
    fit_props = dir([DataPath hmm_suffix 'trace_fit_struct.mat']);
    fit_date = datenum(fit_props(1).date);
    hmm_date = datenum(file_list(1).date);
    if fit_date > hmm_date
        trace_fit_flag = 0;
    end
end

qc_indices = find([nucleus_struct_protein.qc_flag]==1);
particle_index = [nucleus_struct_protein.ParticleID];
% qc_indices = qc_indices(1:10);
% perform soft trace decoding if necessary
if trace_fit_flag            
    trace_fit_struct = struct; 
    mf_trace_indices = qc_indices;
    for inf = 1:numel(inference_results)
        
        A_log = log(inference_results(inf).A_mat);
        v = inference_results(inf).r*Tres;
        sigma = sqrt(inference_results(inf).noise);
        pi0_log = log(inference_results(inf).pi0); 
        eps = 1e-4;
        
        % deduce subset of valid traces
        
        fluo_values = cell(numel(qc_indices),1);
        for i = 1:numel(qc_indices)
           if fluo_dim == 2
                fluo = nucleus_struct_protein(mf_trace_indices(i)).fluo_interp;
            else
                fluo = nucleus_struct_protein(mf_trace_indices(i)).fluo3D_interp;
            end
            start_i = find(~isnan(fluo),1);
            stop_i = find(~isnan(fluo),1,'last');
            fluo = fluo(start_i:stop_i);
            fluo_values{i} = fluo;
        end    
        % viterbi fitting        
    
        
        p = gcp('nocreate'); % If no pool, do not create new one.
        if isempty(p)
            parpool(nWorkers);
        elseif p.NumWorkers~= nWorkers
            delete(p);
            parpool(nWorkers);
        end        
        v_fits = struct;
        parfor f = 1:numel(fluo_values)
%             waitbar(f/numel(fluo_values),h);
            viterbi_out = viterbi (fluo_values{f}, v', sigma, pi0_log,A_log, K, w, alpha);
            fnames = fieldnames(viterbi_out);
            for j = 1:numel(fnames)
                v_fits(f).(fnames{j}) = viterbi_out.(fnames{j});
            end
%             f/numel(fluo_values)
        end
        
        % soft decoding
        disp('conducting soft-decoded trace fits...')
        local_em_outputs = local_em_MS2_reduced_memory (fluo_values, ...
                                v', sigma, pi0_log, A_log, K, w, alpha, 1, eps);              

        trace_fit_struct(inf).viterbi_fits = v_fits;
        
        trace_fit_struct(inf).soft_fits = local_em_outputs.soft_struct.p_z_log_soft;
        trace_fit_struct(inf).particle_index = particle_index(qc_indices);
        trace_fit_struct(inf).inference_id_vec = repelem(inf,numel(particle_index));
        
    end        
    save([DataPath hmm_suffix 'trace_fit_struct.mat'],'trace_fit_struct')
else
    load([DataPath hmm_suffix 'trace_fit_struct.mat'],'trace_fit_struct')
end

%%% now extract corresponding hmm traces
disp('building input/output dataset...')
hmm_input_output = [];
for inf = 1:numel(trace_fit_struct)
    iter = 1;
    soft_fits = trace_fit_struct(inf).soft_fits;
    viterbi_fits = trace_fit_struct(inf).viterbi_fits;
    inference_id_vec = trace_fit_struct(inf).inference_id_vec;
    for i = qc_indices
        % initialize temporary structure to store results
        ParticleID = particle_index(i);    
        temp = struct;
        % extract relevant vectors from protein struct    
        % these quantities have not been interpolated
        if fluo_dim == 3   
            ff_pt = nucleus_struct_protein(i).fluo3D;
            master_fluo = nucleus_struct_protein(i).fluo3D_interp;            
        elseif fluo_dim == 2
            ff_pt = nucleus_struct_protein(i).fluo;
            master_fluo = nucleus_struct_protein(i).fluo_interp;            
        end        
        if protein_dim == 3
            sp_pt = nucleus_struct_protein(i).spot_protein_vec_3d;
            sr_pt = nucleus_struct_protein(i).serial_null_protein_vec_3d; 
        elseif protein_dim == 2
            sp_pt = nucleus_struct_protein(i).spot_protein_vec;
            sr_pt = nucleus_struct_protein(i).serial_null_protein_vec; 
        end
        mcp_pt = nucleus_struct_protein(i).spot_mcp_vec;         
        nn_pt = nucleus_struct_protein(i).edge_null_protein_vec;
        mf_pt_mf = nucleus_struct_protein(i).mf_null_protein_vec;          
        tt_pt = nucleus_struct_protein(i).time;
        
        x_pt = nucleus_struct_protein(i).xPosParticle;  
        y_pt = nucleus_struct_protein(i).yPosParticle;  
        ap_pt = nucleus_struct_protein(i).APPosParticle;  
        
        if sum(~isnan(mf_pt_mf)) > minDP && sum(~isnan(sr_pt)) > minDP && sum(~isnan(sp_pt)) > minDP            
            % extract interpolated fluorescence and time vectors
            master_time = nucleus_struct_protein(i).time_interp;
            
            % extract position vectors (used for selecting nearest neighbor)
            x_nc = double(nucleus_struct_protein(i).xPos);
            y_nc = double(nucleus_struct_protein(i).yPos);
            temp.xPosMean = nanmean(x_nc(~isnan(ff_pt)));
            temp.yPosMean = nanmean(y_nc(~isnan(ff_pt)));
            
            % record time, space, and fluo vars
            start_i = find(~isnan(master_fluo),1);
            stop_i = find(~isnan(master_fluo),1,'last');
            temp.time = master_time(start_i:stop_i);
            temp.fluo = master_fluo(start_i:stop_i);   
            temp.fluo_raw = ff_pt;      

            % extract useful HMM inference parameters             
            inference_id = inference_id_vec(iter); % inference id
            [r,ri] = sort(inference_results(inference_id).r); % enforce rank ordering of states
            z = exp(soft_fits{iter});    
            temp.z_mat = z(ri,:)';    
            temp.r_mat = z(ri,:)'.*r';
            temp.r_inf = r';
            temp.r_vec = sum(temp.r_mat,2)';
            [~,z_vec] = max(temp.z_mat,[],2);
            temp.z_vec = z_vec; 
            
            % extract viterbi fits
            temp.z_viterbi = viterbi_fits(iter).z_viterbi;
            temp.f_viterbi = viterbi_fits(iter).fluo_viterbi;
            
            % make predicted fluo vec (for consistency checks)
            fluo_hmm = conv(temp.r_vec,alpha_kernel);
            temp.fluo_hmm = fluo_hmm(1:numel(temp.r_vec));        

            % checks using mcp channel to ensure that we are correctly matching
            % particles and time frames
            temp.mcp_check = interp1(tt_pt(~isnan(mcp_pt)),mcp_pt(~isnan(mcp_pt)),temp.time);
            temp.fluo_check = interp1(tt_pt(~isnan(ff_pt)),ff_pt(~isnan(ff_pt)),temp.time);
            
            % record raw data vectors
            temp.spot_protein_raw = sp_pt;        
            temp.mf_protein_raw = mf_pt_mf;
            temp.null_protein_raw = nn_pt;
            temp.serial_protein_raw = sr_pt;  
            temp.time_raw = tt_pt;  
            temp.xPos_raw = x_nc;
            temp.yPos_raw = y_nc;
            
            % interpolate protein information such that it coincides with HMM
            % inference results    
            temp.spot_protein = interp1(tt_pt(~isnan(sp_pt)),sp_pt(~isnan(sp_pt)),temp.time);                    
            temp.serial_protein = interp1(tt_pt(~isnan(sr_pt)),sr_pt(~isnan(sr_pt)),temp.time);            
            temp.null_protein = interp1(tt_pt(~isnan(nn_pt)),nn_pt(~isnan(nn_pt)),temp.time);
            temp.mf_protein = interp1(tt_pt(~isnan(mf_pt_mf)),mf_pt_mf(~isnan(mf_pt_mf)),temp.time);
            
            % interpolate position info
            temp.xPos = interp1(tt_pt(~isnan(x_nc)),x_nc(~isnan(x_nc)),temp.time);        
            temp.yPos = interp1(tt_pt(~isnan(y_nc)),y_nc(~isnan(y_nc)),temp.time);
            
            temp.xPosParticle = interp1(tt_pt(~isnan(x_pt)),x_pt(~isnan(x_pt)),temp.time);        
            temp.yPosParticle = interp1(tt_pt(~isnan(y_pt)),y_pt(~isnan(y_pt)),temp.time);
            temp.apPosParticle = interp1(tt_pt(~isnan(ap_pt)),ap_pt(~isnan(ap_pt)),temp.time);
            % generate flag var indicating interpolated obs that are too far from 
            % true points
            input_times = tt_pt(~isnan(sp_pt));
            dt_vec_gap = NaN(size(temp.time));
            for t = 1:numel(dt_vec_gap)
                dt_vec_gap(t) = min(abs(temp.time(t)-input_times));   
            end
            temp.dt_filter_gap = dt_vec_gap > maxDT;            
            % record general info for later use
            temp.ParticleID = ParticleID; 
            temp.Tres = Tres;
            temp.maxDT = maxDT;
            temp.InferenceID = inf;    
            hmm_input_output  = [hmm_input_output temp];
        end
        % increment
        iter = iter + 1;
    end
end

% find nearest neighbor particles
% use name nearest neighbor for each bootstrap instance
n_unique = numel(hmm_input_output) / n_boots;%numel(inference_results);
start_time_vec = NaN(1,n_unique);
stop_time_vec = NaN(1,n_unique);
set_vec = NaN(1,n_unique);
for i = 1:n_unique
    dt_flag = hmm_input_output(i).dt_filter_gap;
    t_vec = hmm_input_output(i).time(~dt_flag);
    start_time_vec(i) = min(t_vec);
    stop_time_vec(i) = max(t_vec);
    set_vec(i) = floor(hmm_input_output(i).ParticleID);
end

% xy nearest neighbor calculations
dist_mat_x = pdist2([hmm_input_output(1:n_unique).xPosMean]',[hmm_input_output(1:n_unique).xPosMean]');
dist_mat_y = pdist2([hmm_input_output(1:n_unique).yPosMean]',[hmm_input_output(1:n_unique).yPosMean]');
dist_mat_r = sqrt(dist_mat_x.^2 + dist_mat_y.^2);

% now find closest match for each nucleus
for i = 1:n_unique
    % require that mat trace is (1) from same set, (2) starts and ends
    % within 3 min of locus trace
    setID = set_vec(i);  
    option_filter = ((start_time_vec-start_time_vec(i)) <= 3*60) & ...
        ((stop_time_vec-stop_time_vec(i)) >= -3*60) & set_vec==setID;        
    
    %%% Spatial Nearest Neighbor   
    time_i = hmm_input_output(i).time;
    dist_vec = dist_mat_r(i,:);               
    dist_vec(~option_filter) = NaN;
    dist_vec(i) = NaN; % remove self
    [best_r, best_ind_dist] = nanmin(dist_vec);
    % record vales 
    time_swap_dist = hmm_input_output(best_ind_dist).time;       
    % fill
    swap_ft = ismember(time_swap_dist,time_i);
    target_ft = ismember(time_i,time_swap_dist);
    s_pt_dist = NaN(size(time_i));
    s_pt_dist(target_ft) = hmm_input_output(best_ind_dist).spot_protein(swap_ft);    
    mf_pt_dist = NaN(size(time_i));
    mf_pt_dist(target_ft) = hmm_input_output(best_ind_dist).mf_protein(swap_ft);    
    r_fluo_dist = NaN(size(time_i));
    r_fluo_dist(target_ft) = hmm_input_output(best_ind_dist).r_vec(swap_ft);
    dt_filter_dist = true(size(time_i));
    dt_filter_dist(target_ft) = hmm_input_output(best_ind_dist).dt_filter_gap(swap_ft);
    
    % assign to ALL copies
    for inf = 1:numel(trace_fit_struct)
        ind = (inf-1)*n_unique + i;
        % record dist
        hmm_input_output(ind).nn_best_r = best_r;
        hmm_input_output(ind).dist_swap_ind = best_ind_dist;
        hmm_input_output(ind).dist_swap_spot_protein = s_pt_dist;
        hmm_input_output(ind).dist_swap_mf_protein = mf_pt_dist;
        hmm_input_output(ind).dist_swap_hmm = r_fluo_dist;
        hmm_input_output(ind).dist_swap_dt_filter_gap = dt_filter_dist;
    end
end

% save results
save([DataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '_f' num2str(fluo_dim)  'D_p' num2str(protein_dim) 'D.mat'],'hmm_input_output')