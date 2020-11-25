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

function main06_incorporate_hmm_results(projectName,varargin)

close all
addpath(genpath('utilities'))

makeLongFormSet = 0;
myCluster = parcluster('local');
nWorkersMax = ceil(myCluster.NumWorkers/2);
bootstrap_flag = 0;
n_bootstraps = 1;

% check for optional inputs
for i = 1:(numel(varargin)-1)  
    if i ~= numel(varargin)        
        eval([varargin{i} '=varargin{i+1};']);                
    end    
end

% get path to results
if ~exist('resultsRoot','var')
    liveProject = LiveProject(projectName);
    resultsRoot = [liveProject.dataPath filesep];
    resultsDir = [resultsRoot 'cpHMM_results' filesep];
else
    resultsRoot = [resultsRoot filesep projectName filesep];
    resultsDir = [resultsRoot filesep 'cpHMM_results' filesep];
end
% get list of all inference subdirectories. By default, we'll generate
% summaries for all non-empty inference sub-directory
infDirList = dir([resultsDir 'w*']);% get list of all inference subdirectories. By default, we'll generate
   
for traceInd = 1:length(infDirList)
  
    % get path to results
    resultsPath = [infDirList(traceInd).folder filesep];
    
    % load inference options
    load([resultsPath infDirList(traceInd).name filesep 'inferenceOptions.mat'])
    if length(fieldnames(inferenceOptions))==1
        inferenceOptions = inferenceOptions.inferenceOptions;
    end
    
    % extract key inference parameters     
    nStates = inferenceOptions.nStates;
    nSteps = inferenceOptions.nSteps;
    alpha = inferenceOptions.alpha;
    Tres = inferenceOptions.Tres;
    timeGrid = 0:Tres:60*60;
    maxDT = 1.2*Tres; % maximum distance from observed data point
    
    % load compiled inference files
    load([resultsPath 'compiledResults_' infDirList(traceInd).name '.mat'],'compiledResults');
    
    % load corresponding trace structure
    if ~inferenceOptions.ProteinBinFlag
        load([resultsRoot filesep 'spot_struct.mat'])
        analysis_traces = spot_struct;
        clear spot_struct;
    else
        load([resultsRoot filesep 'spot_struct_protein.mat'])
        analysis_traces = spot_struct_protein;
        clear spot_struct_protein;
    end
               

    % check for existence of fit structure
    trace_fit_flag = 1;
    inf_files = dir([resultsPath infDirList(traceInd).name filesep 'hmm_results*.mat']);
    if exist([resultsPath 'singleTraceFits_' infDirList(traceInd).name '.mat']) > 0
        fit_props = dir([resultsPath 'singleTraceFits_' infDirList(traceInd).name '.mat']);
        fit_date = datenum(fit_props(1).date);
        hmm_date = datenum(inf_files(1).date);
        if fit_date > hmm_date
            trace_fit_flag = 0;
        end
    end
    
    % get list of particle IDs
    trace_particle_index = [analysis_traces.particleID];
    
    % perform viterbi trace decoding if necessary
    if trace_fit_flag    
        singleTraceFits = performSingleTraceFits(compiledResults, inferenceOptions, bootstrap_flag, ...
                  analysis_traces, trace_particle_index, nWorkersMax, resultsPath);        
    else
        load([resultsPath 'singleTraceFits_' infDirList(traceInd).name '.mat'],'singleTraceFits')
    end
  
    % generate longform dataset
    if makeLongFormSet        
        resultsTable = generateLongFormTable(analysis_traces, timeGrid, singleTraceFits, resultsPath, infDirList(traceInd).name);
    end

end

% generate alpha kernel for estimating predicted HMM fluroescence
alpha_kernel = NaN(1,nSteps);
for i = 1:nSteps
    if i < alpha
        alpha_kernel(i) = ((i-1) / alpha  + .5 * 1/alpha )*Tres;
    elseif i > alpha && (i-1) < alpha
        alpha_kernel(i) = Tres*(1 - .5*(alpha-i+1)*(1-(i-1)/alpha));
    else
        alpha_kernel(i) = Tres;
    end
end

%%% now combine with protein data and raw traces
disp('building input/output dataset...')
hmm_input_output = [];
for traceInd = 1:length(singleTraceFits)
    
    viterbi_fits = singleTraceFits(traceInd).viterbi_fits;
    inference_id_vec = singleTraceFits(traceInd).inference_id_vec;
    inference_particles = singleTraceFits(traceInd).particle_index;
    
    fit_trace_indices = find(ismember(trace_particle_index,inference_particles));
    for i = fit_trace_indices
        % initialize temporary structure to store results
        ParticleID = trace_particle_index(i);   
        traceIndex = find(inference_particles==ParticleID);
        if isempty(traceIndex)
          error('uh oh')
        end
        temp = struct;
        % extract relevant vectors from protein struct    
        % these quantities have not been interpolated
        if fluo_dim == 3   
            ff_pt = analysis_traces(i).fluo3D;
            master_fluo = analysis_traces(i).fluo3D_interp;            
        elseif fluo_dim == 2
            ff_pt = analysis_traces(i).fluo;
            master_fluo = analysis_traces(i).fluo_interp;            
        end        
        if protein_dim == 3
            sp_pt = analysis_traces(i).spot_protein_vec_3d;
            sr_pt = analysis_traces(i).serial_null_protein_vec_3d; 
        elseif protein_dim == 2
            sp_pt = analysis_traces(i).spot_protein_vec;
            sr_pt = analysis_traces(i).serial_null_protein_vec; 
        end
        mcp_pt = analysis_traces(i).spot_mcp_vec;         
        nn_pt = analysis_traces(i).edge_null_protein_vec;
        mf_pt_mf = analysis_traces(i).mf_null_protein_vec;          
        tt_pt = analysis_traces(i).time;

        x_pt = analysis_traces(i).xPosParticle;  
        y_pt = analysis_traces(i).yPosParticle;  
        ap_pt = analysis_traces(i).APPosParticle;  


        if sum(~isnan(mf_pt_mf)) > minDP && sum(~isnan(sr_pt)) > minDP && sum(~isnan(sp_pt)) > minDP              
            % extract interpolated fluorescence and time vectors
            master_time = analysis_traces(i).time_interp;

            % check for mismatch between analysis_traces and
            % analysis_traces...this is due to a dumb mistake on my
            % part
            time_vec_orig = analysis_traces(trace_particle_index_orig==ParticleID).time_interp;
            if ~isequal(master_time,time_vec_orig)
              master_time = time_vec_orig;
              analysis_traces(i).time_interp = master_time;

              master_fluo = analysis_traces(trace_particle_index_orig==ParticleID).fluo_interp;
              analysis_traces(i).fluo_interp = master_fluo;
            end
            % extract position vectors (used for selecting nearest neighbor)
            x_nc = double(analysis_traces(i).xPos);
            y_nc = double(analysis_traces(i).yPos);
            temp.xPosMean = nanmean(x_nc(~isnan(ff_pt)));
            temp.yPosMean = nanmean(y_nc(~isnan(ff_pt)));

            % record time, space, and fluo vars
            start_i = find(~isnan(master_fluo),1);
            stop_i = find(~isnan(master_fluo),1,'last');
            temp.time = master_time(start_i:stop_i);
            temp.fluo = master_fluo(start_i:stop_i);   
            temp.fluo_raw = ff_pt;      

            % extract useful HMM inference parameters             
            inference_id = inference_id_vec(traceIndex); % inference id
            [r,ri] = sort(inference_results(inference_id).r); % enforce rank ordering of states
            z = exp(soft_fits{traceIndex});    
            temp.z_mat = z(ri,:)';    
            temp.r_mat = z(ri,:)'.*r';
            temp.r_inf = r';
            temp.r_vec = sum(temp.r_mat,2)';
            [~,z_vec] = max(temp.z_mat,[],2);
            temp.z_vec = z_vec; 
            if length(z_vec)~=length(temp.fluo)
              error('goddammit')
            end
            % extract viterbi fits
            temp.z_viterbi = viterbi_fits(traceIndex).z_viterbi;
            temp.f_viterbi = viterbi_fits(traceIndex).fluo_viterbi;

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
            temp.InferenceID = traceInd;    
            hmm_input_output  = [hmm_input_output temp];
        end
        % increment
%         iter = iter + 1;
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
    for ind = i:n_unique:length(hmm_input_output)
%         ind = (inf-1)*n_unique + i;
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
save([DataPath 'hmm_input_output_w' num2str(w) '_K' num2str(nSteps) '_f' num2str(fluo_dim)  'D_p' num2str(protein_dim) 'D.mat'],'hmm_input_output')
end