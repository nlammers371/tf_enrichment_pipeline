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
   
for trace_i = 1:length(infDirList)
  
    % get path to results
    resultsPath = [infDirList(trace_i).folder filesep];
    
    % load inference options
    load([resultsPath infDirList(trace_i).name filesep 'inferenceOptions.mat'])
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
    load([resultsPath 'compiledResults_' infDirList(trace_i).name '.mat'],'compiledResults');
    
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
    inf_files = dir([resultsPath infDirList(trace_i).name filesep 'hmm_results*.mat']);
    if exist([resultsPath 'singleTraceFits_' infDirList(trace_i).name '.mat']) > 0
        fit_props = dir([resultsPath 'singleTraceFits_' infDirList(trace_i).name '.mat']);
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
        load([resultsPath 'singleTraceFits_' infDirList(trace_i).name '.mat'],'singleTraceFits')
    end
  
    % generate longform dataset
    if makeLongFormSet        
        resultsTable = generateLongFormTable(analysis_traces, timeGrid, singleTraceFits, resultsPath, infDirList(trace_i).name);
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

% first combine HMM and raw fluoresence trace data points
hmm_input_output = addTraceAnalysisFields(analysis_traces,singleTraceFits);

% next, incorporate local enrichment fields (if they exist)
if exist([resultsRoot filesep 'spot_struct_protein.mat'],'file')
    load([resultsRoot filesep 'spot_struct_protein.mat'],'spot_protein_struct')
    fitParticleVec = [singleTraceFits.particleID];
    traceParticleVec = [spot_struct_protein.particleID];
    
    for trace_i = 1:length(hmm_input_output)
          
        % basic indexing info
        ParticleID = hmm_input_output(trace_i);   
        fitIndex = find(fitParticleVec==ParticleID);
        traceIndex = find(traceParticleVec==ParticleID);
  
        % extract time vectors
        timeFrom = analysis_traces(traceIndex).time;
        timeTo = hmm_input_output(traceIndex).time;
        
        % interpolate protein vectors to match cpHMM spacing
        interpFields = {'spot_protein_vec','serial_null_protein_vec', 'edge_null_protein_vec', 'nuclear_protein_vec',...
          'spot_mcp_vec'};
        for field_i = 1:length(interpFields)
            vecRaw = spot_struct_protein(traceIndex).(interpFields{field_i});
            nanFT = ~isnan(vec);
            if sum(nanFT) > 1
                hmm_input_output(trace_i).(interpFields{field_i}) = interp1(timeFrom(nanFT), vecRaw(nanFT), timeTo);
            elseif sum(nanFT) == 1
                vecTo = NaN(size(timeTo));
                [~,mi] = min(abs(timeTo-timeFrom(nanFT)));
                vecTo(mi) = vec(nanFT);
                hmm_input_output(trace_i).(interpFields{field_i}) = vecTo;
            else
                hmm_input_output(trace_i).(interpFields{field_i}) = NaN(size(timeTo));
            end
        end
                                             
        % generate flag var indicating interpolated obs that are too far from 
        % true points
        input_times = timeFrom(~isnan(hmm_input_output(trace_i).spot_protein_vec));
        dt_vec_gap = NaN(size(hmm_input_output(trace_i).time));
        for t = 1:numel(dt_vec_gap)
            dt_vec_gap(t) = min(abs(hmm_input_output(trace_i).time(t)-input_times));   
        end
        hmm_input_output(trace_i).dt_filter_gap = dt_vec_gap > maxDT;            
        % record general info for later use        
        hmm_input_output(trace_i).Tres = Tres;
        hmm_input_output(trace_i).maxDT = maxDT;
        hmm_input_output(trace_i).InferenceID = fitIndex;            
        hmm_input_output(trace_i).rawSpotStructIndex = traceIndex;            
        
    end

    % find nearest neighbor particles
    % use name nearest neighbor for each bootstrap instance
%     n_unique = numel(hmm_input_output) / n_boots;%numel(inference_results);
%     start_time_vec = NaN(1,n_unique);
%     stop_time_vec = NaN(1,n_unique);
%     set_vec = NaN(1,n_unique);
%     for i = 1:n_unique
%         dt_flag = hmm_input_output(i).dt_filter_gap;
%         t_vec = hmm_input_output(i).time(~dt_flag);
%         start_time_vec(i) = min(t_vec);
%         stop_time_vec(i) = max(t_vec);
%         set_vec(i) = floor(hmm_input_output(i).ParticleID);
%     end
% 
%     % xy nearest neighbor calculations
%     dist_mat_x = pdist2([hmm_input_output(1:n_unique).xPosMean]',[hmm_input_output(1:n_unique).xPosMean]');
%     dist_mat_y = pdist2([hmm_input_output(1:n_unique).yPosMean]',[hmm_input_output(1:n_unique).yPosMean]');
%     dist_mat_r = sqrt(dist_mat_x.^2 + dist_mat_y.^2);
% 
%     % now find closest match for each nucleus
%     for i = 1:n_unique
%         % require that mat trace is (1) from same set, (2) starts and ends
%         % within 3 min of locus trace
%         setID = set_vec(i);  
%         option_filter = ((start_time_vec-start_time_vec(i)) <= 3*60) & ...
%             ((stop_time_vec-stop_time_vec(i)) >= -3*60) & set_vec==setID;        
% 
%         %%% Spatial Nearest Neighbor   
%         time_i = hmm_input_output(i).time;
%         dist_vec = dist_mat_r(i,:);               
%         dist_vec(~option_filter) = NaN;
%         dist_vec(i) = NaN; % remove self
%         [best_r, best_ind_dist] = nanmin(dist_vec);
%         % record vales 
%         time_swap_dist = hmm_input_output(best_ind_dist).time;       
%         % fill
%         swap_ft = ismember(time_swap_dist,time_i);
%         target_ft = ismember(time_i,time_swap_dist);
%         s_pt_dist = NaN(size(time_i));
%         s_pt_dist(target_ft) = hmm_input_output(best_ind_dist).spot_protein(swap_ft);    
%         mf_pt_dist = NaN(size(time_i));
%         mf_pt_dist(target_ft) = hmm_input_output(best_ind_dist).mf_protein(swap_ft);    
%         r_fluo_dist = NaN(size(time_i));
%         r_fluo_dist(target_ft) = hmm_input_output(best_ind_dist).r_vec(swap_ft);
%         dt_filter_dist = true(size(time_i));
%         dt_filter_dist(target_ft) = hmm_input_output(best_ind_dist).dt_filter_gap(swap_ft);
% 
%         % assign to ALL copies
%         for ind = i:n_unique:length(hmm_input_output)
%     %         ind = (inf-1)*n_unique + i;
%             % record dist
%             hmm_input_output(ind).nn_best_r = best_r;
%             hmm_input_output(ind).dist_swap_ind = best_ind_dist;
%             hmm_input_output(ind).dist_swap_spot_protein = s_pt_dist;
%             hmm_input_output(ind).dist_swap_mf_protein = mf_pt_dist;
%             hmm_input_output(ind).dist_swap_hmm = r_fluo_dist;
%             hmm_input_output(ind).dist_swap_dt_filter_gap = dt_filter_dist;
%         end
%     end
end
% save results
save([DataPath 'hmm_input_output_w' num2str(w) '_K' num2str(nSteps) '_f' num2str(fluo_dim)  'D_p' num2str(protein_dim) 'D.mat'],'hmm_input_output')
