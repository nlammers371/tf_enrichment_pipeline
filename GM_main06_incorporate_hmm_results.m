% DESCRIPTION
% Script to conduct HMM inference
%
% ARGUMENTS
% singleTraceFitsDir: master ID variable 
%
% wInf: memory used for inference (GM 3/8/21: This option seems to be
% deprecated)
%
% KInf: number of states used for inference (GM 3/8/21: This option seems to be
% deprecated)
%
% OPTIONS
% dropboxFolder: Path to data folder where you wish to save
%                pipeline-generated data sets and figures. If this
%                var is not specified, output will be saved one level
%                above git repo in the folder structure
% 
% savioFlag: if running on savio change number of cores
%
% inferenceModel: specifies an inference model to use 
%
% controlProject: specifies a project to use as an external control
%
% OUTPUT: hmm_input_output, structure containing vectors of protein and MS2
% intensities, along with corresponding HMM-decoded activity trajectories

function hmm_input_output = GM_main06_incorporate_hmm_results(singleTraceFitsDir,varargin)

close all
warning('off', 'all')
addpath(genpath('utilities'))

currentDir = pwd;
savioFlag = 0;
if contains(currentDir, 'global/')
    savioFlag = 1;
end


% process other options
for i = (1+~savioFlag):2:(numel(varargin)-1)  
    if i ~= numel(varargin)        
        eval([varargin{i} '=varargin{i+1};']);                
    end    
end

load([singleTraceFitsDir 'singleTraceFitInfo.mat'],'singleTraceFitInfo')


projectName = singleTraceFitInfo.projectNameCell;
inferenceModel = singleTraceFitInfo.inferenceModel;
skipViterbiFitFlag = singleTraceFitInfo.skipViterbiFitFlag;
skipSSFitFlag =  singleTraceFitInfo.skipViterbiFitFlag;
makeLongFormSet = singleTraceFitInfo.makeLongFormSet;

myCluster = parcluster('local');
if ~savioFlag
    nWorkersMax = ceil(myCluster.NumWorkers/2);
else
    nWorkersMax = myCluster.NumWorkers; 
end
    
bootstrap_flag = singleTraceFitInfo.bootstrap_flag;
n_bootstraps = singleTraceFitInfo.n_bootstraps;




% get path to results
if ~exist('resultsRoot','var') % GM made some edits to allow for using subdirectories within project folders
    if contains(projectName(end), '/') | contains(projectName(end), '\')
        projectName = projectName(1:end-1);
    end
 
    [~, resultsRoot] = getDataPaths(true, projectName);
    resultsDir = [resultsRoot 'cpHMM_results' filesep];

  
    
else
    resultsRoot = [resultsRoot filesep projectName filesep];
    resultsDir = [resultsRoot filesep 'cpHMM_results' filesep];
end

% get list of all inference subdirectories. By default, we'll generate
% summaries for all non-empty inference sub-directory
if isempty(inferenceModel)
    infDirList = dir([resultsDir 'w*']); % get list of all inference subdirectories. By default, we'll generate
else
    infDirList = dir([resultsDir inferenceModel '*']); 
end
    
if length(length(infDirList)) > 1
    warning('Multiple inference directories found. This script will currently only process the first directory in the list')
end

for inf_i = 1:length(infDirList)
  
    % get path to results
    resultsPath = [infDirList(inf_i).folder filesep];
    
    % load inference options
    load([resultsPath infDirList(inf_i).name filesep 'inferenceOptions.mat'])
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
    load([resultsPath 'compiledResults_' infDirList(inf_i).name '.mat'],'compiledResults');
    
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
    trace_fit_flag_viterbi = 1;
    inf_files = dir([resultsPath infDirList(inf_i).name filesep 'hmm_results*.mat']);
    if exist([resultsPath 'singleTraceFits_' infDirList(inf_i).name '.mat']) > 0
        fit_props = dir([resultsPath 'singleTraceFits_' infDirList(inf_i).name '.mat']);
        fit_date = datenum(fit_props(1).date);
        hmm_date = datenum(inf_files(1).date);
        if fit_date > hmm_date
            trace_fit_flag_viterbi = 0;
        end
    end
    
    % get list of particle IDs
    trace_particle_index = [analysis_traces.particleID];
    
    % perform viterbi trace decoding if necessary
    if trace_fit_flag_viterbi
        singleTraceFits = performSingleTraceFits(compiledResults, inferenceOptions, bootstrap_flag, ...
                  analysis_traces, trace_particle_index, nWorkersMax, resultsPath, infDirList(inf_i).name, skipViterbiFitFlag);        
    else
        load([resultsPath 'singleTraceFits_' infDirList(inf_i).name '.mat'],'singleTraceFits')
    end
    
    % perform "soft" posterior trace decoding if necessary
    % check for existence of fit structure
    trace_fit_flag_ss = 1;
    inf_files = dir([resultsPath infDirList(inf_i).name filesep 'hmm_results*.mat']);
    if exist([resultsPath 'singleTraceFitsSS_' infDirList(inf_i).name '.mat']) > 0
        fit_props = dir([resultsPath 'singleTraceFitsSS_' infDirList(inf_i).name '.mat']);
        fit_date = datenum(fit_props(1).date);
        hmm_date = datenum(inf_files(1).date);
        if fit_date > hmm_date
            trace_fit_flag_ss = 0;
        end
    end
    
    if ~skipSSFitFlag
        if trace_fit_flag_ss
            singleTraceFitsSS = performSingleTraceFitsSoft(compiledResults, inferenceOptions, bootstrap_flag, ...
                      analysis_traces, trace_particle_index, nWorkersMax, resultsPath, infDirList(inf_i).name);        
        else
            load([resultsPath 'singleTraceFitsSS_' infDirList(inf_i).name '.mat'],'singleTraceFitsSS')
        end
    end
    % generate longform dataset
    if makeLongFormSet        
        resultsTable = generateLongFormTable(analysis_traces, timeGrid, singleTraceFits, resultsPath, infDirList(inf_i).name);
    end

end



% next, incorporate local enrichment fields (if they exist)
if exist([resultsRoot filesep 'spot_struct_protein.mat'],'file')
    %%% now combine with protein data and raw traces
    disp('building input/output dataset...')

    % first combine HMM and raw fluoresence trace data points
    hmm_input_output = addTraceAnalysisFields(analysis_traces,singleTraceFits,singleTraceFitsSS,inferenceOptions);

    load([resultsRoot filesep 'spot_struct_protein.mat'],'spot_struct_protein')
    fitParticleVec = [singleTraceFits.particleID];
    traceParticleVec = [spot_struct_protein.particleID];
    
    for trace_i = 1:length(hmm_input_output)
          
        % basic indexing info
        ParticleID = hmm_input_output(trace_i).particleID(1);   
        fitIndex = find(fitParticleVec==ParticleID);
        if isempty(fitIndex)
            fitIndex = NaN;
        end
        traceIndex = find(traceParticleVec==ParticleID);
  
        % extract time vectors
        timeFrom = spot_struct_protein(traceIndex).time;
        timeTo = hmm_input_output(trace_i).time;
        
        % interpolate protein vectors to match cpHMM spacing
        interpFields = {'spot_protein_vec','serial_null_protein_vec', 'edge_null_protein_vec', 'nuclear_protein_vec',...
          'spot_mcp_vec'};
        
        for field_i = 1:length(interpFields)
            vecRaw = spot_struct_protein(traceIndex).(interpFields{field_i});
            nanFT = ~isnan(vecRaw);
            if sum(nanFT) > 1
                hmm_input_output(trace_i).(interpFields{field_i}) = interp1(timeFrom(nanFT), vecRaw(nanFT), timeTo);
            elseif sum(nanFT) == 1
                vecTo = NaN(size(timeTo));
                [~,mi] = min(abs(timeTo-timeFrom(nanFT)));
                vecTo(mi) = vecRaw(nanFT);
                hmm_input_output(trace_i).(interpFields{field_i}) = vecTo;
            else
                hmm_input_output(trace_i).(interpFields{field_i}) = NaN(size(timeTo));
            end
        end
                                             
        % generate flag var indicating interpolated obs that are too far from 
        % true points
        input_times = timeFrom(~isnan(spot_struct_protein(trace_i).spot_protein_vec));
        dt_vec_gap = NaN(size(hmm_input_output(trace_i).time));
        if ~isempty(input_times)
            for t = 1:numel(dt_vec_gap)
                dt_vec_gap(t) = min(abs(hmm_input_output(trace_i).time(t)-input_times));   
            end
        end
        hmm_input_output(trace_i).dt_filter_gap = dt_vec_gap > maxDT;            
        % record general info for later use        
        hmm_input_output(trace_i).Tres = Tres;
        hmm_input_output(trace_i).maxDT = maxDT;
        hmm_input_output(trace_i).InferenceID = fitIndex;            
        hmm_input_output(trace_i).rawSpotStructIndex = traceIndex;            
        
    end
    
    % remove traces with not cpHMM fit info 
    hmm_input_output = hmm_input_output(~isnan([hmm_input_output.InferenceID]));

    % save results
    save([resultsDir 'hmm_input_output.mat'],'hmm_input_output')
end  


