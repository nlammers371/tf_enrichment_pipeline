% Script to scan through parameter space for selected biological variables
function sweepInfo = io_sweep_wrapper(resultsRoot,nParamIncrement,simType,...
                                      paramVals,keep_prediction_flag,varargin)

granularity = 1;
n_traces = 100;
for i = 1:numel(varargin)
   if ischar(varargin{i}) && i <  numel(varargin)
       eval([varargin{i} ' = varargin{i+1};'])
   end
end

% load data
load([resultsRoot 'io_ref_ra.mat'],'io_ref_ra')
load([resultsRoot 'io_ref_wt.mat'],'io_ref_wt')
load([resultsRoot 'io_ref_ON.mat'],'io_ref_ON')

outPath = [resultsRoot date filesep 'sweeps_n' num2str(nParamIncrement) filesep];
mkdir(outPath)
% resultsRoot = [resultsRoot 'temp' filesep];
% mkdir(resultsRoot);

% set basic parameters
sweepInfo = struct;
sweepInfo.nParamIncrement = nParamIncrement;
sweepInfo.granularity = granularity;
sweepInfo.n_traces = n_traces;
% sweepInfoRaw.n_traces_per_ap = 25;
sweepInfo.n_keep = 5;
sweepInfo.rate_max = 1; % nothing faster than a second
sweepInfo.keep_prediction_flag = keep_prediction_flag;
sweepInfo.max_ra_time = 7*60;
sweepInfo.calculate_ap_metrics = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%% Load cpHMM results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% for s = 1:length(simTypeCell)
%     sweepInfo = sweepInfo;
    
%     simType = simTypeCell{s};
    
% save simulation type
sweepInfo.simType = simType;

% load markov system parameter info
sweepInfo = getMarkovSystemInfo(sweepInfo);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Call Sweep function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all force 

% record "true" profile
myCluster = parcluster('local');
max_workers = myCluster.NumWorkers;
sweepInfo.NumWorkers = min([24 max_workers]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate ground truth reference curves       
sweepInfo = addGroundTruthFields(sweepInfo, io_ref_ra, io_ref_wt, io_ref_ON);

% initialize vectors to store results
sweepResults = struct;
[sweepInfo, sweepResults] = initializeFitFields(sweepInfo,sweepResults, paramVals);             
sweepResults = initializeSweepValues(sweepInfo, sweepResults, paramVals);              

% call parallel sweep script
sweepInfo.NumWorkers = min([min([28 max_workers]), length(sweepResults)]);
tic
sweepResults= sweep_par_loop_v3(sweepInfo,sweepResults);    
toc

% recombine 
fnames = fieldnames(sweepResults);
if ~keep_prediction_flag
    for f = 1:length(fnames)
        sweepInfo.(fnames{f}) = vertcat(sweepResults.(fnames{f}));
    end
else
    for f = [1:11 length(fnames)]
        sweepInfo.(fnames{f}) = vertcat(sweepResults.(fnames{f}));
    end
    for f = 12:length(fnames)-1
        sweepInfo.(fnames{f}) = cat(3,sweepResults.(fnames{f}));
    end
end
clear sweepResults    
if isempty(paramVals)
    save([outPath 'sweepInfo_' simType '.mat'],'sweepInfo', '-v7.3');
end