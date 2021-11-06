function inferenceOptions = determineInferenceOptions(modelSpecs,varargin)

  inferenceOptions = struct;
  
  %% PATH to HMM FUNCTIONS
  inferenceOptions.modelPath = './utilities';
  
  %% INFERENCE PARAMETERS
  inferenceOptions.savioFlag = 1;
  inferenceOptions.fluo3DFlag = 0;
  inferenceOptions.automaticBinning = 1;
  inferenceOptions.ProteinBinFlag = 0;
  inferenceOptions.intensityBinVar = '';
  inferenceOptions.FluoBinFlag = 0;
  inferenceOptions.SampleSize = 5000;
  inferenceOptions.singleTraceInference = 0;
  inferenceOptions.maxWorkers = 20;  
  inferenceOptions.ignoreNDP = 0;
  inferenceOptions.alwaysTruncInference = 0;
  
  %% Core inference options (these generally remain fixed)
  inferenceOptions.n_localEM = 25; % set num local runs
  inferenceOptions.nStepsMax = 500; % set max stinferenceOptions.eps per inference
  inferenceOptions.eps = 1e-4; % set convergence criteria
  inferenceOptions.minDPperInf = 1000; % inference will be aborted if fewer present
  %%%%% 
  
  %% MODEL ARCHITECTURE
  if all([isfield(modelSpecs,'nStates'),isfield(modelSpecs,'nSteps'),isfield(modelSpecs,'alphaFrac')])
    inferenceOptions.nStates = modelSpecs.nStates; % number of states
    inferenceOptions.nSteps = modelSpecs.nSteps; % number of time stinferenceOptions.eps needed for elongation 
    inferenceOptions.alphaFrac = modelSpecs.alphaFrac; 
  else
    error('modelSpecs not fully determined. Be sure it includes "nStates" (number of states), "nSteps" (number time steps to traverse gene), and "alphaFrac" (length of MS2 cassette in time steps)')
  end
  inferenceOptions.alpha = inferenceOptions.alphaFrac*inferenceOptions.nSteps;
  
  %% Binning variables
  inferenceOptions.timeBins = {[0 Inf]};
  inferenceOptions.apBins = [-Inf Inf];
  inferenceOptions.dt = [];
  
  % initialize optional grouping variable
  inferenceOptions.additionalGroupIDs = 1;
  inferenceOptions.AdditionalGroupingVariable = '';
  
  %% Truncated inference (starting mid-trace) or normal inference?  
  inferenceOptions.minDP = 2*inferenceOptions.nSteps; % minimum number of data points for inclusion
  inferenceOptions.useQCFlag = 1;
  
  %% Check for user-specified options  
  for i = 1:2:length(varargin)
    if i ~= length(varargin) && ischar(varargin{i})       
      if isfield(inferenceOptions,varargin{i})
        inferenceOptions.(varargin{i}) = varargin{i+1};  
      else
        warning(['Unrecognized input option "' varargin{i} '". Ignoring.'])
      end
    end    
  end
  
  %% Adjust options as needed
  
  % parpool options
  if inferenceOptions.savioFlag
    inferenceOptions.maxWorkers = 24;
  else
    myCluster = parcluster('local');
    inferenceOptions.maxWorkers = 24;%ceil(myCluster.NumWorkers/2);
  end
  
  % adjust if we're doing single trace inference
  if inferenceOptions.singleTraceInference
      inferenceOptions.minDPperInf = 0;
      inferenceOptions.minDP = 30;%
      inferenceOptions.SampleSize = 1;
  end
  
  % assign binary flags to indicate wheter space or time groupings are used
  if isempty(inferenceOptions.apBins)
      inferenceOptions.apBins = [-Inf Inf];
  end
  inferenceOptions.apBinFlag = any(~ismember(inferenceOptions.apBins,[0 Inf -Inf]));
  inferenceOptions.timeBinFlag = any(~ismember([inferenceOptions.timeBins{:}],[0 Inf]));
  
  if inferenceOptions.alwaysTruncInference
      for t = 1:length(inferenceOptions.timeBins)       
          inferenceOptions.truncInference(t) = 1;            
      end
  else
      for t = 1:length(inferenceOptions.timeBins)
        if inferenceOptions.timeBins{t}(1) == 0
          inferenceOptions.truncInference(t) = 0;    
        else
          inferenceOptions.truncInference(t) = 1;    
        end
      end
  end
  
  % set the number of bootstraps
  if ~isfield(inferenceOptions,'nBoots')
    if inferenceOptions.savioFlag
        inferenceOptions.nBoots = 1; % will run multiple instances on savio
    else  
        inferenceOptions.nBoots = 10;
    end
  end
  
  % update grouping flags
  if strcmpi(inferenceOptions.intensityBinVar,'fluo')
      inferenceOptions.ProteinBinFlag = 0;
      inferenceOptions.FluoBinFlag = 1;
  elseif strcmpi(inferenceOptions.intensityBinVar,'')
      inferenceOptions.ProteinBinFlag = 0;
      inferenceOptions.FluoBinFlag = 0;
  end