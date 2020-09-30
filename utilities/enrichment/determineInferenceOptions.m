function inferenceOptions = determineInferenceOptions(varargin)

  inferenceOptions = struct;
  
  %% PATH to HMM FUNCTIONS
  inferenceOptions.modelPath = './utilities';
  
  %% INFERENCE PARAMETERS
  inferenceOptions.savioFlag = 1;
  inferenceOptions.fluo3DFlag = 0;
  inferenceOptions.automaticBinning = 1;
  inferenceOptions.ProteinBinFlag = 1;
  inferenceOptions.dpBootstrap = 1;
  inferenceOptions.SampleSize = 5000;
  inferenceOptions.maxWorkers = 24;  
  
  %% Core inference options (these generally remain fixed)
  inferenceOptions.n_localEM = 25; % set num local runs
  inferenceOptions.nStinferenceOptions.epsMax = 500; % set max stinferenceOptions.eps per inference
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
  inferenceOptions.timeBins = [0 Inf];
  inferenceOptions.apBins = [0 Inf];
  
  %% Truncated inference (starting mid-trace) or normal inference?
  inferenceOptions.truncInference = 0;    
  inferenceOptions.minDP = 2*inferenceOptions.nSteps; % minimum number of data points for inclusion
  
  
  %% Check for user-specified options
  for i = 1:length(varargin)-1
    if i ~= length(varargin) && ischar(varargin{i})       
      if isfield(infrenceOptions,varargin{i})
        inferenceOptions.(varargin{i}) = varargin{i+1};  
      else
        warining(['Unrecognized input option "' varargin{i} '". Ignoring.'])
      end
    end    
  end
  
  %% Adjust options as needed
  if ~inferenceOptions.dpBootstrap
    warning('Bootstrap option not selected. Setting inferenceOptions.nBoots to 1')
    inferenceOptions.nBoots = 1;
  end  
  
  if length(inferenceOptions.truncInference) ~= length(inferenceOptions.timeBins)-1
    warning('Updating inference type options to be consistent with time bins')
    inferenceOptions.truncInference = ones(1,length(inferenceOptions.timeBins)-1);
    if inferenceOptions.timeBins(1) == 0
      inferenceOptions.truncInference(1) = 0;
    end
  end
  
  % parpool options
  if inferenceOptions.savioFlag && ~isfield(inferenceOptions,'maxWorkers')
    inferenceOptions.maxWorkers = 24;
  else
    myCluster = parcluster('local');
    inferenceOptions.maxWorkers = ceil(myCluster.NumWorkers/2);
  end
  
  % assign binary flags to indicate wheter space or time groupings are used
  inferenceOptions.apBinFlag = length(inferenceOptions.apBins)>2;
  inferenceOptions.timeBinFlag = length(inferenceOptions.timeBins)>2;
  
  % set the number of bootstraps
  if ~isfield(inferenceOptions,'nBoots')
    if inferenceOptions.ProteinBinFlag && inferenceOptions.savioFlag
        inferenceOptions.nBoots = 1; % will run multiple instances on savio
    else  
        inferenceOptions.nBoots = 5;
    end
  end
  