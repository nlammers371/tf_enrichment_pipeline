function inferenceOptions = determineInferenceOptions(modelSpecs,varargin)

  inferenceOptions = struct;
  
  %% PATH to HMM FUNCTIONS
  inferenceOptions.modelPath = './utilities';
  
  %% INFERENCE PARAMETERS
  inferenceOptions.savioFlag = 1;
  inferenceOptions.fluo3DFlag = 0;
  inferenceOptions.automaticBinning = 1;
  inferenceOptions.ProteinBinFlag = 1;
  inferenceOptions.intensityBinVar = 'nuclear_protein_vec';
  inferenceOptions.FluoBinFlag = 0;
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
  
  % initialize optional grouping variable
  inferenceOptions.additionalGroupIDs = 1;
  inferenceOptions.AdditionalGroupingVariable = '';
  
  %% Truncated inference (starting mid-trace) or normal inference?
  inferenceOptions.truncInference = 0;    
  inferenceOptions.minDP = 2*inferenceOptions.nSteps; % minimum number of data points for inclusion
  
  
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
  
  if length(inferenceOptions.truncInference) ~= length(inferenceOptions.timeBins)-1
    warning('Updating inference type options to be consistent with time bins')
    inferenceOptions.truncInference = ones(1,length(inferenceOptions.timeBins)-1);
    if inferenceOptions.timeBins(1) == 0
      inferenceOptions.truncInference(1) = 0;
    end
  end
  
  % parpool options
%   if inferenceOptions.savioFlag && ~isfield(inferenceOptions,'maxWorkers')
%     inferenceOptions.maxWorkers = 24;
%   else
%     myCluster = parcluster('local');
%     inferenceOptions.maxWorkers = ceil(myCluster.NumWorkers/2);
%   end
  inferenceOptions.maxWorkers = 24;
  
  % assign binary flags to indicate wheter space or time groupings are used
  inferenceOptions.apBinFlag = length(inferenceOptions.apBins)>2;
  inferenceOptions.timeBinFlag = length(inferenceOptions.timeBins)>2;
  
  % set the number of bootstraps
  if ~isfield(inferenceOptions,'nBoots')
    if inferenceOptions.savioFlag
        inferenceOptions.nBoots = 5; % will run multiple instances on savio
    else  
        inferenceOptions.nBoots = 5;
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