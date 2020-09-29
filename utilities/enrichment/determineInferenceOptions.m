function inferenceOptions = determineInferenceOptions

  inferenceOptions = struct;
  
  % INFERENCE PARAMETERS
  inferenceOptions.savioFlag = 1;
  inferenceOptions.fluo3DFlag = 0;
  inferenceOptions.automaticBinning = true;
  inferenceOptions.ProteinBinFlag = true;
  inferenceOptions.dpBootstrap = 1;
  inferenceOptions.SampleSize = 5000;
  inferenceOptions.maxWorkers = 24;
  
  %%%%% These options generally remain fixed 
  inferenceOptions.n_localEM = 25; % set num local runs
  inferenceOptions.nStinferenceOptions.epsMax = 500; % set max stinferenceOptions.eps per inference
  inferenceOptions.eps = 1e-4; % set convergence criteria
  inferenceOptions.minDPperInf = 1000; % inference will be aborted if fewer present
  
  % MODEL ARCHITECTURE
  inferenceOptions.nStates = 3; % number of states
  inferenceOptions.nSteps = 7; % number of time stinferenceOptions.eps needed for elongation
  
  % Truncated inference (starting mid-trace) or normal inference?
  inferenceOptions.truncInference = false;
  
  % QC parameters
  inferenceOptions.minDP = 2*inferenceOptions.nSteps;