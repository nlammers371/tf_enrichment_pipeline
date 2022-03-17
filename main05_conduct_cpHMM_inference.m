% Script to call primary cpHMM wrapper function
function main05_conduct_cpHMM_inference(varargin)

close all
warning('off','all') %Shut off Warnings
addpath(genpath('utilities'))


% check to see if we are on savio
currentDir = pwd;
savioFlag = 0;
if contains(currentDir,'global/')
  savioFlag = 1;
end

if ~savioFlag && isempty(varargin)
  warning('Cell array of project names must be passed when not running function on Savio')
elseif ~isempty(varargin)
  projectNameCell = varargin{1};
end

% process other options
for i = (1+~savioFlag):2:(numel(varargin)-1)  
    if i ~= numel(varargin)        
        eval([varargin{i} '=varargin{i+1};']);                
    end    
end

manualInferenceInfo = exist('inferenceInfo','var');
customProjectFlag = exist('customProjectPath','var');

% get path to inference files
if savioFlag && ~manualInferenceInfo
  inferenceDir = '~/dat/tf_enrichment/inferenceDirectory/';
elseif ~manualInferenceInfo
  try
    liveProject = LiveEnrichmentProject(projectNameCell{1});
    dataPath = liveProject.dataPath;
  catch
    customProjectFlag = true;
    if savioFlag
      dataPath = ['~/dat/tf_enrichment/' projectNameCell{1} filesep];
    else % NL: this is just for testing purposes
      dataPath = ['S:/Nick/Dropbox (Personal)/ProcessedEnrichmentData/' projectNameCell{1} filesep];
    end
    customProjectPath = dataPath;
  end
  slashes = regexp(dataPath,'/|\');
  dataDir = dataPath(1:slashes(end-1));
  inferenceDir = [dataDir 'inferenceDirectory' filesep];
end

if ~manualInferenceInfo
    % load dataset with inference info
    load([inferenceDir 'inferenceInfo.mat'],'inferenceInfo')
end

if savioFlag 
  projectNameCell = inferenceInfo.projectNameCell;
end

% set core model specs
modelSpecs = inferenceInfo.modelSpecs;

% generate options cell
options = {'savioFlag',savioFlag,'SampleSize',inferenceInfo.SampleSize};
if isfield(inferenceInfo,'apBins')
    options(end+1:end+2) = {'apBins',inferenceInfo.apBins};
end
if isfield(inferenceInfo,'n_localEM')
    options(end+1:end+2) = {'n_localEM',inferenceInfo.n_localEM};
end
if isfield(inferenceInfo,'timeBins')
    options(end+1:end+2) = {'timeBins',inferenceInfo.timeBins};
end
if isfield(inferenceInfo,'useQCFlag')
    options(end+1:end+2) = {'useQCFlag',inferenceInfo.useQCFlag};
end
if isfield(inferenceInfo,'ignoreNDP')
    options(end+1:end+2) = {'ignoreNDP',inferenceInfo.ignoreNDP};
end
if isfield(inferenceInfo,'singleTraceInference')
    options(end+1:end+2) = {'singleTraceInference',inferenceInfo.singleTraceInference};
end
if isfield(inferenceInfo,'alwaysTruncInference')
    options(end+1:end+2) = {'alwaysTruncInference',inferenceInfo.alwaysTruncInference};
end

% intensity binning
if inferenceInfo.FluoBinFlag
  options(end+1:end+2) = {'intensityBinVar','fluo'};
  options(end+1:end+2) = {'FluoBinFlag',1};
elseif inferenceInfo.ProteinBinFlag
  options(end+1:end+2) = {'intensityBinVar','rawNCProteinInterp'};  
  options(end+1:end+2) = {'ProteinBinFlag',1};
end

% additional grouping options
if ~isempty(inferenceInfo.AdditionalGroupingVariable)
  options(end+1:end+2) = {'AdditionalGroupingVariable',inferenceInfo.AdditionalGroupingVariable};  
end

% randomly seeds random number generator
rng('shuffle')

for p = 1:length(projectNameCell)%randsample(1:length(projectNameCell),length(projectNameCell),false)
    
    if ~customProjectFlag
        % Get basic project info and determing file paths
        [InputDataPath, OutputDataPath] = getDataPaths(savioFlag,projectNameCell{p});
    else
        InputDataPath = [customProjectPath  filesep];
        OutputDataPath = [customProjectPath filesep];
    end
    % Call main inference function
    cpHMMInferenceGrouped(InputDataPath,OutputDataPath,modelSpecs,options{:})
end

