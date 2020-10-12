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

% get path to inference files
if savioFlag
  inferenceDir = '~/dat/tf_enrichment/inferenceDirectory/';
else
  liveProject = LiveProject(projectNameCell{1});
  slashes = regexp(liveProject.dataPath,'/|\');
  dataDir = liveProject.dataPath(1:slashes(end-1));
  inferenceDir = [dataDir 'inferenceDirectory' filesep];
end

% load dataset with inference info
load([inferenceDir 'inferenceInfo.mat'],'inferenceInfo')

if savioFlag 
  projectNameCell = inferenceInfo.projectNameCell;
end

% set core model specs
modelSpecs = inferenceInfo.modelSpecs;

% generate options cell
options = {'savioFlag',savioFlag,'timeBins',inferenceInfo.timeBins,'SampleSize',inferenceInfo.SampleSize};

% intensity binning
if inferenceInfo.FluoBinFlag
  options(end+1:end+2) = {'intensityBinVar','fluo'};
elseif inferenceInfo.ProteinBinFlag
  options(end+1:end+2) = {'intensityBinVar','nuclear_protein_vec'};  
end

% additional grouping options
if ~isempty(inferenceInfo.AdditionalGroupingVariable)
  options(end+1:end+2) = {'AdditionalGroupingVariable',inferenceInfo.AdditionalGroupingVariable};  
end
rng('shuffle')

for p = randsample(1:length(projectNameCell),length(projectNameCell),false)
  
    % Get basic project info and determing file paths
    [InputDataPath, OutputDataPath] = getDataPaths(savioFlag,projectNameCell{p});
   
    % Call main inference function
    cpHMMInferenceGrouped(InputDataPath,OutputDataPath,modelSpecs,options{:})
end

