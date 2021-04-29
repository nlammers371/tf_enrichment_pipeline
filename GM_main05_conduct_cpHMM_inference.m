% Script to call primary cpHMM wrapper function
function GM_main05_conduct_cpHMM_inference(inferenceDir, varargin)
% MUST BE SUBMITTED FROM THE DIRECTORY CONTAINING run_cpHMM.sh TO WORK
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


% GM 3/1/21: commented this out to be able to submit multiple simultaneous
% jobs
% get path to inference files

% if savioFlag && ~manualInferenceInfo
%   inferenceDir = pwd;
% elseif ~manualInferenceInfo
%   liveProject = LiveEnrichmentProject(projectNameCell{1});
%   slashes = regexp(liveProject.dataPath,'/|\');
%   dataDir = liveProject.dataPath(1:slashes(end-1));
%   inferenceDir = [dataDir 'inferenceDirectory' filesep];
% end

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
if isfield(inferenceInfo, 'dt')
    if ~isempty(inferenceInfo.dt)
        options(end+1:end+2) = {'dt', inferenceInfo.dt};
    end
end
if isfield(inferenceInfo,'apBins')
    options(end+1:end+2) = {'apBins',inferenceInfo.apBins};
end
if isfield(inferenceInfo,'timeBins')
    options(end+1:end+2) = {'timeBins',inferenceInfo.timeBins};
end
if isfield(inferenceInfo,'useQCFlag')
    options(end+1:end+2) = {'useQCFlag',inferenceInfo.useQCFlag};
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

for p = randsample(1:length(projectNameCell),length(projectNameCell),false)

    if ~customProjectFlag
        % Get basic project info and determing file paths
        [InputDataPath, OutputDataPath] = getDataPaths(savioFlag,projectNameCell{p});
    else
        InputDataPath = [customProjectPath projectNameCell{p} filesep];
        OutputDataPath = [customProjectPath projectNameCell{p} filesep];
    end
    % Call main inference function
    cpHMMInferenceGrouped(InputDataPath,OutputDataPath,modelSpecs,options{:})
end

