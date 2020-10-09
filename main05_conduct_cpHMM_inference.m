% Script to call primary cpHMM wrapper function

clear
close all
warning('off','all') %Shut off Warnings
addpath(genpath('utilities'))

% set project identifier
% projectName = 'EveGtSL';
projectNameCell = {'EveS1Null','EveGtSL','EveGtSL-S1Null','EveWT'};

% set inference options
ProteinBinFlag = 0;
timeBins = [60*25 60*60];

% set core model specs
modelSpecs.nStates = 2; % number of states in system
modelSpecs.nSteps = 4; % number of steps to traverse gene
modelSpecs.alphaFrac = 1302/6444;

% check to see if we are on savio
currentDir = pwd;
savioFlag = 0;
if contains(currentDir,'global/')
  savioFlag = 1;
end
% generate options cell
options = {'savioFlag',savioFlag,'timeBins',timeBins,'intensityBinVar','fluo','AdditionalGroupingVariable','ectopicFlag','SampleSize',2500};
rng('shuffle')

for p = randsample(1:length(projectNameCell),length(projectNameCell),false)
    % Get basic project info and determing file paths
    [InputDataPath, OutputDataPath] = getDataPaths(savioFlag,projectNameCell{p});
   
    % Call main inference function
    cpHMMInferenceGrouped(InputDataPath,OutputDataPath,modelSpecs,options{:})
end

