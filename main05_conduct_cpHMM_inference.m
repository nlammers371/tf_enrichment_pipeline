% Script to call primary cpHMM wrapper function

clear
close all
warning('off','all') %Shut off Warnings
addpath(genpath('utilities'))

% set project identifier
projectName = '2xDl-Ven_snaBAC-mCh';

% set inference options
ProteinBinFlag = 0;
timeBins = [0 60*15 60*30 60*45];

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

% Get basic project info and determing file paths
[InputDataPath, OutputDataPath] = getDataPaths(savioFlag,projectName);

% Call main inference function
cpHMMInferenceGrouped(projectName,InputDataPath,OutputDataPath,modelSpecs,'savioFlag',savioFlag,'timeBins',timeBins)