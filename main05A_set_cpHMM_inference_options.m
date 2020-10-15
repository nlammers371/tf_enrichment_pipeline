% Script to call primary cpHMM wrapper function

clear
close all
warning('off','all') %Shut off Warnings
addpath(genpath('utilities'))


inferenceInfo = struct;
% set project identifiers (only applicable if running this on savio)
inferenceInfo.projectNameCell = {'EveS1Null','EveGtSL','EveGtSL-S1Null','EveWt'};

% set inference options
inferenceInfo.ProteinBinFlag = 0;
inferenceInfo.FluoBinFlag = 1;
inferenceInfo.timeBins = [60*25 60*60];

% set core model specs
inferenceInfo.modelSpecs.nStates = 3; % number of states in system
inferenceInfo.modelSpecs.nSteps = 7; % number of steps to traverse gene
inferenceInfo.modelSpecs.alphaFrac = 1302/6444;

% other info
inferenceInfo.AdditionalGroupingVariable = 'ectopicFlag';%'Stripe'
inferenceInfo.SampleSize = 2500;

% Get basic project info and determing file paths
liveProject = LiveProject(inferenceInfo.projectNameCell{1});

% save
slashes = regexp(liveProject.dataPath,'/|\');
dataDir = liveProject.dataPath(1:slashes(end-1));
inferenceDir = [dataDir 'inferenceDirectory' filesep];
mkdir(inferenceDir)
save([inferenceDir 'inferenceInfo.mat'],'inferenceInfo')

% copy bash file to inference directory
copyfile('run_cpHMM.sh',[inferenceDir 'run_cpHMM.sh'])