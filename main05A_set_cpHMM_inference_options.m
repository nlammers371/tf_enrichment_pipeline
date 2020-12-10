% Script to call primary cpHMM wrapper function

clear
close all
warning('off','all') %Shut off Warnings
addpath(genpath('utilities'))


inferenceInfo = struct;

% set project identifiers (only applicable if running this on savio)
inferenceInfo.projectNameCell = {'2xDl-Ven_hbP2P-mCh'};

% set inference options
inferenceInfo.ProteinBinFlag = 1;
inferenceInfo.FluoBinFlag = 0;
%inferenceInfo.timeBins = {[0 60*10],[60*10 60*40]};
inferenceInfo.timeBins = {[0 Inf]}; % should be longer than 15min
inferenceInfo.apBins = [0 Inf];

% set core model specs
inferenceInfo.modelSpecs.nStates = 3; % number of states in system
inferenceInfo.modelSpecs.nSteps = 7; % number of steps to traverse gene
inferenceInfo.modelSpecs.alphaFrac = 1275 / 4670;%1302/6444;%

% other info
inferenceInfo.AdditionalGroupingVariable = '';%'Stripe'
inferenceInfo.SampleSize = 3000;
inferenceInfo.useQCFlag = true;

% Get basic project info and determing file paths
liveProject = LiveEnrichmentProject(inferenceInfo.projectNameCell{1});

% save
slashes = regexp(liveProject.dataPath,'/|\');
dataDir = liveProject.dataPath(1:slashes(end-1));
inferenceDir = [dataDir 'inferenceDirectory' filesep];
mkdir(inferenceDir)
save([inferenceDir 'inferenceInfo.mat'],'inferenceInfo')

% copy bash file to inference directory
copyfile('run_cpHMM.sh',[inferenceDir 'run_cpHMM.sh'])