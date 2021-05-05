% Script to call primary cpHMM wrapper function

clear
close all
warning('off','all') %Shut off Warnings
addpath(genpath('utilities'))


inferenceInfo = struct;

% set project identifiers (only applicable if running this on savio)
%inferenceInfo.projectNameCell = {'optokni_eve4+6_MCP-GFP_Het'};
%inferenceInfo.projectNameCell = {'optokni_eve4+6_MCP-GFP_Homo'};
inferenceInfo.projectNameCell = {'optokni_eve4+6_WT'};

% set inference options
inferenceInfo.ProteinBinFlag = 0;
inferenceInfo.FluoBinFlag = 0;

%optokni_eve4+6_MCP-GFP_Homo
%inferenceInfo.apBins = linspace(55,67.5,7);
%inferenceInfo.timeBins = {[0 50*60]};

%optokni_eve4+6_MCP-GFP_Het
inferenceInfo.apBins = linspace(55,67.5,5);
inferenceInfo.timeBins = {[0 50*60]};

%optokni_eveBAC_WT
%inferenceInfo.timeBins = {[6*60 50*60]}; % should be longer than 15min
%inferenceInfo.apBins = linspace(-0.12,0.12,8);%linspace(-.2,.2,10);
%inferenceInfo.apBins = ([-0.12 0.12]);

% set core model specs
inferenceInfo.modelSpecs.nStates = 2; % number of states in system
inferenceInfo.modelSpecs.nSteps = 7; % number of steps to traverse gene
inferenceInfo.modelSpecs.alphaFrac = 1302/6444;

% other info
inferenceInfo.AdditionalGroupingVariable = '';%'Stripe'
inferenceInfo.SampleSize = 4500;
inferenceInfo.useQCFlag = true;

inferenceInfo.n_localEM = 75;

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