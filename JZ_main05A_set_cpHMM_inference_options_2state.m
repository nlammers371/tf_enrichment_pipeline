% Script to call primary cpHMM wrapper function

clear
close all
warning('off','all') %Shut off Warnings
addpath(genpath('utilities'))


inferenceInfo = struct;

% set project identifiers (only applicable if running this on savio)
%inferenceInfo.projectNameCell = {'optokni_eve4+6_MCP-GFP_Het'};
%inferenceInfo.projectNameCell = {'optokni_eve4+6_MCP-GFP_Homo'};
%inferenceInfo.projectNameCell = {'optokni_eve4+6_WT'};
%inferenceInfo.projectNameCell = {'optokni_eve4+6_WT_FULL'};
inferenceInfo.projectNameCell = {'optokni_eve4+6_ON_CONST'};
%inferenceInfo.projectNameCell = {'optokni_eve4+6_ON_LOW_FULL'};
%inferenceInfo.projectNameCell = {'optokni_eve4+6_ON_LOW'};

% set inference options
inferenceInfo.ProteinBinFlag = 1;
inferenceInfo.FluoBinFlag = 0;

%Inference for HIGH and LOW comparison (v1)
inferenceInfo.timeBins = {[15*60 30*60]}; % should be longer than 15min
inferenceInfo.apBins = [-0.03 0.03];

%Inference for HIGH and LOW comparison (v2)
%inferenceInfo.timeBins = {[15*60 30*60]}; % should be longer than 15min
%inferenceInfo.apBins = [-0.02 0.02];

%Inference for WT, HIGH and LOW comparison (v2)
%inferenceInfo.timeBins = {[10*60 30*60]}; % should be longer than 15min
%inferenceInfo.apBins = [-0.025 0.025];%linspace(-.2,.2,10);

%Long inference, paper figure, optokni_eve4+6_WT and ON_CONST
%inferenceInfo.timeBins = {[6*60 50*60]}; % should be longer than 15min
%inferenceInfo.apBins = linspace(-0.12,0.12,10);%linspace(-.2,.2,10);

%Shorter range, test: optokni_eve4+6_WT and ON_CONST
%inferenceInfo.timeBins = {[10*60 25*60]}; % should be longer than 15min
%inferenceInfo.apBins = linspace(-0.12,0.12,10);%linspace(-.2,.2,10);

%Temporal inference
%inferenceInfo.timeBins = {[10*60 25*60],[12.5*60 27.5*60],[15*60 30*60],[17.5*60 32.5*60],[20*60 35*60],[22.5*60 37.5*60],[25*60 40*60]}; % should be longer than 15min
% 11 bin
%inferenceInfo.timeBins = {[0*60 15*60],[2.5*60 17.5*60],[5*60 20*60],[7.5*60 22.5*60],[10*60 25*60],[12.5*60 27.5*60],[15*60 30*60],[17.5*60 32.5*60],[20*60 35*60],[22.5*60 37.5*60],[25*60 40*60]};
%inferenceInfo.apBins = linspace(-0.12,0.12,10);%linspace(-.2,.2,10);
%inferenceInfo.apBins = ([-0.025 0.025]);

%% old parameters
%optokni_eve4+6_ON_CONST
%inferenceInfo.timeBins = {[6*60 50*60]}; % should be longer than 15min
%inferenceInfo.apBins = linspace(-0.12,0.12,8);%linspace(-.2,.2,10);

%optokni_eve4+6_MCP-GFP_Homo
%inferenceInfo.apBins = linspace(55,67.5,7);
%inferenceInfo.timeBins = {[0 50*60]};

%optokni_eve4+6_MCP-GFP_Het
%inferenceInfo.apBins = linspace(55,67.5,5);
%inferenceInfo.timeBins = {[0 50*60]};

%optokni_eve4+6_WT
%inferenceInfo.timeBins = {[6*60 50*60]}; % should be longer than 15min
%inferenceInfo.apBins = linspace(-0.12,0.12,8);%linspace(-.2,.2,10);
%inferenceInfo.apBins = ([-0.12 0.12]);

% set core model specs
inferenceInfo.modelSpecs.nStates = 2; % number of states in system
inferenceInfo.modelSpecs.nSteps = 7; % number of steps to traverse gene
inferenceInfo.modelSpecs.alphaFrac = 1302/6444;

% other info
inferenceInfo.AdditionalGroupingVariable = '';%'Stripe'
inferenceInfo.SampleSize = 2500;
inferenceInfo.useQCFlag = true;

inferenceInfo.n_localEM = 25;

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