% Script to call primary cpHMM wrapper function

clear
close all
warning('off','all') %Shut off Warnings
addpath(genpath('utilities'))


inferenceInfo = struct;

% set project identifiers (only applicable if running this on savio)
%inferenceInfo.projectNameCell = {'Rbp1-GFP_eveBAC-mCh'}; % {'2xDl-Ven_hbP2P-mCh'};
inferenceInfo.projectNameCell = {'optokni_eve4+6_WT'};
%inferenceInfo.projectNameCell = {'optokni_eveBAC_ON'};
%inferenceInfo.projectNameCell = {'optokni_eveBAC_ON'};
%inferenceInfo.projectNameCell = {'optokni_eve4+6_MCP-GFP_Homo'};

% set inference options
inferenceInfo.ProteinBinFlag = 0;
inferenceInfo.FluoBinFlag = 0;
%inferenceInfo.timeBins = {[0 60*10],[60*10 60*40]};
%inferenceInfo.timeBins = {[7.5*60 22.5*60],[10*60 25*60],[12.5*60 27.5*60],[15*60 30*60],[17.5*60 32.5*60],[20*60 35*60],[22.5*60 37.5*60]};
%inferenceInfo.timeBins = {[7.5*60 37.5*60]}; % should be longer than 15min

%optokni_eve4+6_MCP-GFP_Homo
%inferenceInfo.apBins = linspace(50,80,16);
%inferenceInfo.apBins = linspace(52.5,70,11);
%inferenceInfo.apBins = ([52.5,70]);
%inferenceInfo.timeBins = {[0 50*60]};

%optokni_eveBAC_WT
inferenceInfo.timeBins = {[6*60 50*60]}; % should be longer than 15min
inferenceInfo.apBins = linspace(-0.12,0.12,11);%linspace(-.2,.2,10);
%inferenceInfo.apBins = ([-0.12 0.12]);

%optokni_eveBAC_ON
%inferenceInfo.timeBins = {[0*60 60*60]}; % should be longer than 15min
%inferenceInfo.apBins = linspace(-0.12,0.12,11);%linspace(-.2,.2,10);

% set core model specs
inferenceInfo.modelSpecs.nStates = 2; % number of states in system
inferenceInfo.modelSpecs.nSteps = 7; % number of steps to traverse gene
inferenceInfo.modelSpecs.alphaFrac = 1302/6444;

% other info
inferenceInfo.AdditionalGroupingVariable = '';%'Stripe'
inferenceInfo.SampleSize = 4500;
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