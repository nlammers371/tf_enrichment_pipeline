% Script to call primary cpHMM wrapper function

clear
close all
warning('off','all') %Shut off Warnings
addpath(genpath('utilities'))


inferenceInfo = struct;

% set project identifiers (only applicable if running this on savio)
inferenceInfo.projectNameCell = {'hbBAC-MS2-20C'}; % {'2xDl-Ven_hbP2P-mCh'};

% set inference options
inferenceInfo.ProteinBinFlag = 0;
inferenceInfo.FluoBinFlag = 0;
inferenceOptions.fluo3DFlag = 1;
%inferenceInfo.timeBins = {[0 60*10],[60*10 60*40]};
timeBins = cell(1, 8);
for i = 1:length(timeBins)
    timeBins{i} = [((i-1)*5)*60, ((i-1)*5+15)*60];
end
inferenceInfo.timeBins = timeBins; % should be longer than 15min
inferenceInfo.apBins = linspace(22.5, 45, 10);%linspace(-.2,.2,10);

% set core model specs
inferenceInfo.modelSpecs.nStates = 3; %3; % number of states in system
inferenceInfo.modelSpecs.nSteps = 7; % 7; % number of steps to traverse gene
inferenceInfo.modelSpecs.alphaFrac =  1302/6444;%1275 / 4670;%

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
inferenceDir2 = [inferenceDir, inferenceInfo.projectNameCell{1}, filesep];
mkdir(inferenceDir2)
ProjectParamsString = ['w', num2str(inferenceInfo.modelSpecs.nSteps), '_K', ...
    num2str(inferenceInfo.modelSpecs.nStates),'_p', num2str(inferenceInfo.ProteinBinFlag),...
    '_ap', num2str(length(inferenceInfo.apBins)-1),'_t', num2str(length(inferenceInfo.timeBins))];

inferenceDir3 = [inferenceDir2, ProjectParamsString, filesep];
mkdir(inferenceDir3)
save([inferenceDir3 'inferenceInfo.mat'],'inferenceInfo')


% copy bash file to inference directory
copyfile('GM_run_cpHMM.sh',[inferenceDir3 'run_cpHMM.sh'])