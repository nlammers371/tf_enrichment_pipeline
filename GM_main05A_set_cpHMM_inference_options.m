% Script to call primary cpHMM wrapper function

clear
close all
warning('off','all') %Shut off Warnings
addpath(genpath('utilities'))


inferenceInfo = struct;

ProjectName = 'hbBAC-MS2-17_5C-Approved';

inferenceInfo.projectNameCell = {ProjectName};

% set inference options
inferenceInfo.ProteinBinFlag = 0;
inferenceInfo.FluoBinFlag = 0;
inferenceInfo.fluo3DFlag = 0;
%inferenceInfo.timeBins = {[0 60*10],[60*10 60*40]};
timeBins = cell(1, 6);
BinSepMin = 12.4;
BinWidthMin =37.1;
for i = 1:length(timeBins)
    timeBins{i} = [((i-1)*BinSepMin)*60, ((i-1)*BinSepMin+BinWidthMin)*60];
    disp(['i = ' num2str(i), ', Bin Max: ', num2str(timeBins{i}(2)/60)])
end
inferenceInfo.timeBins = timeBins; % should be longer than 15min\
%inferenceInfo.timeBins = {[0 Inf]};
inferenceInfo.apBins = linspace(10, 60, 21);%linspace(-.2,.2,10);
inferenceInfo.dt = 45; % leave empty [] to use default 

% set core model specs
inferenceInfo.modelSpecs.nStates = 3; %3; % number of states in system
inferenceInfo.modelSpecs.nSteps = 7; % 7; % number of steps to traverse gen
inferenceInfo.modelSpecs.alphaFrac =  1302/6444;%1275 / 4670;%

% other info
inferenceInfo.AdditionalGroupingVariable = '';%'Stripe'
inferenceInfo.SampleSize = 4500;
inferenceInfo.useQCFlag = true;

% Get basic project info and determing file paths
if contains(ProjectName(end), '/') | contains(ProjectName(end), '\') 
    ProjectName = ProjectName(1:end-1);
end

if ~contains(ProjectName, '/') & ~contains(ProjectName, '\') 
    liveProject = LiveEnrichmentProject(inferenceInfo.projectNameCell{1});
else
     liveProject = LiveEnrichmentProject(fileparts(inferenceInfo.projectNameCell{1}));
end

% save
slashes = regexp(liveProject.dataPath,'/|\');
dataDir = liveProject.dataPath(1:slashes(end-1));
inferenceDir = [dataDir 'inferenceDirectory' filesep];
mkdir(inferenceDir)
inferenceDir2 = [inferenceDir, ProjectName, filesep];
mkdir(inferenceDir2)
ProjectParamsString = ['w', num2str(inferenceInfo.modelSpecs.nSteps), '_K', ...
    num2str(inferenceInfo.modelSpecs.nStates),'_p', num2str(inferenceInfo.ProteinBinFlag),...
    '_ap', num2str(length(inferenceInfo.apBins)-1),'_t', num2str(length(inferenceInfo.timeBins))];
if ~isempty(inferenceInfo.dt)
     ProjectParamsString = [ProjectParamsString,'_dt', num2str(inferenceInfo.dt)];
end
if inferenceInfo.fluo3DFlag == 1
    ProjectParamsString = [ProjectParamsString,'_f3d'];
end



inferenceDir3 = [inferenceDir2, ProjectParamsString, filesep];
mkdir(inferenceDir3)
save([inferenceDir3 'inferenceInfo.mat'],'inferenceInfo')


% copy bash file to inference directory
copyfile('GM_run_cpHMM.sh',[inferenceDir3 'run_cpHMM.sh'])