% Script to call primary cpHMM wrapper function

clear
close all
warning('off','all') %Shut off Warnings
addpath(genpath('utilities'))

% set basic paths
DataRoot = 'C:\Users\nlamm\Dropbox (Personal)\ProcessedEnrichmentData\';
if ~exist(DataRoot)
  DataRoot = 'S:\Nick\Dropbox\ProcessedEnrichmentData\';
end

project_prefix = '20210928_Oct4_raw_traces_nz';
projectList = dir([DataRoot project_prefix '*']);

master_struct = struct;

project_index = 2;
for p = project_index%1:length(projectList)     

    % load spots struct
    DataPath = [DataRoot filesep projectList(p).name filesep];
    load([DataPath 'spot_struct.mat']);
    
    inferenceInfo = struct;

    % set project identifiers (only applicable if running this on savio)
    inferenceInfo.projectNameCell = {projectList(p).name}; % {'2xDl-Ven_hbP2P-mCh'};

    % set inference options
    inferenceInfo.ProteinBinFlag = 0;
    inferenceInfo.FluoBinFlag = 0;
    %inferenceInfo.timeBins = {[0 60*10],[60*10 60*40]};
    inferenceInfo.timeBins = {[0 45]*60, [135 195]*60}; % should be >= than 15min
    inferenceInfo.apBins = [];%linspace(-.2,.2,10);

    % set core model specs
    inferenceInfo.modelSpecs.nStates = 3; % number of states in system
    inferenceInfo.modelSpecs.nSteps = spot_struct(1).nStepsEst; % number of steps to traverse gene
    inferenceInfo.modelSpecs.alphaFrac =  spot_struct(1).alpha_frac;%1275 / 4670;%

    % other info
    inferenceInfo.AdditionalGroupingVariable = '';%'Stripe'
    inferenceInfo.SampleSize = 10000;
    inferenceInfo.useQCFlag = false;
    inferenceInfo.ignoreNDP = true;
    inferenceInfo.n_localEM = 25;

    % Get basic project info and determing file paths
%     liveProject = LiveEnrichmentProject(inferenceInfo.projectNameCell{1});

    % save
    slashes = regexp(DataPath,'/|\');
    dataDir = DataPath(1:slashes(end-1));
    inferenceDir = [dataDir 'inferenceDirectory' filesep];
    mkdir(inferenceDir)
    save([inferenceDir 'inferenceInfo.mat'],'inferenceInfo')

    % copy bash file to inference directory
    copyfile('run_cpHMM.sh',[inferenceDir 'run_cpHMM.sh'])
end