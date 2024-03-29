% Script to call primary cpHMM wrapper function

clear
close all
warning('off','all') %Shut off Warnings
addpath(genpath('utilities'))

% set basic paths
DataRoot = 'C:\Users\nlamm\Dropbox (Personal)\ProcessedEnrichmentData\';
if ~exist(DataRoot)
  DataRoot = 'S:\Nick\Dropbox (Personal)\ProcessedEnrichmentData\';
end

project_prefix = '20220912_KO_experiments';
projectList = dir([DataRoot project_prefix '*']);
master_struct = struct;

%% project_index = 1;
for p = 1:2%length(projectList)     

    % load spots struct
    DataPath = [DataRoot filesep projectList(p).name filesep];
    load([DataPath 'spot_struct.mat']);
    
    inferenceInfo = struct;

    % set project identifiers (only applicable if running this on savio)
    inferenceInfo.projectNameCell = {projectList(p).name}; % {'2xDl-Ven_hbP2P-mCh'};

    % set inference options
    inferenceInfo.ProteinBinFlag = 0;
    inferenceInfo.FluoBinFlag = 0;
    inferenceInfo.singleTraceInference = 0;
    inferenceInfo.alwaysTruncInference = 1;
    
    %inferenceInfo.timeBins = {[0 60*10],[60*10 60*40]};
    inferenceInfo.timeBins = {[0 60]*60}; % should be >= than 15min
    inferenceInfo.apBins = [];%linspace(-.2,.2,10);

    % set core model specs
    inferenceInfo.modelSpecs.nStates = 2; % number of states in system
    if strcmp(project_prefix,'20220701_Oct4_opto')
        inferenceInfo.modelSpecs.nSteps = 3; % number of steps to traverse gene
    else
        inferenceInfo.modelSpecs.nSteps = spot_struct(1).nStepsEst;
    end
    inferenceInfo.modelSpecs.alphaFrac =  spot_struct(1).alpha_frac;%1275 / 4670;%

    % other info
    inferenceInfo.AdditionalGroupingVariable = 'KO_flag';%'Stripe'
    inferenceInfo.SampleSize = 5000;
    inferenceInfo.useQCFlag = false;
    inferenceInfo.ignoreNDP = true;
    inferenceInfo.n_localEM = 25;
    inferenceInfo.nBoots = 25;
    % Get basic project info and determing file paths
%     liveProject = LiveEnrichmentProject(inferenceInfo.projectNameCell{1});

    % save
    slashes = regexp(DataPath,'/|\');
    dataDir = DataPath(1:slashes(end-1));
    inferenceDir = [dataDir 'inferenceDirectory' filesep];
    mkdir(inferenceDir)
    save([inferenceDir 'inferenceInfo.mat'],'inferenceInfo')
    
    main05_conduct_cpHMM_inference({projectList(p).name},'inferenceInfo',inferenceInfo,'customProjectPath',DataPath)
    % copy bash file to inference directory
%     copyfile('run_cpHMM.sh',[inferenceDir 'run_cpHMM.sh'])
end