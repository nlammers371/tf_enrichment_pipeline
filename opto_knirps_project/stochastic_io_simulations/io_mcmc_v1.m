% Script to scan through parameter space for selected biological variables
clear
close all

addpath(genpath('../utilities'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%% Load experimental data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

resultsRoot = 'S:\Nick\Dropbox\ProcessedEnrichmentData\parameterSweeps\';
if ~isfolder(resultsRoot)
  resultsRoot = 'C:\Users\nlamm\Dropbox (Personal)\ProcessedEnrichmentData\parameterSweeps\';
end

% load data
load([resultsRoot 'io_ref_ra.mat'])
load([resultsRoot 'io_ref_wt.mat'])

% resultsRoot = [resultsRoot 'temp' filesep];
% mkdir(resultsRoot);

% set basic parameters
sweepInfoRaw = struct;
sweepInfoRaw.n_chains = 1;
sweepInfoRaw.nParamIncrement = 250; % sets granularity of space to explore
sweepInfoRaw.granularity_ra = 1;
sweepInfoRaw.granularity_wt = 1;
sweepInfoRaw.n_traces_ra = 100;
sweepInfoRaw.n_traces_per_ap = 25;
sweepInfoRaw.n_keep = 5;
sweepInfoRaw.rate_max = 1; % nothing faster than a second
sweepInfoRaw.keep_prediction_flag = false;
sweepInfoRaw.max_ra_time = 7*60;
sweepInfoRaw.calculate_ap_metrics = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%% Load cpHMM results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% designate simulation type
simTypeCell = {'out_only','in_only','koff_only_2','kon_only_2'};
tfDependentParamCell = {'koff','kon','ks','ka'};

for s = 1:length(simTypeCell)
    sweepInfo = sweepInfoRaw;
    
    simType = simTypeCell{s};
    
    % save simulation type
    sweepInfo.simType = simType;
    
    % load markov system parameter info
    sweepInfo = getMarkovSystemInfo(sweepInfo);
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%% Call Sweep function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    close all force 

    % record "true" profile
    myCluster = parcluster('local');
    max_workers = myCluster.NumWorkers;
    sweepInfo.NumWorkers = min([24 max_workers]);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate ground truth reference curves       
    sweepInfo = addGroundTruthFields(sweepInfo, io_ref_ra, io_ref_wt);

    % initialize vectors to store results
%     sweepResults = struct;
    sweepInfo = initializeFitFieldsMCMC(sweepInfo);
%     sweepResults = initializeSweepValues(sweepInfo, sweepResults);              
    
    % call parallel sweep script
    tic
    sweepResults = run_mcmc_sampling(sweepInfo,sweepResults);    
    toc

    % recombine 
    fnames = fieldnames(sweepResults);
    for f = 1:length(fnames)
        sweepInfo.(fnames{f}) = vertcat(sweepResults.(fnames{f}));
    end
    
    clear sweepResults    

    save([resultsRoot 'sweepInfo_' simType '.mat'],'sweepInfo', '-v7.3');
    
    % identify best thousand or best 1% (whichever is less) and run
    % sweeps that save key info
    sweepInfoBest = sweepInfoRaw;
    
    % save simulation type
    sweepInfoBest.simType = simType;    
    
    % load markov system parameter info
    sweepInfoBest = getMarkovSystemInfo(sweepInfoBest);  
    
    % generate ground truth reference curves       
    sweepInfoBest = addGroundTruthFields(sweepInfoBest, io_ref_ra, io_ref_wt);
    
    % generate array of weights to apply to different metrics. This will
    % allow us to extract different versions of "optimal" networks
    sweepInfoBest.fit_fields_to_use = {'ra_full_fit_R2','mean_fluo_fit_R2','off_time_fit_R2'};
    n_fit_fields = length(sweepInfoBest.fit_fields_to_use);
    sweepInfoBest.n_raw = min([size(sweepInfo.param_val_vec,1) max([100 min([round(0.01*size(sweepInfo.param_val_vec,1)), 1e3])])]);
    sweepInfoBest.resultsRoot = resultsRoot;
    
    % now iterate throug and find the best-scoring parameter set for each set
    % of weights
    sweepInfoBest = sweepBestPerformers(sweepInfoBest,sweepInfo);   

    save([resultsRoot 'sweepInfoBest_' simType '.mat'],'sweepInfoBest', '-v7.3');
end