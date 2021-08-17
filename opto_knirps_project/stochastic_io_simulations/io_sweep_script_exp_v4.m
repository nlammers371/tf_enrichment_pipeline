% Script to scan through parameter space for selected biological variables
clear
close all

addpath(genpath('../utilities'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%% Load experimental data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

projectNameCell = 'optokni_eve4+6_ON'; 
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
sweepInfoRaw.nParamIncrement = 4;
sweepInfoRaw.granularity = 2;
sweepInfoRaw.n_traces = 100;
sweepInfoRaw.n_keep = 5;
sweepInfoRaw.rate_max = 1; % nothing faster than a second
sweepInfoRaw.keep_prediction_flag = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%% Load cpHMM results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% designate simulation type
simTypeCell = {'koff_only_2','kon_only_2','out_only','in_only'};
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
    reactivation_cdf = io_ref_ra.reactivation_time_cdf;
    keep_flags = true(size(reactivation_cdf));%<=0.95;
    reactivation_time = io_ref_ra.reactivation_time_axis;
    
    sweepInfo.reactivation_cdf = reactivation_cdf(keep_flags);
    sweepInfo.reactivation_cdf_full = io_ref_ra.reactivation_time_cdf_full(keep_flags);
    sweepInfo.reactivation_time = reactivation_time(keep_flags);
    sweepInfo.off_frame_ref = io_ref_ra.off_frame_ref;
    
    % mean fluorescence vs. AP
    sweepInfo.mean_fluo_ap = io_ref_wt.fluo_vec_mean;
    
    % observed off times
    sweepInfo.off_time_ap = io_ref_wt.off_time_vec_mean;
    sweepInfo.ap_axis_mean = io_ref_wt.ap_axis_mean;               
    
    % save TF profiles and time vec for RA type
    sweepInfo.tf_profile_array_ra = io_ref_ra.knirps_array;               
    sweepInfo.time_axis_ra = io_ref_ra.time_vec';
    
    % save TF profiles and time vec for WT type
    % calculate a sensible start time
    start_time = ceil(nanmean(io_ref_wt.on_time_vec) / 60)*60;
    [~,start_i] = min(abs(io_ref_wt.time_axis-start_time));
    
    sweepInfo.tf_profile_array_wt = io_ref_wt.knirps_array(start_i:end,:);               
    sweepInfo.wt_start_time = start_time;
    sweepInfo.ap_profile_vec_tf = io_ref_wt.mean_ap;               
    sweepInfo.time_axis_wt = io_ref_wt.time_axis(start_i:end)';        

    % initialize vectors to store results
    [sweepInfo, sweepResults] = initializeSweepValues(sweepInfo);          
        
    % call parallel sweep script
    tic
    sweepResults= sweep_par_loop_v3(sweepInfo,sweepResults,resultsRoot);    
    toc

    % recombine 
    fnames = fieldnames(sweepResults);
    for f = 1:length(fnames)
        sweepInfo.(fnames{f}) = vertcat(sweepResults.(fnames{f}));
    end
    
    clear sweepResults    

    save([resultsRoot 'sweepInfo_' simType '.mat'],'sweepInfo', '-v7.3');
%     
%     clear sweepInfo
%     save([resultsRoot 'gillespie_' simType '.mat'],'gillespie_struct', '-v7.3');
end