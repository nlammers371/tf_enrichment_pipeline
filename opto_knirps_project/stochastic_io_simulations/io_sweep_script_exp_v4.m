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
sweepInfoRaw.nParamIncrement = 20;
sweepInfoRaw.granularity = 2;
sweepInfoRaw.n_traces = 100;
sweepInfoRaw.n_keep = 5;
sweepInfoRaw.rate_max = 1; % nothing faster than a second
% sweepInfoRaw.time_full = ((1:length(io_ref_struct.time_vec))*systemParams.deltaT)/60;% size(io_ref_struct.on_off_array,1);
% sweepInfoRaw.seq_length = length(systemParams.time_full);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%% Load cpHMM results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% designate simulation type
simTypeCell = {'koff_only_2','kon_only_2','out_only','in_only'};
tfDependentParamCell = {'koff','kon','ks','ka'};

for s = 1%length(simTypeCell)
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
    sweepInfo.tf_profile_array_wt = io_ref_wt.knirps_array;               
    sweepInfo.ap_profile_vec_tf = io_ref_wt.mean_ap;               
    sweepInfo.time_axis_wt = io_ref_wt.time_axis';        

    % initialize vectors to store results
    sweepInfo = initializeSweepValues(sweepInfo);   

    % track fits to observables of interest
    sweepInfo.ra_fit = NaN(sweepInfo.nIterations,1);
    sweepInfo.mean_fluo_fit = NaN(sweepInfo.nIterations,1);
    sweepInfo.off_time_fit = NaN(sweepInfo.nIterations,1);
    
    % call parallel sweep script
    tic
    sweepTemp = sweep_par_loop_v3(sweepInfo,resultsRoot);    
    toc

    % recombine    
    sweepInfo.fluo_fit_array = [sweepTemp.fluo_fit_array];
    sweepInfo.p_on_fit_array = [sweepTemp.p_on_fit_array];
    sweepInfo.fluo_raw_fit_array = [sweepTemp.fluo_raw_fit_array];
    
    sweepInfo.off_fluo_array = [sweepTemp.off_fluo_array];
    sweepInfo.reactivation_time_vec = [sweepTemp.reactivation_time_vec];
    sweepInfo.fluo_raw_fit_array = [sweepTemp.fluo_raw_fit_array];
    sweepInfo.tf_rate_trends = [sweepTemp.tf_dependent_rate];
    sweepInfo.fluo_obs_only_array =  [sweepTemp.fluo_obs_only_array];
    sweepInfo.reactivation_time_cdf_array = vertcat(sweepTemp.reactivation_time_cdf)';
    sweepInfo.reactivation_time_axis = sweepTemp(1).reactivation_time_cdf';
    
    clear sweepTemp
    
%     sweepInfo.objective_val_p_on = vertcat(sweepTemp.objective_val_p_on);
%     sweepInfo.objective_val_fluo = vertcat(sweepTemp.objective_val_fluo);
%     sweepInfo.objective_val_fluo_raw = vertcat(sweepTemp.objective_val_fluo_raw);
% %     gillespie_struct = [sweepTemp.gillespie];

    save([resultsRoot 'sweepInfo_' simType '.mat'],'sweepInfo', '-v7.3');
%     
%     clear sweepInfo
%     save([resultsRoot 'gillespie_' simType '.mat'],'gillespie_struct', '-v7.3');
end