% Script to scan through parameter space for selected biological variables
clear
close all

addpath(genpath('../utilities'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%% Load experimental data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% designate simulation type
simType = 'koff_only_2';%'out_only','in_only','kon_only_2';
projectName = 'optokni_eve4+6_ON';

try
  liveProject = LiveEnrichmentProject(projectName);
  resultsRoot = [liveProject.dataPath filesep];
catch  
  resultsRoot = 'C:\Users\nlamm\Dropbox (Personal)\ProcessedEnrichmentData\optokni_eve4+6_ON\';
end
% load data
load([resultsRoot 'io_ref_struct.mat'])
load([resultsRoot 'sweepInfo_' simType '.mat'])

% find best param vals to use for initialization
[~,mi_pon] = min(sweepInfo.objective_val_p_on);
bestParamVals = sweepInfo.param_fit_array(mi_pon,:);

mcmcInfo = sweepInfo;

% set basic parameters
mcmcInfo.nParamIncrement = 15;
mcmcInfo.granularity = 5;
mcmcInfo.n_traces = 100;
mcmcInfo.n_keep = 10;

mcmcInfo.nIterations = 1e2;
mcmcInfo.prop_sigma = 0.05;
mcmcInfo.nBoots = 100;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Initialize Fields %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all force 

% record "true" profile
myCluster = parcluster('local');
max_workers = myCluster.NumWorkers;

mcmcInfo.NumWorkers = min([24 max_workers]);
mcmcInfo.nFit = sum(mcmcInfo.fitFlags);


% set mcmc sampling hyperparameters
mcmcInfo.param_bounds = sweepInfo.param_bounds([1 end],:);%linspace(0.5, 20, mcmcInfo.nParamIncrement);
mcmcInfo.param_bounds(1,:) = mcmcInfo.param_bounds(1,:)/10;
mcmcInfo.param_bounds(2,:) = mcmcInfo.param_bounds(2,:)*10;
mcmcInfo.fitFlags = mcmcInfo.fitFlags==1;

% initialize vectors to store results
mcmcInfo.param_fit_array = NaN(mcmcInfo.nIterations,length(mcmcInfo.paramList));
mcmcInfo.param_fit_array(:,~mcmcInfo.fitFlags) = repmat(bestParamVals(~mcmcInfo.fitFlags),size(mcmcInfo.param_fit_array,1),1);
mcmcInfo.param_fit_array(1,mcmcInfo.fitFlags) = bestParamVals(:,mcmcInfo.fitFlags);

% track profile
mcmcInfo.p_on_fit_array = NaN(mcmcInfo.nIterations,length(mcmcInfo.p_on_true));
mcmcInfo.fluo_fit_array = NaN(mcmcInfo.nIterations,length(mcmcInfo.p_on_true));

% track function fit
mcmcInfo.objective_val_p_on = NaN(mcmcInfo.nIterations,1);
mcmcInfo.move_flags = zeros(mcmcInfo.nIterations,1);
mcmcInfo.objective_val_fluo = NaN(mcmcInfo.nIterations,1);

% generate reference array for indexing
mcmcInfo.ref_array = repmat(reshape(1:size(mcmcInfo.tf_profile_array,1),1,1,[]),mcmcInfo.nBoots,mcmcInfo.n_traces);
mcmcInfo.ref_array = mcmcInfo.ref_array(:,:,mcmcInfo.t_filter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%% Perform MCMC fits %%%%%%%%%%%%%%%%%%%%%%%%%%%%
temperature = 25;

ref_array = mcmcInfo.ref_array;
fit_flags = mcmcInfo.fitFlags;
param_bounds = log10(mcmcInfo.param_bounds);
prop_sigma = mcmcInfo.prop_sigma;

% generate initial prediction and logL score
mcmcInfo.step = 1;
simInfoPD = io_prediction_wrapper_v2(mcmcInfo);    
% estimate error bars with bootstrap resampling
logL_vec = calculate_bootstrap_logL(mcmcInfo,simInfoPD);
mcmcInfo.objective_val_p_on(1) = sum(logL_vec);

for mcmc_step = 2:mcmcInfo.nIterations
    % increment
    mcmcInfo.step = mcmc_step;
    
    % get current guess
    current_vals = log10(mcmcInfo.param_fit_array(mcmc_step-1,fit_flags));
    
    % calculate bounds
    lb_array = (param_bounds(1,fit_flags)-current_vals(:,fit_flags))/prop_sigma;
    ub_array = (param_bounds(2,fit_flags)-current_vals(:,fit_flags))/prop_sigma;
    
    % generate variants
    new_params_prop = 10.^(current_vals+reshape(prop_sigma*trandn(lb_array,ub_array),[],sum(fit_flags)));
    
    % generate new profile prediction
    mcmcInfo.param_fit_array(mcmc_step,fit_flags) = new_params_prop;
    simInfoPD = io_prediction_wrapper_v2(mcmcInfo);
    
    % estimate error bars with bootstrap resampling
    [logL_vec,simInfoPD] = calculate_bootstrap_logL(mcmcInfo,simInfoPD);
    new_logL = sum(logL_vec);
    % perform MH step
    move_flag = exp((sum(logL_vec) - mcmcInfo.objective_val_p_on(mcmc_step-1))/temperature)> rand();
    
    if move_flag
        mcmcInfo.objective_val_p_on(mcmc_step) = new_logL;
        % document move
        mcmcInfo.move_flags(mcmc_step) = 1;
    else
        mcmcInfo.objective_val_p_on(mcmc_step) = mcmcInfo.objective_val_p_on(mcmc_step-1);
        % reset param values
        mcmcInfo.param_fit_array(mcmc_step,fit_flags) = mcmcInfo.param_fit_array(mcmc_step-1,fit_flags);        
        % document move
        mcmcInfo.move_flags(mcmc_step) = 0;
    end
end
%%
% recombine
mcmcInfo.param_val_array = vertcat(sweepTemp.param_fit_array);
mcmcInfo.fluo_fit_array = [sweepTemp.fluo_fit_array];
mcmcInfo.p_on_fit_array = [sweepTemp.p_on_fit_array];
mcmcInfo.objective_val_p_on = vertcat(sweepTemp.objective_val_p_on);
mcmcInfo.objective_val_fluo = vertcat(sweepTemp.objective_val_fluo);
gillespie_struct = [sweepTemp.gillespie];

