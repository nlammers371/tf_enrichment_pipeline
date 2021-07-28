% Script to scan through parameter space for selected biological variables
clear
close all

addpath(genpath('../utilities'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%% Load experimental data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

projectName = 'optokni_eve4+6_ON'; 
try
  liveProject = LiveEnrichmentProject(projectName);
  resultsRoot = [liveProject.dataPath filesep];
catch  
  resultsRoot = 'C:\Users\nlamm\Dropbox (Personal)\ProcessedEnrichmentData\optokni_eve4+6_ON\';
end
% load data
load([resultsRoot 'io_ref_struct.mat'])


% set basic parameters
sweepInfo = struct;
sweepInfo.nParamIncrement = 25;
sweepInfo.granularity = 2;
sweepInfo.n_traces = 100;
sweepInfo.n_keep = 5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify key system characteristics
systemParams = struct;
systemParams.memory = 7;
systemParams.deltaT = io_ref_struct.deltaT(1);
systemParams.t_MS2 = 1.4;
systemParams.rate_max = 1; % nothing faster than a second
systemParams.time_full = ((1:length(io_ref_struct.time_vec))*systemParams.deltaT)/60;% size(io_ref_struct.on_off_array,1);
systemParams.seq_length = length(systemParams.time_full);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%% Load cpHMM results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% designate simulation type
simTypeCell = {'out_only','in_only','koff_only_2','kon_only_2'};
for s = 1:length(simTypeCell)
    simType = simTypeCell{s};

    % specify 2 state network architecture (eventually this will be drawn from
    % actual fits)
    systemParams.R2 = [-.92  1/1.07; 
                        .92 -1/1.07]/60;
    % systemParams.r2 = [0 2.7]*1e5/60; % loading rate for each state

    % estimate r for now
    pon = systemParams.R2(2,1) / (systemParams.R2(2,1) + systemParams.R2(1,2));
    f_mean = nanmean(nanmean(io_ref_struct.fluo_raw_array(1:30,:),2),1);
    fluo_kernel = ms2_loading_coeff (systemParams.t_MS2, systemParams.memory);
    r2 = f_mean / sum(fluo_kernel) / pon;

    % systemParams.r2 = [0 4]*1e4; % loading rate for each state (in units of au per time step)
    systemParams.r2 = [0 r2];
    systemParams.pi0 = [0.5 0.5];
    if contains(simType,'_on') && ~contains(simType,'2')
        systemParams.pi0 = [0 0];
    end
    systemParams.noise = 5e3; % NL: this isn't really doing anything atm

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%% Call MCMC function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    close all force 

    % record "true" profile
    myCluster = parcluster('local');
    max_workers = myCluster.NumWorkers;
    sweepInfo.NumWorkers = min([24 max_workers]);
    sweepInfo.systemParams = systemParams;
    
    % generate ground truth vectors
    sweepInfo.p_on_true = nanmean(io_ref_struct.on_off_array,2);
    sweepInfo.fluo_true = nanmean(io_ref_struct.fluo_array,2);
    sweepInfo.fluo_true_raw = nanmean(io_ref_struct.fluo_raw_array,2);
    
    % save simulation type
    sweepInfo.simType = simType;
    
    % generate input TF profile array
    sweepInfo.tf_profile_array = io_ref_struct.knirps_array_norm;            
    
    % save time info
    sweepInfo.time_vec = io_ref_struct.time_vec;    
    sweepInfo.t_filter = io_ref_struct.t_filter;    
    sweepInfo.exp_time_bounds = io_ref_struct.time_bounds;
    
    % set list of parameters to sample 
    sweepInfo.paramList = {'HC','KD','F_min','K_out','K_in','k0'};
    sweepInfo.fitFlags = [1 1 0 1 1 0];
    sweepInfo.trueVals = [7,6e5,io_ref_struct.F_min_fit,1,1/60,0]; % NL: these are guesses...defaults that will be used if not flagged for fitting
    if contains(sweepInfo.simType,'2')      
        sweepInfo.fitFlags(strcmp(sweepInfo.paramList,'K_in')) = 0;
        sweepInfo.fitFlags(strcmp(sweepInfo.paramList,'K_out')) = 0;  
        if contains(sweepInfo.simType,'koff') 
            sweepInfo.fitFlags(strcmp(sweepInfo.paramList,'k0')) = 1;  
            sweepInfo.trueVals(strcmp(sweepInfo.paramList,'k0')) = systemParams.R2(1,2);  
        end
    end    
    sweepInfo.nFit = sum(sweepInfo.fitFlags);
    sweepInfo.nIterations = sweepInfo.nParamIncrement^sweepInfo.nFit;

    % set mcmc sampling hyperparameters
    sweepInfo.param_values = NaN(2,sweepInfo.nParamIncrement);
    sweepInfo.param_bounds(:,1) = linspace(0.5, 20, sweepInfo.nParamIncrement);
    sweepInfo.param_bounds(:,2) = linspace(1e5, 2e6, sweepInfo.nParamIncrement);
    sweepInfo.param_bounds(:,3) =  linspace(1e4, 1e5, sweepInfo.nParamIncrement);
    sweepInfo.param_bounds(:,4) = logspace(-3,0,sweepInfo.nParamIncrement);%[1e-2 60]';
    sweepInfo.param_bounds(:,5) = logspace(-3,0,sweepInfo.nParamIncrement);
    sweepInfo.param_bounds(:,6) = logspace(-3,0,sweepInfo.nParamIncrement);

    % initialize vectors to store results
    sweepInfo.param_fit_array = NaN(sweepInfo.nIterations,length(sweepInfo.paramList));
    iter = 1;
    for i = find(sweepInfo.fitFlags)
        vec_orig = sweepInfo.param_bounds(:,i);
        vec1 = repmat(vec_orig,sweepInfo.nIterations / sweepInfo.nParamIncrement^iter,1);
        vec2 = repelem(vec1,sweepInfo.nIterations  / sweepInfo.nParamIncrement^(sweepInfo.nFit-iter+1));
        sweepInfo.param_fit_array(:,i) = vec2;

        iter = iter + 1;
    end
    sweepInfo.param_fit_array(:,~sweepInfo.fitFlags) = repmat(sweepInfo.trueVals(~sweepInfo.fitFlags),size(sweepInfo.param_fit_array,1),1);

    % track profile
    sweepInfo.p_on_fit_array = NaN(sweepInfo.nIterations,length(sweepInfo.p_on_true));
    sweepInfo.fluo_fit_array = NaN(sweepInfo.nIterations,length(sweepInfo.p_on_true));

    % track function fit
    sweepInfo.objective_val_p_on = NaN(sweepInfo.nIterations,1);
    sweepInfo.objective_val_fluo = NaN(sweepInfo.nIterations,1);
    sweepInfo.objective_val_fluo_raw = NaN(sweepInfo.nIterations,1);

    % call parallel sweep script
    tic
    sweepTemp = sweep_par_loop(sweepInfo);
    toc

    % recombine
    sweepInfo.param_val_array = vertcat(sweepTemp.param_fit_array);
    sweepInfo.fluo_fit_array = [sweepTemp.fluo_fit_array];
    sweepInfo.p_on_fit_array = [sweepTemp.p_on_fit_array];
    sweepInfo.fluo_raw_fit_array = [sweepTemp.fluo_raw_fit_array];
    sweepInfo.objective_val_p_on = vertcat(sweepTemp.objective_val_p_on);
    sweepInfo.objective_val_fluo = vertcat(sweepTemp.objective_val_fluo);
    sweepInfo.objective_val_fluo_raw = vertcat(sweepTemp.objective_val_fluo_raw);
%     gillespie_struct = [sweepTemp.gillespie];

    save([resultsRoot 'sweepInfo_' simType '.mat'],'sweepInfo', '-v7.3');
%     save([resultsRoot 'gillespie_' simType '.mat'],'gillespie_struct', '-v7.3');
end