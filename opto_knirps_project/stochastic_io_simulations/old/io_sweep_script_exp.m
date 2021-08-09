% Script to scan through parameter space for selected biological variables
clear
close all

addpath(genpath('../utilities'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%% Load experimental data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

projectName = 'optokni_eve4+6_ON'; 

liveProject = LiveEnrichmentProject(projectName);
resultsRoot = [liveProject.dataPath filesep];

% load data
load([resultsRoot 'io_ref_struct.mat'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%% Load cpHMM results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% designate simulation type
simType = 'in_only_off';%'out_only_off';

% specify 2 state network architecture (eventually this will be drawn from
% actual fits)
systemParams = struct;
systemParams.R2 = [-.92  1/1.07; 
                    .92 -1/1.07]/60;
% systemParams.r2 = [0 2.7]*1e5/60; % loading rate for each state
systemParams.r2 = [0 4]*1e4; % loading rate for each state
systemParams.pi0 = [0.5 0.5];
if contains(simType,'on')
    systemParams.pi0 = [0 0];
end
systemParams.noise = 5e3; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify key system characteristics

systemParams.memory = 7;
systemParams.deltaT = io_ref_struct.deltaT;
systemParams.t_MS2 = 1.4;
systemParams.time_full = (systemParams.deltaT:systemParams.deltaT:io_ref_struct.time_vec(end)*60)/60;% size(io_ref_struct.on_off_array,1);
systemParams.seq_length = length(systemParams.time_full);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Call MCMC function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all force 

sweepInfo = struct;

% record "true" profile
sweepInfo.systemParams = systemParams;
sweepInfo.p_on_true = nanmean(io_ref_struct.on_off_array,2);
sweepInfo.fluo_true = nanmean(io_ref_struct.fluo_array,2);
sweepInfo.simType = simType;
sweepInfo.tf_profile_array_true = io_ref_struct.knirps_array_norm;
sweepInfo.tf_profile_array = NaN(systemParams.seq_length,size(sweepInfo.tf_profile_array_true,2));
start_i = find(systemParams.time_full<io_ref_struct.time_vec(1),1,'last');
sweepInfo.tf_profile_array(1:start_i,:) = repmat(sweepInfo.tf_profile_array_true(1,:),start_i,1);
sweepInfo.tf_profile_array(start_i+1:end,:) = sweepInfo.tf_profile_array_true;
sweepInfo.time_vec = io_ref_struct.time_vec;
sweepInfo.t_filter = ismember(systemParams.time_full,sweepInfo.time_vec);
sweepInfo.t_start = sweepInfo.time_vec(1);
sweepInfo.t_stop = sweepInfo.time_vec(end);

% set basic parameters
sweepInfo.nParamIncrement = 1e1;
sweepInfo.granularity = 5;
sweepInfo.n_traces = 100;
sweepInfo.n_keep = 10;

% set list of parameters to sample 
sweepInfo.paramList = {'HC','KD','F_min','K_out','K_in'};
sweepInfo.fitFlags = [1 1 0 1 1];
sweepInfo.trueVals = [7,6e5,3e4,1,1/60]; % NL: these are guesses
sweepInfo.nFit = sum(sweepInfo.fitFlags);
sweepInfo.nIterations = sweepInfo.nParamIncrement^sweepInfo.nFit;

% set mcmc sampling hyperparameters
sweepInfo.param_values = NaN(2,sweepInfo.nParamIncrement);
sweepInfo.param_bounds(:,1) = linspace(0.5, 10, sweepInfo.nParamIncrement);
sweepInfo.param_bounds(:,2) = linspace(1e5, 1e6, sweepInfo.nParamIncrement);
sweepInfo.param_bounds(:,3) =  linspace(1e4, 1e5, sweepInfo.nParamIncrement);
sweepInfo.param_bounds(:,4) = logspace(-3,0,sweepInfo.nParamIncrement);%[1e-2 60]';
sweepInfo.param_bounds(:,5) = logspace(-3,0,sweepInfo.nParamIncrement);

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
sweepInfo.param_fit_array(:,~sweepInfo.fitFlags) = sweepInfo.trueVals(~sweepInfo.fitFlags);

% track profile
sweepInfo.p_on_fit_array = NaN(sweepInfo.nIterations,length(sweepInfo.p_on_true));
sweepInfo.fluo_fit_array = NaN(sweepInfo.nIterations,length(sweepInfo.p_on_true));

% track function fit
sweepInfo.objective_val_p_on = NaN(sweepInfo.nIterations,1);
sweepInfo.objective_val_fluo = NaN(sweepInfo.nIterations,1);

% % initialize temporary structure to permit parallelization
% sweepTemp = struct;
% for i = 1:sweepInfo.nIterations
%     sweepTemp(i).param_fit_array = sweepInfo.param_fit_array(i,:);
%     sweepTemp(i).p_on_fit_array = NaN(size(sweepInfo.p_on_fit_array(1,:)));
%     sweepTemp(i).fluo_fit_array = NaN(size(sweepInfo.p_on_fit_array(1,:)));
% end

% initialize array to track acceptance events
wb = waitbar(0,'conducting mcmc inference...');

% conduct mcmc sampling
for sweep_step = 1:sweepInfo.nIterations
    waitbar(sweep_step/sweepInfo.nIterations,wb);
    
    % update step
    sweepInfo.step = sweep_step;
    
    % cpnduct simulation
    simInfoPD = io_prediction_wrapper_v2(sweepInfo);
    
    % store mean profiles    
    sweepInfo.p_on_fit_array(sweep_step,:) = simInfoPD.p_on_array;
    sweepInfo.fluo_fit_array(sweep_step,:) = simInfoPD.fluo_array;
    
    % calculat error
    ref_array_pon = sweepInfo.p_on_true;    
    test_array_pon = simInfoPD.p_on_array;
    sweepInfo.objective_val_p_on(sweep_step) = sum((ref_array_pon-test_array_pon).^2);    
    
    ref_array_fluo = sweepInfo.fluo_true;    
    test_array_fluo = simInfoPD.fluo_array;
    sweepInfo.objective_val_fluo(sweep_step) = sum((ref_array_fluo-test_array_fluo).^2);  
    
    % store select trace sim results
    gill_small = rmfield(simInfoPD.gillespie,{'io_ref_out','initiation_rate_array','fluo_ms2_array_noise','fluo_ms2_array_noise_full'});
    keep_indices = randsample(1:sweepInfo.n_traces,min([sweepInfo.n_keep sweepInfo.n_traces]),false);
    gill_small.initiation_rate_array = gill_small.promoter_state_array(:,keep_indices);
    gill_small.fluo_ms2_array = gill_small.fluo_ms2_array(:,keep_indices);
    gill_small.initiation_rate_array = gill_small.fluo_ms2_array_full(:,keep_indices);
    gill_small.io_ref_in = permute(gill_small.io_ref_in(:,1,keep_indices),[1 3 2]);
    sweepInfo.gillespie{sweep_step} = gill_small;
end
delete(wb);

save([resultsRoot 'sweepInfo_' simType '.mat'],'sweepInfo', '-v7.3');