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
simType = 'out_only_on';%'out_only_off';

% specify 2 state network architecture (eventually this will be drawn from
% actual fits)
systemParams = struct;
systemParams.R2 = [-2.5  1; 
                    2.5 -1]/60;
systemParams.r2 = [0 2.7]*1e5/60; % loading rate for each state
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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Call MCMC function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all force 

mcmcInfo = struct;

% record "true" profile
mcmcInfo.systemParams = systemParams;
mcmcInfo.p_on_true = nanmean(io_ref_struct.on_off_array,2);
mcmcInfo.simType = simType;
mcmcInfo.tf_profile_array_true = io_ref_struct.knirps_array_norm;
mcmcInfo.tf_profile_array = NaN(systemParams.seq_length,size(mcmcInfo.tf_profile_array_true,2));
start_i = find(systemParams.time_full<io_ref_struct.time_vec(1),1,'last');
mcmcInfo.tf_profile_array(1:start_i,:) = repmat(mcmcInfo.tf_profile_array_true(1,:),start_i,1);
mcmcInfo.tf_profile_array(start_i+1:end,:) = mcmcInfo.tf_profile_array_true;
mcmcInfo.time_vec = io_ref_struct.time_vec;
mcmcInfo.t_filter = ismember(systemParams.time_full,mcmcInfo.time_vec);
mcmcInfo.t_start = mcmcInfo.time_vec(1);
mcmcInfo.t_stop = mcmcInfo.time_vec(end);

% set basic parameters
mcmcInfo.nSteps = 1e2;
mcmcInfo.granularity = 10;
mcmcInfo.n_traces = 50;

% set list of parameters to sample 
mcmcInfo.paramList = {'HC','KD','F_min','K_out','K_in'};
mcmcInfo.fitFlags = [1 1 1 1 1];
mcmcInfo.trueVals = [6,6e5,3e4,1,1/60]; % NL: these are guesses

% set mcmc sampling hyperparameters
mcmcInfo.param_bounds = NaN(2,length(mcmcInfo.paramList));
mcmcInfo.param_bounds(:,1) = [0.5 10]';
mcmcInfo.param_bounds(:,2) = [1e5 1e6]';
mcmcInfo.param_bounds(:,3) = [1e4 1e5]';
mcmcInfo.param_bounds(:,4) = [1e-2 60]';
mcmcInfo.param_bounds(:,5) = [1e-2 60]';

% initialize vectors to store results
mcmcInfo.param_fit_array = NaN(mcmcInfo.nSteps,length(mcmcInfo.paramList));

% track profile
mcmcInfo.p_on_fit_array = NaN(mcmcInfo.nSteps,length(mcmcInfo.p_on_true));

% track function fit
mcmcInfo.objective_val = NaN(mcmcInfo.nSteps,1);

% initialize array to track acceptance events
mcmcInfo.param_move_array = NaN(mcmcInfo.nSteps,length(mcmcInfo.paramList));
wb = waitbar(0,'conducting mcmc inference...');
% conduct mcmc sampling
for mcmc_step = 1:mcmcInfo.nSteps
    waitbar(mcmc_step/mcmcInfo.nSteps,wb);
    
    mcmcInfo.step = mcmc_step;
    if mcmc_step~=1
        mcmcInfo.p_on_fit_array(mcmc_step,:) = mcmcInfo.p_on_fit_array(mcmc_step-1,:);
        mcmcInfo.param_fit_array(mcmc_step,:) = mcmcInfo.param_fit_array(mcmc_step-1,:);
        mcmcInfo.param_move_array(mcmc_step,:) = 0;
        mcmcInfo.objective_val(mcmc_step) = mcmcInfo.objective_val(mcmc_step-1);    
    end
    for p = 1:length(mcmcInfo.paramList)
        param_name = mcmcInfo.paramList{p};
        param_bounds = mcmcInfo.param_bounds(:,p);
        if mcmc_step == 1
            % if first step, draw samples using bounds
            if mcmcInfo.fitFlags(p)
                mcmcInfo.param_fit_array(mcmc_step,p) = rand()*diff(param_bounds)+param_bounds(1);          
            else
                mcmcInfo.param_fit_array(mcmc_step,p) = mcmcInfo.trueVals(p);
            end
            mcmcInfo.param_move_array(mcmc_step,p) = 1;            
        else            
            % extract current val
            param_val_curr = mcmcInfo.param_fit_array(mcmc_step,p);            
            % get proposal            
            if mcmcInfo.fitFlags(p)
                param_val_prop = make_mh_proposal(param_val_curr,param_bounds);           
                % generate a prediction
                mcmcInfo.param_fit_array(mcmc_step,p) = param_val_prop;
                simInfoPD = io_prediction_wrapper_v2(mcmcInfo);
                
                % perform MH step
                ref_array = mcmcInfo.p_on_true+1e-6;
                ref_array = ref_array/sum(ref_array);
                test_array = simInfoPD.p_on_array+1e-6;
                test_array = test_array/sum(test_array);
                objective_prop = sum((ref_array-test_array).^2);%sum(ref_array.*log(ref_array./test_array));
    
    %             objective_prop = sum((ref_array-test_array).^2);
%                 move_flag = exp(mcmcInfo.objective_val(mcmc_step)-objective_prop) > rand();
                move_flag = mcmcInfo.objective_val(mcmc_step)/objective_prop > rand();
                if move_flag
                    mcmcInfo.param_fit_array(mcmc_step,p) = param_val_prop;
                    mcmcInfo.param_move_array(mcmc_step,p) = 1;
                    mcmcInfo.objective_val(mcmc_step) = objective_prop;
                    mcmcInfo.p_on_fit_array(mcmc_step,:) = simInfoPD.p_on_array;
                end            
            end
            
            
        end        
        
    end
    
    if mcmc_step == 1
        % get first predicted profile
        simInfoPD = io_prediction_wrapper_v2(mcmcInfo);
        mcmcInfo.p_on_fit_array(1,:) = simInfoPD.p_on_array;

        ref_array = mcmcInfo.p_on_true+1e-6;
        test_array = simInfoPD.p_on_array+1e-6;
        objective_prop = sum((ref_array-test_array).^2);
%         objective_prop = sum(ref_array.*log(ref_array./test_array));
        mcmcInfo.objective_val(mcmc_step) = objective_prop;
    end
end

delete(wb);