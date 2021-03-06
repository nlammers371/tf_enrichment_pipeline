% Script to scan through parameter space for selected biological variables
clear
close all

addpath(genpath('../utilities'));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Simulate "true" data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% designate simulation type
simType = 'out_only_on';%'out_only_off';

% specify 2 state network architecture (eventually this will be drawn from
% actual fits)
systemParams = struct;
systemParams.R2 = [-2  1; 
                    2 -1]/60;
systemParams.r2 = [0 1]*1e4; % loading rate for each state
systemParams.pi0 = [0.5 0.5];
if contains(simType,'on')
    systemParams.pi0 = [0 0];
end
systemParams.noise = 5e3; 

% designate key variables to sweep
KD = 3e4;
HC = 5;
K_out = 1/60;%linspace(1,60,10);
if contains(simType,'out')
    K_out = 60;
end
K_in = 1/60;
% F_min = 1e4;
%%
% call stochastic simulation function
tic
simInfoTrue = io_sim_function(simType,systemParams,KD,HC,K_out(end),K_in,...
                                          'n_traces',1000,'granularity',1);
toc

% use output to generate predicted cumulative OFF (or ON) curve(s)
simInfoTrue.F_min = 5e3;%logspace(3,log(5e4));
simInfoTrue = calculate_cumulative_dist(simInfoTrue);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Call MCMC function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all force 

mcmcInfo = struct;
% record "true" profile
mcmcInfo.systemParams = systemParams;
mcmcInfo.p_on_true = simInfoTrue.p_on_array;
mcmcInfo.simType = simType;
mcmcInfo.trueParams = simInfoTrue;
mcmcInfo.tf_profile_array = simInfoTrue.tf_profile_array;

% set basic parameters
mcmcInfo.nSteps = 5e2;
mcmcInfo.granularity = 10;
mcmcInfo.n_traces = 50;

% set list of parameters to sample 
mcmcInfo.paramList = {'HC','KD','F_min','K_out','K_in'};
mcmcInfo.fitFlags = [1 1 1 0 1];
mcmcInfo.trueVals = [HC,KD,F_min,K_out,K_in];

% set mcmc sampling hyperparameters
mcmcInfo.param_bounds = NaN(2,length(mcmcInfo.paramList));
mcmcInfo.param_bounds(:,1) = [0.5 10]';
mcmcInfo.param_bounds(:,2) = [1e4 1e5]';
mcmcInfo.param_bounds(:,3) = [1e3 1e4]';
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
                simInfoPD = io_prediction_wrapper(mcmcInfo);
                % perform MH step
                ref_array = diff(mcmcInfo.p_on_true)+1e-6;
                test_array = diff(simInfoPD.p_on_array)+1e-6;
                objective_prop = sum(ref_array.*log(ref_array./test_array));
    %             objective_prop = sum((ref_array-test_array).^2);
                move_flag = exp((mcmcInfo.objective_val(mcmc_step))-objective_prop) > rand();
    %             move_flag = mcmcInfo.objective_val(mcmc_step)/objective_prop > rand();
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
        simInfoPD = io_prediction_wrapper(mcmcInfo);
        mcmcInfo.p_on_fit_array(1,:) = simInfoPD.p_on_array;

        ref_array = diff(mcmcInfo.p_on_true)+1e-6;
        test_array = diff(simInfoPD.p_on_array)+1e-6;
%             objective_prop = sum((ref_array-test_array).^2);
        objective_prop = sum(ref_array.*log(ref_array./test_array));
        mcmcInfo.objective_val(mcmc_step) = objective_prop;
    end
end

delete(wb);