function sweepInfo = run_mcmc_sampling(sweepInfo)

% make directory to store temporary files
% tempSavePath = [savePath filesep sweepInfo.simType '_tempSweepFiles' filesep];
% mkdir(tempSavePath)

% initialize parallel pools
if sweepInfo.n_chains > 5
    initializePool(sweepInfo)
end

% initialize stuff for waitbar
WB = waitbar(0,'conducting parameter sweeps...');
% D = parallel.pool.DataQueue;    
% afterEach(D, @nUpdateWaitbar);
% 
% N = sweepInfo.nIterations;
% p = 1;
    
% iterate through different param values
for chain = 1
    for sweep_step = 1:n_iters_max%sweepInfo.n_iters_max%sweepInfo.n_iters_max
        % increment waitbar
        waitbar(sweep_step/sweepInfo.n_iters_max,WB);
        sweepInfo.sweep_step = sweep_step;
        
        % first, we need to generate new values to test
        if sweep_step==1
            [sweepInfo, prop_params] = initialGuessMCMC(sweepInfo);
        else
            [sweepInfo, prop_params] = proposeMoveMCMC(sweepInfo);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Check goodness of fit
        sweepTemp = struct;
        sweepTemp.param_val_vec = prop_params;
        
        % conduct RA simulations
        sweepTemp = io_prediction_wrapper_ra(sweepInfo,sweepTemp);
        
        % conduct WT simulations
        sweepTemp = io_prediction_wrapper_wt(sweepInfo,sweepTemp);

        % initialize move flag        
        move_flag = false;
        % compare likelihoods
        if sweep_step > 1
            prev_logL = sweepInfo.overall_fit(sweep_step-1);
            prop_logL = sweepTemp.still_on_fit + sweepTemp.ra_fit + sweepTemp.fluo_time_fit;
            move_flag = exp(prop_logL-prev_logL) >= rand();
        end
        
        % update waitbar
        if sweep_step == 1 || move_flag
            sweepInfo.parameterFits(sweep_step,:) = sweepTemp.param_val_vec;
            sweepInfo.ra_fit(sweep_step) = sweepTemp.ra_fit;
            sweepInfo.fluo_time_fit(sweep_step) = sweepTemp.fluo_time_fit;
            sweepInfo.still_on_fit(sweep_step) = sweepTemp.still_on_fit;
            sweepInfo.overall_fit(sweep_step) = sweepTemp.still_on_fit + sweepTemp.ra_fit + sweepTemp.fluo_time_fit;
        else
            sweepInfo.parameterFits(sweep_step,:) = sweepInfo.parameterFits(sweep_step-1,:);
            sweepInfo.ra_fit(sweep_step) = sweepInfo.ra_fit(sweep_step-1);
            sweepInfo.fluo_time_fit(sweep_step) = sweepInfo.fluo_time_fit(sweep_step-1);
            sweepInfo.still_on_fit(sweep_step) = sweepInfo.still_on_fit(sweep_step-1);
            sweepInfo.overall_fit(sweep_step) = sweepInfo.still_on_fit(sweep_step-1) + sweepInfo.ra_fit(sweep_step-1) + sweepInfo.fluo_time_fit(sweep_step-1);
        end
    end
end
% delete pool (necessary to clear from RAM)
% delete(gcp) 
 
delete(WB);

% sweepTempFull = reassembleTempResults(tempSavePath,sweepInfo.nIterations);

% helper function
% function nUpdateWaitbar(~)
%   waitbar(p/N, WB);
%   p = p + 1;
% end


end