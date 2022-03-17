function [sweepInfo, move_prop] = proposeMoveMCMC(sweepInfo)

    % get fit flags 
    fitFlags = sweepInfo.fitFlags==1;
    
    % get bounds
    paramBounds = 10.^sweepInfo.param_bounds;
    raw_params = sweepInfo.parameterFits(sweepInfo.sweep_step-1,:);
    curr_params = raw_params;
    
    sigmaVec = raw_params(fitFlags).*sweepInfo.moveSigma(fitFlags);
    
    % draw initial values from clipped gaussians
    lb_array = (paramBounds(1,fitFlags)-curr_params(fitFlags))./sigmaVec(fitFlags);
    ub_array = (paramBounds(2,fitFlags)-curr_params(fitFlags))./sigmaVec(fitFlags);
    
    % draw parameters
    move_prop = NaN(size(fitFlags));
    move_prop(fitFlags) = (curr_params(fitFlags)+sigmaVec(fitFlags).*trandn(lb_array,ub_array)');
    move_prop(~fitFlags) = curr_params(~fitFlags);