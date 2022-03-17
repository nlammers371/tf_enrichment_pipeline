function simInfoPD = io_prediction_wrapper(mcmcInfo)

    paramList = mcmcInfo.paramList;
    HC = mcmcInfo.param_fit_array(mcmcInfo.step,strcmp(paramList,'HC'));
    KD = mcmcInfo.param_fit_array(mcmcInfo.step,strcmp(paramList,'KD'));
    K_out = mcmcInfo.param_fit_array(mcmcInfo.step,strcmp(paramList,'K_out'));
    K_in = mcmcInfo.param_fit_array(mcmcInfo.step,strcmp(paramList,'K_in'));
    F_min = mcmcInfo.param_fit_array(mcmcInfo.step,strcmp(paramList,'F_min'));
    
    tf_profile_array = mcmcInfo.tf_profile_array;
    n_traces = mcmcInfo.n_traces;
    granularity = mcmcInfo.granularity;
    
    % call stochastic simulation function
    simInfoPD = io_sim_function_v2(mcmcInfo.simType,mcmcInfo.systemParams,KD,HC,K_out,K_in,tf_profile_array,...
                                                                n_traces,granularity);

    % use output to generate predicted cumulative OFF (or ON) curve(s)
    simInfoPD.F_min = F_min;%logspace(3,log(5e4));
    simInfoPD = calculate_cumulative_dist(simInfoPD);