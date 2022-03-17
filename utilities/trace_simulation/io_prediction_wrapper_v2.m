function simInfoPD = io_prediction_wrapper_v2(mcmcInfo)

    paramList = mcmcInfo.paramList;
    HC = mcmcInfo.param_fit_array(mcmcInfo.step,strcmp(paramList,'HC'));
    KD = mcmcInfo.param_fit_array(mcmcInfo.step,strcmp(paramList,'KD'));
    K_out = mcmcInfo.param_fit_array(mcmcInfo.step,strcmp(paramList,'K_out'));
    K_in = mcmcInfo.param_fit_array(mcmcInfo.step,strcmp(paramList,'K_in'));
    k0 = mcmcInfo.param_fit_array(mcmcInfo.step,strcmp(paramList,'k0'));
    F_min = mcmcInfo.param_fit_array(mcmcInfo.step,strcmp(paramList,'F_min'));
    
    tf_profile_array = mcmcInfo.tf_profile_array;
    n_traces = mcmcInfo.n_traces;
    granularity = mcmcInfo.granularity;
    
    % call stochastic simulation function
    simInfoPD = io_sim_function_v2(mcmcInfo.simType,mcmcInfo.systemParams,KD,HC,K_out,K_in,k0,tf_profile_array,...
                                                                n_traces,granularity);

    % use output to generate predicted curves
    fluo_array = simInfoPD.gillespie.fluo_ms2_array;
    
    % instanaeous fraction ON    
    fluo_array_zeros = fluo_array;
    fluo_array_zeros(fluo_array_zeros<F_min) = 0;
    simInfoPD.p_on_array = nanmean(fluo_array_zeros>0,2);
    
    % average overall fluorscence (including "OFF" frames)    
    simInfoPD.fluo_array = nanmean(fluo_array_zeros,2);
    
    % make raw version
    fluo_array_obs_only = simInfoPD.gillespie.fluo_ms2_array;
    fluo_array_obs_only(fluo_array_obs_only<F_min) = NaN;
    simInfoPD.fluo_array_obs_only = nanmean(fluo_array_obs_only,2);
    
    % calculate stats for fraction of traces that actually turn off
    
    
    
    simInfoPD.p_on_array = simInfoPD.p_on_array(mcmcInfo.t_filter,:);
    simInfoPD.fluo_array = simInfoPD.fluo_array(mcmcInfo.t_filter);
    simInfoPD.fluo_raw_array = simInfoPD.fluo_array_raw(mcmcInfo.t_filter);
