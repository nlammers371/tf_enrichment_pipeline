function logL_vec = calculate_bootstrap_logL(mcmcInfo,simInfoPD)

    selection_array = floor(rand(mcmcInfo.nBoots,mcmcInfo.n_traces)*mcmcInfo.n_traces)*size(simInfoPD.gillespie.fluo_ms2_array,1) + mcmcInfo.ref_array;    
    trace_samples = simInfoPD.gillespie.fluo_ms2_array(selection_array);
    pon_samples = trace_samples >= simInfoPD.F_min;
    pon_mean = permute(nanmean(nanmean(pon_samples,2),1),[1 3 2]);
    pon_ste = permute(nanstd(nanmean(pon_samples,2),[],1),[1 3 2]);
    
    % calculate likelihood of experimental data given prediction
    ref_array_pon = mcmcInfo.p_on_true; 
    logL_vec = -0.5*((ref_array_pon'-pon_mean)./pon_ste).^2 - 0.5*log(2*pi*pon_ste.^2);
    simInfoPD.pon_mean = pon_mean;
    simInfoPD.pon_ste = pon_ste;