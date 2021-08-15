function sweepInfo = io_prediction_wrapper_wt(sweepInfo)

    paramList = sweepInfo.paramList;
    sweepInfo.HC = sweepInfo.param_fit_array(sweepInfo.step,strcmp(paramList,'HC'));
    sweepInfo.KD = sweepInfo.param_fit_array(sweepInfo.step,strcmp(paramList,'KD'));
    sweepInfo.ks = sweepInfo.param_fit_array(sweepInfo.step,strcmp(paramList,'ks'));
    sweepInfo.ka = sweepInfo.param_fit_array(sweepInfo.step,strcmp(paramList,'ka'));
    sweepInfo.k0 = sweepInfo.param_fit_array(sweepInfo.step,strcmp(paramList,'k0'));
    sweepInfo.F_min = sweepInfo.param_fit_array(sweepInfo.step,strcmp(paramList,'F_min'));   
        
%     % final model-building step
%     sweepInfo = generate_full_model(sweepInfo);                                                     
    
    % randomly draw tf profiles
    ap_axis_mean = sweepInfo.ap_axis_mean;
    dAP = ap_axis_mean(2)-ap_axis_mean(1);
    ap_bins = [ap_axis_mean-dAP/2 ap_axis_mean(end)+dAP/2]; 
    ap_profile_vec_tf = sweepInfo.ap_profile_vec_tf;
    ap_groups = discretize(ap_profile_vec_tf,ap_bins);
    
    n_traces_per_bin = ceil(sweepInfo.n_traces/length(ap_axis_mean));
    trace_id_vec = NaN(1,n_traces_per_bin*length(ap_axis_mean));
    trace_ap_vec = NaN(1,n_traces_per_bin*length(ap_axis_mean));
    for a = 1:length(ap_axis_mean)
        ind_ref = (a-1)*n_traces_per_bin+1:a*n_traces_per_bin;
        trace_id_vec(ind_ref) = randsample(find(ap_groups==a),n_traces_per_bin,true);
        trace_ap_vec(ind_ref) = a;
    end
            
    tf_profile_array = permute(sweepInfo.tf_profile_array_wt(:,trace_id_vec),[1 3 2]);
    
    %% call simulation function
    granularity_orig = sweepInfo.granularity;
    sweepInfo.granularity = sweepInfo.granularity*5;
    gillespie = synthetic_rate_gillespie_io_v3(sweepInfo,tf_profile_array);
    sweepInfo.granularity = granularity_orig;
    
    % use output to generate predicted curves
    fluo_array = gillespie.fluo_ms2_array;
    
    % instanaeous fraction ON    
    fluo_array_zeros = fluo_array;
    fluo_array_zeros(fluo_array_zeros<sweepInfo.F_min) = 0;    
    
    % calculate stats for fraction of traces that actually turn off
    index_vec = (1:size(fluo_array_zeros,1))';
    
    % calculate first and last active frames for each trace
    active_indices = 1*(fluo_array_zeros>0) .* index_vec;    
    active_indices(active_indices==0) = Inf;
    first_i_vec = min(active_indices);
    last_i_vec = max(active_indices.*~isinf(active_indices));
    all_off_flags = all(fluo_array_zeros==0);
    average_array = false(size(fluo_array_zeros));
    for a = find(~all_off_flags)      
        average_array(first_i_vec(a):last_i_vec(a),a) = true;  
    end
              
    %% calculate average fluorescence over trace lifetime
    off_time_vec_mean = NaN(size(ap_axis_mean));
    off_time_vec_ste = NaN(size(ap_axis_mean));
    fluo_vec_mean = NaN(size(ap_axis_mean));
    fluo_vec_ste = NaN(size(ap_axis_mean));
    for a = 1:length(ap_axis_mean)
        ap_ft = average_array;
        ap_ft(:,trace_ap_vec~=a) = false;
        
        mean_fluo = bootstrp(100,@(x)mean(x),fluo_array_zeros(ap_ft));
        fluo_vec_mean(a) = mean(mean_fluo);
        fluo_vec_ste(a) = std(mean_fluo);
        
        mean_off_times = bootstrp(100,@(x)mean(x),last_i_vec(trace_ap_vec==a));
        off_time_vec_mean(a) = mean(mean_off_times);
        off_time_vec_ste(a) = std(mean_off_times);
    end
                 

    sweepInfo.reactivation_time_cdf_predicted = ra_count_interp;
    sweepInfo.reactivation_time_axis_predicted = ra_time_vec;