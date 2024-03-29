function sweepResults = io_prediction_wrapper_ON(sweepInfo,sweepResults)

    paramList = sweepInfo.paramList;
    sweepInfo.HC = sweepResults.param_val_vec(strcmp(paramList,'HC'));
    sweepInfo.KD = sweepResults.param_val_vec(strcmp(paramList,'KD'));
    sweepInfo.ks = sweepResults.param_val_vec(strcmp(paramList,'ks'));
    sweepInfo.ka = sweepResults.param_val_vec(strcmp(paramList,'ka'));      
    sweepInfo.kon = sweepResults.param_val_vec(strcmp(paramList,'kon'));
    sweepInfo.koff = sweepResults.param_val_vec(strcmp(paramList,'koff'));
    
    % final model-building step
    sweepInfo = generate_full_model(sweepInfo);                                                     
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % randomly draw tf profiles   
    ap_bins = [sweepInfo.ap_limits_still_on(1) sweepInfo.ap_limits_still_on(2)];
    ap_profile_vec_tf = sweepInfo.ap_profile_vec_tf_ON;
    ap_filter = ap_profile_vec_tf<=ap_bins(2)&ap_profile_vec_tf>=ap_bins(1);    
    trace_id_vec = randsample(find(ap_filter),sweepInfo.n_traces,true);                
    tf_profile_array = permute(sweepInfo.tf_profile_array_ON(:,trace_id_vec),[1 3 2]);
    
    %% call simulation function    
    gillespie = synthetic_rate_gillespie_io_v3(sweepInfo,tf_profile_array);        
    
    % use output to generate predicted curves
    fluo_array = gillespie.fluo_ms2_array;
    
    % instanaeous fraction ON    
    fluo_array_zeros = fluo_array;
    
    % use reference curve to determine which time points are detected and
    % which are not    
    df = sweepInfo.fluo_ref_curve(2) - sweepInfo.fluo_ref_curve(1);
    ref_vec = [sweepInfo.fluo_ref_curve-df/2 sweepInfo.fluo_ref_curve(end)+df];
    fluo_groups = discretize(fluo_array_zeros,ref_vec);
    miss_probs = sweepInfo.p_miss_ref_vec (fluo_groups);
    miss_status = miss_probs > rand(size(miss_probs));
    fluo_array_zeros(miss_status) = 0;    

       
    % calculate stats for fraction of traces that actually turn off
    index_vec = (1:size(fluo_array_zeros,1))';
    
    % calculate first and last active frames for each trace
    active_indices = 1*(fluo_array_zeros>0) .* index_vec;    
    active_indices2 = active_indices;
    active_indices2(active_indices2==0) = Inf;
    first_i_vec = min(active_indices2);
    last_i_vec = max(active_indices);
    all_off_flags = all(fluo_array_zeros==0);
    average_array = false(size(fluo_array_zeros));
    for a = find(~all_off_flags)      
        average_array(first_i_vec(a):last_i_vec(a),a) = true;  
    end
    
    %% calculate fraction off response curve             
    still_on_array = NaN(size(fluo_array_zeros));
    for i = 1:size(fluo_array_zeros,2)
        if ~isinf(first_i_vec(i))
            still_on_array(first_i_vec(i):end,i) = 1;
            still_on_array(last_i_vec(i):end,i) = 0;
        end
    end         
       
    %% calculate average fluorescence over trace lifetime 
    fluo_nan = fluo_array_zeros;
    still_on_flags = still_on_array;
    still_on_flags(isnan(still_on_flags)) = 0;
    fluo_nan(~still_on_flags) = NaN;    
    time_filter = ismember(round(sweepInfo.time_axis_wt/60,2),round(sweepInfo.time_axis_mf,2));
    mean_fluo_time_predicted = nanmean(fluo_nan(time_filter,:),2);
    mean_fluo_time_predicted(isnan(mean_fluo_time_predicted)) = 0;      
           
    % calcualte differences
    delta_fluo_time = mean_fluo_time_predicted-sweepInfo.fluo_time_mean_ON;
    
    % calculate simple RMS differences          
    sweepResults.fluo_time_ON_fit_R2 = mean(delta_fluo_time.^2);        
    
    % calculate log likelihood of experimental trends assuming gaussian
    % errors
    logL_fluo_time = (delta_fluo_time./sweepInfo.fluo_time_ste).^2;% 
    sweepResults.fluo_time_ON_fit = -sqrt(mean(logL_fluo_time));  
    
    % save extra simulation info if option is flagged
    if sweepInfo.keep_prediction_flag
      
        fluo_out = fluo_array_zeros;
        fluo_out(~average_array) = NaN;
        
        % basic AP trends
        sweepResults.fluo_time_predicted_ON = mean_fluo_time_predicted;        
        
        % traces
        sweepResults.ms2_traces_observed_ON = fluo_out;
        sweepResults.ms2_traces_true_ON = fluo_array;
        sweepResults.knirps_traces_ON = sweepInfo.tf_profile_array_wt(:,trace_id_vec);        
        sweepResults.tf_dependent_curves_ON = permute(gillespie.rate_curve_in,[1 3 2]);        
    end