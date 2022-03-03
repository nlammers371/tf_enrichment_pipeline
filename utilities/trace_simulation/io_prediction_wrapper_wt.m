function sweepResults = io_prediction_wrapper_wt(sweepInfo,sweepResults)

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
    ap_profile_vec_tf = sweepInfo.ap_profile_vec_tf;
    ap_filter = ap_profile_vec_tf<=ap_bins(2)&ap_profile_vec_tf>=ap_bins(1);    
    trace_id_vec = randsample(find(ap_filter),sweepInfo.n_traces,true);                
    tf_profile_array = permute(sweepInfo.tf_profile_array_wt(:,trace_id_vec),[1 3 2]);
    
    %% call simulation function    
    % generate random inputs for tf and parameters    
    gillespie = synthetic_rate_gillespie_io_v3(sweepInfo,tf_profile_array);        
    
    % use output to generate predicted curves
    fluo_array = gillespie.fluo_ms2_array;
    
    % instanaeous fraction ON    
    fluo_array_zeros = fluo_array;
    
    % use reference curve to determine which time points are detected and
    % which are not
    if isfield(sweepInfo,'fluo_ref_curve')
        df = sweepInfo.fluo_ref_curve(2) - sweepInfo.fluo_ref_curve(1);
        ref_vec = [sweepInfo.fluo_ref_curve-df/2 sweepInfo.fluo_ref_curve(end)+df];
        fluo_groups = discretize(fluo_array_zeros,ref_vec);
        miss_probs = sweepInfo.p_miss_ref_vec (fluo_groups);
        miss_status = miss_probs > rand(size(miss_probs));
        fluo_array_zeros(miss_status) = 0;    
    else
        fluo_array_zeros = ones(size(fluo_array));
    end
       
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
    still_on_vec = still_on_array(:);
    knirps_array = permute(tf_profile_array,[1 3 2]);
    knirps_groups = discretize(knirps_array(:),sweepInfo.knirps_bins_still_on);
    predicted_fraction_still_on = zeros(size(sweepInfo.knirps_axis_still_on));
    for k = 1:length(predicted_fraction_still_on)
        val = nanmean(still_on_vec(knirps_groups==k));
        if ~isnan(val)
            predicted_fraction_still_on(k) = val;
        end
    end        
    predicted_fraction_still_on(1) = 1;    
    %% calculate average fluorescence over trace lifetime 
    fluo_nan = fluo_array_zeros;
    still_on_flags = still_on_array;
    still_on_flags(isnan(still_on_flags)) = 0;
    fluo_nan(~still_on_flags) = NaN;    
    time_filter = ismember(round(sweepInfo.time_axis_wt/60,2),round(sweepInfo.time_axis_mf,2));
    mean_fluo_time_predicted = nanmean(fluo_nan(time_filter,:),2);
    mean_fluo_time_predicted(isnan(mean_fluo_time_predicted)) = 0;      
           
    % calcualte differences
    delta_fluo_time = mean_fluo_time_predicted-sweepInfo.fluo_time_mean;
%     delta_fluo_time = delta_fluo_time ./ nanmean(sweepInfo.fluo_time_mean);

    delta_still_on = predicted_fraction_still_on - sweepInfo.fraction_still_on';
%     delta_still_on = delta_still_on ./ nanmean(sweepInfo.fraction_still_on);
    
    % calculate simple RMS differences          
    sweepResults.fluo_time_fit_R2 = mean(delta_fluo_time.^2);    
    sweepResults.still_on_fit_R2 = mean(delta_still_on.^2);
    
    % calculate log likelihood of experimental trends assuming gaussian
    % errors
    logL_fluo_time = (delta_fluo_time./sweepInfo.fluo_time_ste).^2;% 
    sweepResults.fluo_time_fit = -sqrt(mean(logL_fluo_time));
    
    logL_still_on = (delta_still_on./sweepInfo.fraction_still_on_ste').^2;% + log(2*pi*fluo_vec_ste.^2));
    sweepResults.still_on_fit = -sqrt(mean(logL_still_on));
    
    % save extra simulation info if option is flagged
    if sweepInfo.keep_prediction_flag
      
        fluo_out = fluo_array_zeros;
        fluo_out(~average_array) = NaN;
        
        % basic AP trends
        sweepResults.fluo_time_predicted = mean_fluo_time_predicted;
        sweepResults.p_still_on_predicted = predicted_fraction_still_on;
        
        % traces
        sweepResults.ms2_traces_observed_wt = fluo_out;
        sweepResults.ms2_traces_true_wt = fluo_array;
        sweepResults.knirps_traces_wt = permute(gillespie.tf_ref_in,[1 3 2]);        
        sweepResults.tf_dependent_curves_wt = permute(gillespie.rate_curve_in,[1 3 2]);        
    end