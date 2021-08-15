function sweepInfo = io_prediction_wrapper_ra(sweepInfo)

    paramList = sweepInfo.paramList;
    sweepInfo.HC = sweepInfo.param_fit_array(sweepInfo.step,strcmp(paramList,'HC'));
    sweepInfo.KD = sweepInfo.param_fit_array(sweepInfo.step,strcmp(paramList,'KD'));
    sweepInfo.ks = sweepInfo.param_fit_array(sweepInfo.step,strcmp(paramList,'ks'));
    sweepInfo.ka = sweepInfo.param_fit_array(sweepInfo.step,strcmp(paramList,'ka'));
    sweepInfo.k0 = sweepInfo.param_fit_array(sweepInfo.step,strcmp(paramList,'k0'));
    sweepInfo.F_min = sweepInfo.param_fit_array(sweepInfo.step,strcmp(paramList,'F_min'));   
        
    % final model-building step
    sweepInfo = generate_full_model(sweepInfo);

    % extract info
    n_traces = sweepInfo.n_traces;
%     granularity = sweepInfo.granularity;                                                          
    
    % randomly draw tf profiles
    s_indices = randsample(1:size(sweepInfo.tf_profile_array_ra,2),n_traces,true);
    tf_profile_array = permute(sweepInfo.tf_profile_array_ra(:,s_indices),[1 3 2]);

    %% call simulation function
    gillespie = synthetic_rate_gillespie_io_v3(sweepInfo,tf_profile_array);
    
    % use output to generate predicted curves
    fluo_array = gillespie.fluo_ms2_array;
    
    % instanaeous fraction ON    
    fluo_array_zeros = fluo_array;
    fluo_array_zeros(fluo_array_zeros<sweepInfo.F_min) = 0;    
    
    % calculate stats for fraction of traces that actually turn off
    perturbation_frame = ceil(size(tf_profile_array,1)/2);
    window_size = perturbation_frame-1;
    index_vec = -window_size:window_size;
    off_frames = sweepInfo.off_frame_ref;
    
    % was trace on before and after perturbation?
    active_indices = 1*(fluo_array_zeros>0) .* index_vec';
    first_i_vec = min(active_indices);
%     last_i_vec = max(active_indices);
    max_time = max(sweepInfo.reactivation_time);
    before_flags = first_i_vec < 0;
%     after_flags = last_i_vec > 0;
    
    % was trace OFF immediately preceding ON perturbation?
    off_flags = all(fluo_array_zeros(ismember(index_vec,off_frames),:)==0);

    % ra_flags 
    ra_indices = find(off_flags&before_flags);    
    reactivation_time_vec = NaN(1,size(fluo_array_zeros,2));
    
    % set all leading and trailing zeros to NaN
    for i = ra_indices        
        fluo_vec = fluo_array_zeros(:,i);

        % find last "OFF" frame
        temp = fluo_vec;
        temp(index_vec<0) = 0;        
        last_i = find(temp>0,1);        
        if isempty(last_i)
            reactivation_time_vec(i) = max_time+100;
        else
            switch_time = index_vec(last_i);
            % estimate reactivation time                    
            reactivation_time_vec(i) = switch_time*sweepInfo.deltaT;  
        end                                                                                             
    end
            
    % construct empirical cdf for ractivation
    ra_times = reactivation_time_vec(~isnan(reactivation_time_vec));           
    ra_time_vec = sweepInfo.reactivation_time;
    
    % sort   
    [ra_times_sorted,~] = sort(ra_times);  
    
    % generate count vec
    ra_count_raw = (0:length(ra_times))/length(ra_times); 
    
    % generate dummy time vector so that matlab won't throw a fit
    bs_time_vec = sort([0 ra_times_sorted+rand(size(ra_times_sorted))*1e-6]);      

    ra_count_interp = interp1(bs_time_vec,ra_count_raw,ra_time_vec,'previous');                  

    sweepInfo.reactivation_time_cdf_predicted = ra_count_interp;
    sweepInfo.reactivation_time_axis_predicted = ra_time_vec;