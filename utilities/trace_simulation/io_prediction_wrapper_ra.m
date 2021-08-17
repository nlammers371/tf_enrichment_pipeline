function sweepResults = io_prediction_wrapper_ra(sweepInfo,sweepResults)

    paramList = sweepInfo.paramList;
    sweepInfo.HC = sweepResults.param_val_vec(strcmp(paramList,'HC'));
    sweepInfo.KD = sweepResults.param_val_vec(strcmp(paramList,'KD'));
    sweepInfo.ks = sweepResults.param_val_vec(strcmp(paramList,'ks'));
    sweepInfo.ka = sweepResults.param_val_vec(strcmp(paramList,'ka'));
    sweepInfo.k0 = sweepResults.param_val_vec(strcmp(paramList,'k0'));
    sweepInfo.F_min = sweepResults.param_val_vec(strcmp(paramList,'F_min'));   
        
    % final model-building step
    sweepInfo = generate_full_model(sweepInfo);

    % extract info
    n_traces = sweepInfo.n_traces;                                                      
    
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
    max_time = max(sweepInfo.reactivation_time);
    before_flags = first_i_vec < 0;
    n_ever_on = sum(max(active_indices)>0);
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
    nBoots = 100;
    ra_count_raw = (0:length(ra_times))/length(ra_times); 
    option_vec = 1:length(ra_times_sorted);
    ra_array = NaN(nBoots,length(ra_time_vec));
    ra_array_full = NaN(nBoots,length(ra_time_vec));
    % conduct bootstrapping
    for n = 1:nBoots
        boot_indices = randsample(option_vec,length(option_vec),true);
        % generate dummy time vector so that matlab won't throw a fit
        bs_time_vec = sort([0 ra_times_sorted(boot_indices)+rand(size(ra_times_sorted))*1e-6]);      
        if length(bs_time_vec) > 1
            ra_array(n,:) = interp1(bs_time_vec,ra_count_raw,ra_time_vec,'previous');          
            ra_full = ra_array(n,:) + (n_ever_on-length(option_vec))/length(option_vec);
            ra_array_full(n,:) = ra_full * length(option_vec) / n_ever_on;
        else
            ra_array(n,:) = zeros(size(ra_time_vec));                  
        end
    end
    ra_cdf_mean = nanmean(ra_array);
    ra_cdf_ste = nanstd(ra_array)+1e-3;
    
    ra_cdf_full_mean = nanmean(ra_array_full);
    ra_cdf_full_ste = nanstd(ra_array_full)+1e-3;
    
    % calculate likelihood score
    logL_ra_cdf = -0.5*(((ra_cdf_mean-sweepInfo.reactivation_cdf)./ra_cdf_ste).^2);% + log(2*pi*ra_cdf_ste.^2));
    sweepResults.ra_fit = mean(logL_ra_cdf);
    
    logL_ra_full_cdf = -0.5*(((ra_cdf_full_mean-sweepInfo.reactivation_cdf_full)./ra_cdf_full_ste).^2);% + log(2*pi*ra_cdf_ste.^2));
    sweepResults.ra_full_fit = mean(logL_ra_full_cdf);
    
    if sweepInfo.keep_prediction_flag
        sweepResults.ra_time_cdf_predicted = ra_cdf_mean;
        sweepResults.ra_time_cdf_full_predicted = ra_cdf_full_mean;
        sweepResults.tf_dependent_curve_ra = nanmean(gillespie.rate_curve_in,3)';
    end