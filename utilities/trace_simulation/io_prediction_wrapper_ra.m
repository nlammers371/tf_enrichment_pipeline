function sweepResults = io_prediction_wrapper_ra(sweepInfo,sweepResults)

    paramList = sweepInfo.paramList;
    sweepInfo.HC = sweepResults.param_val_vec(strcmp(paramList,'HC'));
    sweepInfo.KD = sweepResults.param_val_vec(strcmp(paramList,'KD'));
    sweepInfo.ks = sweepResults.param_val_vec(strcmp(paramList,'ks'));
    sweepInfo.ka = sweepResults.param_val_vec(strcmp(paramList,'ka'));
    sweepInfo.kon = sweepResults.param_val_vec(strcmp(paramList,'kon'));
    sweepInfo.koff = sweepResults.param_val_vec(strcmp(paramList,'koff'));
              
    % final model-building step
    sweepInfo = generate_full_model(sweepInfo);

    % extract info
    n_traces = sweepInfo.n_traces;                                                      
%     sweepInfo.granularity = sweepInfo.granularity;
    
    % randomly draw tf profiles
    s_indices = randsample(1:size(sweepInfo.tf_profile_array_ra,2),n_traces,true);
    tf_profile_array = permute(sweepInfo.tf_profile_array_ra(:,s_indices),[1 3 2]);

    %% call simulation function
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
    perturbation_frame = ceil(size(tf_profile_array,1)/2);
    window_size = perturbation_frame-1;
    index_vec = -window_size:window_size;
    off_frames = sweepInfo.off_frame_ref;
    
    % was trace on before and after perturbation?
    active_indices = 1*(fluo_array_zeros>0) .* index_vec';
    first_i_vec = min(active_indices);
    max_time = max(sweepInfo.reactivation_time);
    before_flags = first_i_vec < 0;
    
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
    ra_time_vec = sweepInfo.reactivation_time;
             
    % generate count vec
    nBoots = 100;
    
    option_vec = 1:length(reactivation_time_vec);
    ra_array = NaN(nBoots,length(ra_time_vec));
    ra_array_full = NaN(nBoots,length(ra_time_vec));
    
    % conduct bootstrapping
    for n = 1:nBoots
        boot_indices = randsample(option_vec,length(option_vec),true);
        reactivation_time_vec_boot = reactivation_time_vec(boot_indices);
        ra_time_vec_boot = reactivation_time_vec_boot(~isnan(reactivation_time_vec_boot));
        ra_count_raw = (0:length(ra_time_vec_boot))/length(ra_time_vec_boot); 
         % sort   
        [ra_times_sorted,~] = sort(ra_time_vec_boot);
        % generate dummy time vector so that matlab won't throw a fit
        bs_time_vec = sort([0 ra_times_sorted+rand(size(ra_times_sorted))*1e-6]);      
        if length(bs_time_vec) > 1
            try
                if max(bs_time_vec) < ra_time_vec(end)
                    bs_time_vec(end+1) = ra_time_vec(end);
                    ra_count_raw(end+1) = ra_count_raw(end);
                end
                ra_array(n,:) = interp1(bs_time_vec,ra_count_raw,ra_time_vec,'previous');    

                ra_full = ra_array(n,:) + sum(isnan(reactivation_time_vec_boot))/sum(~isnan(reactivation_time_vec_boot));%(n_ever_on-length(option_vec))/length(option_vec);
                ra_array_full(n,:) = ra_full*sum(~isnan(reactivation_time_vec_boot))/length(reactivation_time_vec_boot);%length(option_vec) / n_ever_on;
            catch
                ra_array(n,:) = zeros(size(ra_time_vec));    
                ra_array_full(n,:) = zeros(size(ra_time_vec));  
            end
        else
            ra_array(n,:) = zeros(size(ra_time_vec));    
            ra_array_full(n,:) = ones(size(ra_time_vec));    
        end
    end
    
    % predicted stats
    ra_cdf_pd_mean = nanmean(ra_array);
%     ra_cdf_pd_ste = nanstd(ra_array)+1e-3;
    
    ra_cdf_full_mean = nanmean(ra_array_full);
%     ra_cdf_full_ste = nanstd(ra_array_full)+1e-3;
    
    % actual stats
    ra_cdf_exp_ste = sweepInfo.reactivation_cdf_ste;
    ra_cdf_exp_full_ste = sweepInfo.reactivation_cdf_full_ste;
    
    % calculate likelihood score
    delta_cdf = ra_cdf_pd_mean-sweepInfo.reactivation_cdf;
    logL_ra_cdf = (delta_cdf./ra_cdf_exp_ste).^2;
    sweepResults.ra_fit = -sqrt(mean(logL_ra_cdf));
    
    delta_cdf_full = ra_cdf_full_mean-sweepInfo.reactivation_cdf_full;
    logL_ra_full_cdf = -0.5*(delta_cdf_full./ra_cdf_exp_full_ste).^2;
    sweepResults.ra_full_fit = mean(logL_ra_full_cdf);
    
    % calculate simple R2 metrics
    delta_cdf = delta_cdf ./ mean(sweepInfo.reactivation_cdf);
    sweepResults.ra_fit_R2 = sum(delta_cdf.^2);
    delta_cdf_full = delta_cdf_full ./ mean(sweepInfo.reactivation_cdf_full);
    sweepResults.ra_full_fit_R2 = sum(delta_cdf_full.^2);
    
    if sweepInfo.keep_prediction_flag
        fluo_out = fluo_array_zeros;
        for i = 1:size(fluo_out,2)
            start_i = find(fluo_out(:,i)>0,1);
            stop_i = find(fluo_out(:,i)>0,1,'last');
            fluo_out(1:start_i-1,i) = NaN;
            fluo_out(stop_i+1:end,i) = NaN;
        end
    
        % save basic info
        sweepResults.ra_time_cdf_predicted = ra_cdf_pd_mean;
        sweepResults.ra_time_cdf_full_predicted = ra_cdf_full_mean;
        
        % save trace details
        sweepResults.ms2_traces_observed_ra = fluo_out;
        sweepResults.ms2_traces_true_ra = fluo_array;
        sweepResults.knirps_traces_ra = permute(tf_profile_array,[1 3 2]);
        sweepResults.reactivation_time_vec = reactivation_time_vec;
        sweepResults.tf_dependent_curves_ra = permute(gillespie.rate_curve_in,[1 3 2]); 
    end