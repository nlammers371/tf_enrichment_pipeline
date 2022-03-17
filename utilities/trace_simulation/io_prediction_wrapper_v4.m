function simInfoPD = io_prediction_wrapper_v3(seepInfo)

    paramList = seepInfo.paramList;
    HC = seepInfo.param_fit_array(seepInfo.step,strcmp(paramList,'HC'));
    KD = seepInfo.param_fit_array(seepInfo.step,strcmp(paramList,'KD'));
    K_out = seepInfo.param_fit_array(seepInfo.step,strcmp(paramList,'ks'));
    K_in = seepInfo.param_fit_array(seepInfo.step,strcmp(paramList,'ka'));
    k0 = seepInfo.param_fit_array(seepInfo.step,strcmp(paramList,'k0'));
    F_min = seepInfo.param_fit_array(seepInfo.step,strcmp(paramList,'F_min'));
    
    tf_profile_array = seepInfo.tf_profile_array;
    n_traces = seepInfo.n_traces;
    granularity = seepInfo.granularity;
    
    % call stochastic simulation function
    simInfoPD = io_sim_function_v2(seepInfo.simType,seepInfo.systemParams,KD,HC,K_out,K_in,k0,tf_profile_array,...
                                                                n_traces,granularity);

                                                              
    

    s_indices = randsample(1:size(tf_profile_array,2),simInfo.n_traces,true);
    simInfo.tf_profile_array = permute(tf_profile_array(:,s_indices),[1 3 2]);

    %% call simulation function
    simInfo.gillespie = synthetic_rate_gillespie_io_v2(simInfo);
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
    perturbation_frame = ceil(simInfoPD.seq_length/2);
    window_size = perturbation_frame-1;
    index_vec = -window_size:window_size;
    off_frames = -6:-1;
    
    % was trace on before and after perturbation?
    active_indices = 1*(fluo_array_zeros>0) .* index_vec';
    first_i_vec = min(active_indices);
    last_i_vec = max(active_indices);
    before_flags = first_i_vec < 0;
    after_flags = last_i_vec > 0;
    
    % was trace OFF immediately preceding ON perturbation?
    off_flags = all(fluo_array_zeros(ismember(index_vec,off_frames),:)==0);
        
    fluo_array_raw = fluo_array_zeros;
    all_zero_flags = all(fluo_array_raw==0);
    fluo_array_raw(:,all_zero_flags) = NaN;
    
    off_fluo_frames = -30:-1;
    off_fluo_array = NaN(length(off_fluo_frames),size(fluo_array_zeros,2));
    reactivation_time_vec = NaN(1,size(fluo_array_zeros,2));
    
    % set all leading and trailing zeros to NaN
    for i = 1:size(fluo_array_raw,2)
        fluo_vec_raw = fluo_array_raw(:,i);
        start_i = find(fluo_vec_raw>0,1);
        fluo_vec_raw(1:start_i-1) = NaN;
        fluo_vec_raw(find(fluo_vec_raw>0,1,'last')+1:end) = NaN;
        
        if before_flags(i) && off_flags(i)         
        
            first_i = find(fluo_vec_raw(index_vec<0)>0,1,'last')+1;
            temp = fluo_vec_raw;
            temp(index_vec<0) = 0;
            last_i = find(temp>0,1)-1;
            if isempty(last_i)
                last_i = length(fluo_vec_raw);
            end  

            fluo_vec_raw(first_i:last_i) = NaN;                 

            % store "off fluo"
            off_frames = index_vec(start_i:first_i) - index_vec(first_i) - 1;
            to_i = ismember(off_fluo_frames,off_frames);
            from_i = ismember(off_frames,off_fluo_frames);
            off_fluo_array(to_i,i) = fluo_vec_raw(from_i);  

            if after_flags(i) % require that trace turns back on for reactivation times                                                              
                % estimate reactivation time            
                switch_time = min(index_vec(~isnan(fluo_vec_raw') & index_vec>=0));
                reactivation_time_vec(i) = switch_time*simInfoPD.deltaT;
            end
        end
        
        fluo_array_raw(:,i) = fluo_vec_raw;   
    end
    
    % store additional results
    simInfoPD.reactivation_time_vec = reactivation_time_vec';
    simInfoPD.fluo_array_raw = nanmean(fluo_array_raw,2);
    simInfoPD.off_fluo_array = nanmean(off_fluo_array,2);
    
    % construct empirical cdf for ractivation
    ra_times = simInfoPD.reactivation_time_vec(~isnan(simInfoPD.reactivation_time_vec))';
%     max_ra_time = round(max(ra_times)*simInfoPD.deltaT);   
    max_time = 50*60;
    ra_time_vec = 0:simInfoPD.deltaT:max_time;
%     ra_count_interp = NaN(size(ra_time_vec));
    if length(ra_times) > 20
      
      [ra_times_sorted,~] = sort(ra_times);      
      ra_count_raw = (0:length(ra_times))/length(ra_times);          
      bs_time_vec = [0 ra_times_sorted+rand(size(ra_times_sorted))*1e-6 max_time];      
      [bs_time_sorted,~] = sort(bs_time_vec);
      ra_count_interp = interp1(bs_time_sorted,[ra_count_raw 1],ra_time_vec);            
      
    else      
      ra_count_interp = NaN(size(ra_time_vec));
    end
    simInfoPD.reactivation_time_cdf = ra_count_interp;
    simInfoPD.reactivation_time_axis = ra_time_vec;