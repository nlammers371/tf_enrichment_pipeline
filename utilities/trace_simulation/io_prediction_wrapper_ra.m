function simInfoPD = io_prediction_wrapper_ra(sweepInfo)

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
    granularity = sweepInfo.granularity;                                                          
    
    % randomly draw tf profiles
    s_indices = randsample(1:size(sweepInfo.tf_profile_array_ra,2),n_traces,true);
    tf_profile_array = permute(sweepInfo.tf_profile_array_ra(:,s_indices),[1 3 2]);

    %% call simulation function
    gillespie = synthetic_rate_gillespie_io_v3(simInfo);
    
    % use output to generate predicted curves
    fluo_array = gillespie.fluo_ms2_array;
    
    % instanaeous fraction ON    
    fluo_array_zeros = fluo_array;
    fluo_array_zeros(fluo_array_zeros<sweepInfo.F_min) = 0;    
    
    % calculate stats for fraction of traces that actually turn off
    perturbation_frame = ceil(sweepInfo.seq_length/2);
    window_size = perturbation_frame-1;
    index_vec = -window_size:window_size;
    off_frames = sweepInfo.off_frame_ref;
    
    % was trace on before and after perturbation?
    active_indices = 1*(fluo_array_zeros>0) .* index_vec';
    first_i_vec = min(active_indices);
    last_i_vec = max(active_indices);
    before_flags = first_i_vec < 0;
    after_flags = last_i_vec > 0;
    
    % was trace OFF immediately preceding ON perturbation?
    off_flags = all(fluo_array_zeros(ismember(index_vec,off_frames),:)==0);

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