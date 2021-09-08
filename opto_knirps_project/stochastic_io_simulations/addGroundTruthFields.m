function sweepInfo = addGroundTruthFields(sweepInfo,io_ref_ra,io_ref_wt)

    % RA fields. Take only times up to designated point
    if ~isempty(io_ref_ra)
        use_flags = io_ref_ra.reactivation_time_axis<=sweepInfo.max_ra_time;
        sweepInfo.reactivation_cdf = io_ref_ra.reactivation_time_cdf(use_flags);
        sweepInfo.reactivation_cdf_ste = io_ref_ra.reactivation_time_cdf_ste(use_flags);
        sweepInfo.reactivation_cdf_full = io_ref_ra.reactivation_time_cdf_full(use_flags);
        sweepInfo.reactivation_cdf_full_ste = io_ref_ra.reactivation_time_cdf_ste_full(use_flags);
        sweepInfo.reactivation_time = io_ref_ra.reactivation_time_axis(use_flags);

        sweepInfo.off_frame_ref = io_ref_ra.off_frame_ref;
    end
    
    % detection limit info
    sweepInfo.detection_limit = io_ref_wt.F_min_fit;
    if isfield(io_ref_wt, 'fluo_ref_curve')
        sweepInfo.fluo_ref_curve  = io_ref_wt.fluo_ref_curve;
        sweepInfo.p_miss_ref_vec  = io_ref_wt.p_miss_ref_vec;
    end
    
    % fluo trend over time
    sweepInfo.time_axis_mf = io_ref_wt.time_axis_mf;
    sweepInfo.fluo_time_mean = io_ref_wt.fluo_time_mean;
    sweepInfo.fluo_time_ste = io_ref_wt.fluo_time_ste;
    sweepInfo.time_axis_wt = io_ref_wt.time_axis;
    
    if isfield(io_ref_wt,'knirps_axis_still_on')      
        sweepInfo.knirps_axis_still_on = io_ref_wt.knirps_axis_still_on;
        sweepInfo.knirps_bins_still_on = io_ref_wt.knirps_bins_still_on;
        sweepInfo.fraction_still_on_ste = io_ref_wt.fraction_still_on_ste;
        sweepInfo.fraction_still_on = io_ref_wt.fraction_still_on_mean;
        sweepInfo.ap_limits_still_on = io_ref_wt.ap_limits_still_on;
    end
    if sweepInfo.calculate_ap_metrics
        % mean fluorescence vs. AP
        sweepInfo.mean_fluo_ap = io_ref_wt.fluo_vec_mean;
        sweepInfo.mean_fluo_ap_ste = io_ref_wt.fluo_vec_ste;

        % observed off times
        sweepInfo.off_time_ap = io_ref_wt.off_time_vec_mean;
        sweepInfo.off_time_ap_ste = io_ref_wt.off_time_vec_ste;
        sweepInfo.ap_axis_mean = io_ref_wt.ap_axis_mean;               
    end
    
    % save TF profiles and time vec for RA type
    if ~isempty(io_ref_ra)
        sweepInfo.tf_profile_array_ra = io_ref_ra.knirps_array;               
        sweepInfo.time_axis_ra = io_ref_ra.time_vec';
    end
    % save TF profiles and time vec for WT type
    % calculate a sensible start time
    start_time = ceil(nanmean(io_ref_wt.on_time_vec) / 60)*60;
    [~,start_i] = min(abs(io_ref_wt.time_axis-start_time));
    
    sweepInfo.tf_profile_array_wt = io_ref_wt.knirps_array(start_i:end,:);               
    sweepInfo.wt_start_time = start_time;
    sweepInfo.ap_profile_vec_tf = io_ref_wt.mean_ap;               
    sweepInfo.time_axis_wt = io_ref_wt.time_axis(start_i:end)';