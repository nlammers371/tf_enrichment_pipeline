function sweepTemp = sweep_par_loop(sweepInfo)

% temporalily switch to "slice-able" format
sweepTemp = struct;
for i = 1:sweepInfo.nIterations
    sweepTemp(i).param_fit_array = sweepInfo.param_fit_array(i,:);
    sweepTemp(i).p_on_fit_array = NaN(size(sweepInfo.p_on_fit_array(1,:)));
    sweepTemp(i).fluo_fit_array = NaN(size(sweepInfo.p_on_fit_array(1,:)));
    
%     sweepTemp(i).systemParams = sweepInfo.systemParams;
%     sweepTemp(i).paramList = sweepInfo.paramList;
%     sweepTemp(i).granularity = sweepInfo.granularity;
%     sweepTemp(i).n_traces = sweepInfo.n_traces;
%     sweepTemp(i).t_filter = sweepInfo.t_filter;
%     sweepTemp(i).tf_profile_array = sweepInfo.tf_profile_array;
%     sweepTemp(i).fluo_true = sweepInfo.fluo_true;
%     sweepTemp(i).fluo_true_raw = sweepInfo.fluo_true_raw;
%     sweepTemp(i).p_on_true = sweepInfo.p_on_true;
%     sweepTemp(i).simType = sweepInfo.simType;
end

% initialize array to track acceptance events

% NumWorkers = 24;    
pool = gcp('nocreate');
if isempty(pool)
  parpool(sweepInfo.NumWorkers);  
elseif  pool.NumWorkers ~= sweepInfo.NumWorkers      
  parpool(sweepInfo.NumWorkers);  
end  

h = waitbar(0,'conducting parameter sweeps...');
D = parallel.pool.DataQueue;    
afterEach(D, @nUpdateWaitbar);

N = sweepInfo.nIterations;
p = 1;

parfor sweep_step = 1:sweepInfo.nIterations
%     waitbar(sweep_step/sweepInfo.nIterations,wb);
    
    % update step
    sweepInfoTemp = sweepInfo;
    sweepInfoTemp.step = 1;
    sweepInfoTemp.param_fit_array = sweepInfo.param_fit_array(sweep_step,:);
%     sweepTemp(sweep_step).step = 1;%sweep_step;
    
    % conduct simulation
    simInfoPD = io_prediction_wrapper_v2(sweepInfoTemp);
    
    % store mean profiles    
    sweepTemp(sweep_step).p_on_fit_array = simInfoPD.p_on_array;
    sweepTemp(sweep_step).fluo_fit_array = simInfoPD.fluo_array;
    sweepTemp(sweep_step).fluo_raw_fit_array = simInfoPD.fluo_raw_array;
    
    % calculat error relative to fraction on profile
    ref_array_pon = sweepInfoTemp.p_on_true(sweepInfo.t_filter);    
    test_array_pon = simInfoPD.p_on_array;
    sweepTemp(sweep_step).objective_val_p_on = sum((ref_array_pon-test_array_pon).^2);    
    
    % calculat error relative to overall fluorescence profile
    ref_array_fluo = sweepInfoTemp.fluo_true(sweepInfo.t_filter);    
    test_array_fluo = simInfoPD.fluo_array;
    sweepTemp(sweep_step).objective_val_fluo = sum((ref_array_fluo-test_array_fluo).^2);  
    
    % calculat error relative to active spot fluorescence profile
    ref_array_fluo_raw = sweepInfoTemp.fluo_true_raw(sweepInfo.t_filter); 
    test_array_fluo_raw =  simInfoPD.fluo_raw_array;
    sweepTemp(sweep_step).objective_val_fluo_raw = sum((ref_array_fluo_raw-test_array_fluo_raw).^2);  
    
    % store select trace sim results
%     gill_small = rmfield(simInfoPD.gillespie,{'io_ref_out','initiation_rate_array','fluo_ms2_array_noise','fluo_ms2_array_noise_full'});
%     keep_indices = randsample(1:sweepInfo.n_traces,min([sweepInfo.n_keep sweepInfo.n_traces]),false);
%     gill_small.initiation_rate_array = gill_small.promoter_state_array(:,keep_indices);
%     gill_small.fluo_ms2_array = gill_small.fluo_ms2_array(:,keep_indices);
%     gill_small.initiation_rate_array = gill_small.fluo_ms2_array_full(:,keep_indices);
%     gill_small.output_response_ref = permute(gill_small.io_ref_in(:,1,keep_indices),[1 3 2]);
%     gill_small.output_rate_response_ref = permute(gill_small.rate_curve_in(:,1,keep_indices),[1 3 2]);
%     gill_small.input_tf_ref = permute(gill_small.tf_ref_in(:,1,keep_indices),[1 3 2]);    
%     sweepTemp(sweep_step).gillespie = gill_small;
    
    % update waitbar
    send(D, sweep_step);
end

delete(h);

% helper function
function nUpdateWaitbar(~)
  waitbar(p/N, h);
  p = p + 1;
end


end