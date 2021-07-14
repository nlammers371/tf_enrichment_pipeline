function sweepTemp = sweep_par_loop(sweepInfo)


sweepTemp = struct;
for i = 1:sweepInfo.nIterations
    sweepTemp(i).param_fit_array = sweepInfo.param_fit_array(i,:);
    sweepTemp(i).p_on_fit_array = NaN(size(sweepInfo.p_on_fit_array(1,:)));
    sweepTemp(i).fluo_fit_array = NaN(size(sweepInfo.p_on_fit_array(1,:)));
    sweepTemp(i).systemParams = sweepInfo.systemParams;
    sweepTemp(i).paramList = sweepInfo.paramList;
    sweepTemp(i).granularity = sweepInfo.granularity;
    sweepTemp(i).n_traces = sweepInfo.n_traces;
    sweepTemp(i).t_filter = sweepInfo.t_filter;
    sweepTemp(i).tf_profile_array = sweepInfo.tf_profile_array;
    sweepTemp(i).fluo_true = sweepInfo.fluo_true;
    sweepTemp(i).p_on_true = sweepInfo.p_on_true;
    sweepTemp(i).simType = sweepInfo.simType;
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
tic
parfor sweep_step = 1:sweepInfo.nIterations
%     waitbar(sweep_step/sweepInfo.nIterations,wb);
    
    % update step
    sweepTemp(sweep_step).step = 1;%sweep_step;
    
    % cpnduct simulation
    simInfoPD = io_prediction_wrapper_v2(sweepTemp(sweep_step));
    
    % store mean profiles    
    sweepTemp(sweep_step).p_on_fit_array = simInfoPD.p_on_array;
    sweepTemp(sweep_step).fluo_fit_array = simInfoPD.fluo_array;
    
    % calculat error
    ref_array_pon = sweepTemp(sweep_step).p_on_true;    
    test_array_pon = simInfoPD.p_on_array;
    sweepTemp(sweep_step).objective_val_p_on = sum((ref_array_pon-test_array_pon).^2);    
    
    ref_array_fluo = sweepTemp(sweep_step).fluo_true;    
    test_array_fluo = simInfoPD.fluo_array;
    sweepTemp(sweep_step).objective_val_fluo = sum((ref_array_fluo-test_array_fluo).^2);  
    
    % store select trace sim results
    gill_small = rmfield(simInfoPD.gillespie,{'io_ref_out','initiation_rate_array','fluo_ms2_array_noise','fluo_ms2_array_noise_full'});
    keep_indices = randsample(1:sweepInfo.n_traces,min([sweepInfo.n_keep sweepInfo.n_traces]),false);
    gill_small.initiation_rate_array = gill_small.promoter_state_array(:,keep_indices);
    gill_small.fluo_ms2_array = gill_small.fluo_ms2_array(:,keep_indices);
    gill_small.initiation_rate_array = gill_small.fluo_ms2_array_full(:,keep_indices);
    gill_small.io_ref_in = permute(gill_small.io_ref_in(:,1,keep_indices),[1 3 2]);
    sweepTemp(sweep_step).gillespie = gill_small;
    
    % update waitbar
    send(D, sweep_step);
end
toc
delete(h);

% helper function
function nUpdateWaitbar(~)
  waitbar(p/N, h);
  p = p + 1;
end


end