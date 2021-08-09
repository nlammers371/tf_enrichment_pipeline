function sweepTempFull = sweep_par_loop_v2(sweepInfo,savePath)

% make directory to store temporary files
tempSavePath = [savePath filesep sweepInfo.simType '_tempSweepFiles' filesep];
mkdir(tempSavePath)

% initialize array to track acceptance events

% NumWorkers = 24;    
pool = gcp('nocreate');
if isempty(pool)
  parpool(sweepInfo.NumWorkers);  
elseif  pool.NumWorkers ~= sweepInfo.NumWorkers     
  delete(pool)
  parpool(sweepInfo.NumWorkers);  
end  

h = waitbar(0,'conducting parameter sweeps...');
D = parallel.pool.DataQueue;    
afterEach(D, @nUpdateWaitbar);

N = sweepInfo.nIterations;
p = 1;
tic
for sweep_step = 1:sweepInfo.nIterations
    
%     waitbar(sweep_step/sweepInfo.nIterations,wb);
    
    % update step
    sweepInfoTemp = sweepInfo;
    sweepInfoTemp.step = 1;
    sweepInfoTemp.param_fit_array = sweepInfo.param_fit_array(sweep_step,:);
%     sweepTemp(sweep_step).step = 1;%sweep_step;
    
    % conduct simulation
    simInfoPD = io_prediction_wrapper_v3(sweepInfoTemp);
    
    sweepTemp = struct;
    % store mean profiles    
    sweepTemp.p_on_fit_array = simInfoPD.p_on_array;
    sweepTemp.fluo_fit_array = simInfoPD.fluo_array;
    sweepTemp.fluo_raw_fit_array = simInfoPD.fluo_array_raw;
    sweepTemp.fluo_obs_only_array = simInfoPD.fluo_array_obs_only;
    sweepTemp.off_fluo_array = simInfoPD.off_fluo_array;
    sweepTemp.reactivation_time_vec = simInfoPD.reactivation_time_vec;
    sweepTemp.mean_ra = mean(~isnan(simInfoPD.reactivation_time_vec));
    sweepTemp.reactivation_time_cdf = simInfoPD.reactivation_time_cdf;
    sweepTemp.reactivation_time_axis = simInfoPD.reactivation_time_axis;

    sweepTemp.tf_dependent_rate = nanmean(simInfoPD.gillespie.rate_curve_in,3);
    
    % cave temporary file
    saveTempResults(sweep_step,sweepTemp,tempSavePath);
    
    % update waitbar
    send(D, sweep_step);
    
end
toc
% delete pool
delete(gcp) 

delete 
delete(h);

sweepTempFull = reassembleTempResults(tempSavePath,sweepInfo.nIterations);

% helper function
function nUpdateWaitbar(~)
  waitbar(p/N, h);
  p = p + 1;
end


end