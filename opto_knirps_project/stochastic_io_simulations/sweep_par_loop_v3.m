function sweepResults = sweep_par_loop_v3(sweepInfo,sweepResults,savePath)

% make directory to store temporary files
tempSavePath = [savePath filesep sweepInfo.simType '_tempSweepFiles' filesep];
mkdir(tempSavePath)

% initialize parallel pools
initializePool(sweepInfo)

% initialize stuff for waitbar
WB = waitbar(0,'conducting parameter sweeps...');
D = parallel.pool.DataQueue;    
afterEach(D, @nUpdateWaitbar);

N = sweepInfo.nIterations;
p = 1;
    
% iterate through different param values
nIterations = sweepInfo.nIterations;
parfor sweep_step = 1:nIterations
%     waitbar(sweep_step/nIterations,WB);                
    
    % conduct RA simulations
    sweepResults(sweep_step) = io_prediction_wrapper_ra(sweepInfo,sweepResults(sweep_step));
    
    % conduct WT simulations
    sweepResults(sweep_step) = io_prediction_wrapper_wt(sweepInfo,sweepResults(sweep_step));
         
    % update waitbar
    send(D, sweep_step);
end
toc
% delete pool (necessary to clear from RAM)
delete(gcp) 
 
delete(WB);

% sweepTempFull = reassembleTempResults(tempSavePath,sweepInfo.nIterations);

% helper function
function nUpdateWaitbar(~)
  waitbar(p/N, WB);
  p = p + 1;
end


end