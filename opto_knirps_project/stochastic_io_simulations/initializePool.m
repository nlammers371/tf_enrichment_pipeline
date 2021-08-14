function initializePool(sweepInfo)

    pool = gcp('nocreate');
    if isempty(pool)
      parpool(sweepInfo.NumWorkers);  
    elseif  pool.NumWorkers ~= sweepInfo.NumWorkers     
      delete(pool)
      parpool(sweepInfo.NumWorkers);  
    end  
   