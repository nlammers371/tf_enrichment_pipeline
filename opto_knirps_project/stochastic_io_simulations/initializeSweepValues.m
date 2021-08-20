function sweepResults = initializeSweepValues(sweepInfo, sweepResults) 
  
  % initialize new array to store actual results  
  param_fit_array = NaN(sweepInfo.nIterations,length(sweepInfo.paramList));
  iter = 1;
  for i = find(sweepInfo.fitFlags)
      vec_orig = sweepInfo.param_bounds(:,i);
      vec1 = repmat(vec_orig,sweepInfo.nIterations / sweepInfo.nParamIncrement^iter,1);
      vec2 = repelem(vec1,sweepInfo.nIterations / sweepInfo.nParamIncrement^(sweepInfo.nFit-iter+1));     
      param_fit_array(:,i) = vec2;

      iter = iter + 1;
  end
  param_fit_array(:,~sweepInfo.fitFlags) = repmat(sweepInfo.trueVals(~sweepInfo.fitFlags),size(param_fit_array,1),1);
  % assign  to structure
  for j = 1:size(param_fit_array,1)
       sweepResults(j).param_val_vec = param_fit_array(j,:);
  end
  
 
  