function sweepInfo = initializeSweepValues(sweepInfo)

  % set list of parameters to sample 
  F_min_guess = 3e4;
  sweepInfo.paramList = {'HC','KD','F_min','ks','ka','k0','kon','koff'};
  sweepInfo.fitFlags = [1 1 0 1 1 0 0 0];
  sweepInfo.trueVals = [7,6e5,F_min_guess,1,1,0,sweepInfo.R2(2,1),sweepInfo.R2(1,2)]; 

  if contains(sweepInfo.simType,'2')      
      sweepInfo.fitFlags(strcmp(sweepInfo.paramList,'ks')) = 0;
      sweepInfo.trueVals(strcmp(sweepInfo.paramList,'ks')) = 0;
      sweepInfo.fitFlags(strcmp(sweepInfo.paramList,'ka')) = 0;
      sweepInfo.trueVals(strcmp(sweepInfo.paramList,'ka')) = 0;
      if contains(sweepInfo.simType,'koff') 
          sweepInfo.fitFlags(strcmp(sweepInfo.paramList,'k0')) = 1;  
          sweepInfo.trueVals(strcmp(sweepInfo.paramList,'k0')) = sweepInfo.R2(1,2);  
      end
  end    
  sweepInfo.nFit = sum(sweepInfo.fitFlags);
  sweepInfo.nIterations = sweepInfo.nParamIncrement^sweepInfo.nFit;

  % set mcmc sampling hyperparameters
  sweepInfo.param_values = NaN(2,sweepInfo.nParamIncrement);
  % hill coefficient
  sweepInfo.param_bounds(:,1) = linspace(0.5, 20, sweepInfo.nParamIncrement);
  % KD
  sweepInfo.param_bounds(:,2) = linspace(1, 20, sweepInfo.nParamIncrement);
  % Spot detection threshold
  sweepInfo.param_bounds(:,3) =  linspace(1e4, 1e5, sweepInfo.nParamIncrement);    
  % constrain rates to same range
  sweepInfo.param_bounds(:,4) = logspace(-3,0,sweepInfo.nParamIncrement);%[1e-2 60]';
  sweepInfo.param_bounds(:,5) = logspace(-3,0,sweepInfo.nParamIncrement);
  sweepInfo.param_bounds(:,6) = logspace(-3,0,sweepInfo.nParamIncrement);
  
  sweepInfo.param_fit_array = NaN(sweepInfo.nIterations,length(sweepInfo.paramList));
  iter = 1;
  for i = find(sweepInfo.fitFlags)
      vec_orig = sweepInfo.param_bounds(:,i);
      vec1 = repmat(vec_orig,sweepInfo.nIterations / sweepInfo.nParamIncrement^iter,1);
      vec2 = repelem(vec1,sweepInfo.nIterations  / sweepInfo.nParamIncrement^(sweepInfo.nFit-iter+1));
      sweepInfo.param_fit_array(:,i) = vec2;

      iter = iter + 1;
  end
  sweepInfo.param_fit_array(:,~sweepInfo.fitFlags) = repmat(sweepInfo.trueVals(~sweepInfo.fitFlags),size(sweepInfo.param_fit_array,1),1);