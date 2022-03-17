function sweepInfo= initializeFitFieldsMCMC(sweepInfo)
  
%   sweepResults = struct;
  
  % set list of parameters to sample   
  sweepInfo.paramList = {'HC','KD','ks','ka','kon','koff'};
  sweepInfo.fitFlags = [1 1 1 1 1 1];
  sweepInfo.priorMean = [6,5,0.1,0.1,sweepInfo.R2_orig(2,1),sweepInfo.R2_orig(1,2)];
  sweepInfo.priorSigma = 10.^[1,log10(3),1,1,1,1];
  sweepInfo.moveSigma = [.01,.01,.01,.01,.01,.01];
  
  if contains(sweepInfo.simType,'2')      
      ka_index = strcmp(sweepInfo.paramList,'ka');
      ks_index = strcmp(sweepInfo.paramList,'ks');
      sweepInfo.fitFlags([ks_index ka_index]) = 0; 
      sweepInfo.priorMean([ks_index ka_index]) = 0; 
      sweepInfo.priorSigma([ks_index ka_index]) = 0; 
      if contains(sweepInfo.simType,'koff')
          sweepInfo.fitFlags(strcmp(sweepInfo.paramList,'kon')) = 0;
%           sweepInfo.priorSigma([ks_index ka_index]) = 0; 
      elseif contains(sweepInfo.simType,'kon')
          sweepInfo.fitFlags(strcmp(sweepInfo.paramList,'koff')) = 0;
      end   
  else
      % koff is always fixed for 3 state models
      sweepInfo.fitFlags(strcmp(sweepInfo.paramList,'koff')) = 0;          
  end
  sweepInfo.nFit = sum(sweepInfo.fitFlags);
  
 
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  % Set parameter bounds
  
  % hill coefficient
  sweepInfo.param_bounds(:,1) = log10([0.5, 20]);
  
  % KD
  sweepInfo.param_bounds(:,2) = [0, 2]; 
  
  % constrain rates to same range
  sweepInfo.param_bounds(:,3:6) = repmat([-3 0],4,1)';%[1e-2 60]';    
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % initialize arrays to store results
  sweepInfo.parameterFits = NaN(sweepInfo.n_iters_max,length(sweepInfo.paramList));
  sweepInfo.ra_fit = NaN(sweepInfo.n_iters_max,1);
  sweepInfo.fluo_time_fit = NaN(sweepInfo.n_iters_max,1);
  sweepInfo.still_on_fit = NaN(sweepInfo.n_iters_max,1);
  sweepInfo.overall_fit = NaN(sweepInfo.n_iters_max,1);
  
  
%   % track fits to observables of interest
%   for i = 1:sweepInfo.nIterations
%       sweepResults(i).ra_fit = NaN;%(sweepInfo.nIterations,1);
%       sweepResults(i).ra_full_fit = NaN;
%       sweepResults(i).fluo_time_fit = NaN;%(sweepInfo.nIterations,1);
%       sweepResults(i).still_on_fit = NaN;%(sweepInfo.nIterations,1);
%             
%       sweepResults(i).ra_fit_R2 = NaN;
%       sweepResults(i).ra_full_fit_R2 = NaN;
%       sweepResults(i).fluo_time_fit_R2 = NaN;
%       sweepResults(i).still_on_fit_R2 = NaN;
%       
%       if sweepInfo.keep_prediction_flag
%           sweepResults(i).ra_time_cdf_predicted = NaN;%(1,length(sweepInfo.reactivation_time));
%           sweepResults(i).ra_time_cdf_full_predicted = NaN;%(1,length(sweepInfo.reactivation_time));
% 
%           sweepResults(i).p_still_on_predicted = NaN;
%           sweepResults(i).fluo_time_predicted = NaN;%(1,length(sweepInfo.ap_axis_mean));          
%           
%           sweepResults(i).ms2_traces_observed_wt = NaN;
%           sweepResults(i).ms2_traces_true_wt = NaN;
%           sweepResults(i).knirps_traces_wt = NaN;          
%           sweepResults(i).tf_dependent_curves_wt = NaN;
%           
%           sweepResults(i).ms2_traces_observed_ra = NaN;
%           sweepResults(i).ms2_traces_true_ra = NaN;
%           sweepResults(i).knirps_traces_ra = NaN;
%           sweepResults(i).reactivation_time_vec = NaN;
%           sweepResults(i).tf_dependent_curves_ra = NaN; 
%           
%       end
%   end