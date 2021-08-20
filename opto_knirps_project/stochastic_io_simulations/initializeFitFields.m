function [sweepInfo, sweepResults] = initializeFitFields(sweepInfo,sweepResults)
  
%   sweepResults = struct;
  
  % set list of parameters to sample   
  sweepInfo.paramList = {'HC','KD','F_min','ks','ka','k0','kon','koff'};
  sweepInfo.fitFlags = [1 1 0 1 1 0 0 0];
  sweepInfo.trueVals = [7,6e5,sweepInfo.detection_limit,1,1,0,sweepInfo.R2(2,1),sweepInfo.R2(1,2)]; 
  
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
  if strcmp(sweepInfo.simType,'match_exp')
      sweepInfo.nIterations = 1;
      sweepInfo.nParamIncrement = 1;
  elseif ~isfield(sweepInfo, 'nIterations')
      sweepInfo.nIterations = sweepInfo.nParamIncrement^sweepInfo.nFit;
  end

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
  
 % track fits to observables of interest
  for i = 1:sweepInfo.nIterations
      sweepResults(i).ra_fit = NaN;%(sweepInfo.nIterations,1);
      sweepResults(i).ra_full_fit = NaN;
      sweepResults(i).mean_fluo_fit = NaN;%(sweepInfo.nIterations,1);
      sweepResults(i).off_time_fit = NaN;%(sweepInfo.nIterations,1);
      
      
      sweepResults(i).ra_fit_R2 = NaN;
      sweepResults(i).ra_full_fit_R2 = NaN;
      sweepResults(i).mean_fluo_fit_R2 = NaN;
      sweepResults(i).off_time_fit_R2 = NaN;
      sweepResults(i).pon_fit_R2 = NaN;
      
      if sweepInfo.keep_prediction_flag
          sweepResults(i).ra_time_cdf_predicted = NaN;%(1,length(sweepInfo.reactivation_time));
          sweepResults(i).ra_time_cdf_full_predicted = NaN;%(1,length(sweepInfo.reactivation_time));

          sweepResults(i).p_still_on_predicted = NaN;
          sweepResults(i).mean_fluo_predicted = NaN;%(1,length(sweepInfo.ap_axis_mean));
          sweepResults(i).off_time_predicted = NaN;%(1,length(sweepInfo.ap_axis_mean));
          
          sweepResults(i).ms2_traces_observed_wt = NaN;
          sweepResults(i).ms2_traces_true_wt = NaN;
          sweepResults(i).knirps_traces_wt = NaN;
          sweepResults(i).trace_ap_vec_wt = NaN;
          sweepResults(i).tf_dependent_curves_wt = NaN;
          
          sweepResults(i).ms2_traces_observed_ra = NaN;
          sweepResults(i).ms2_traces_true_ra = NaN;
          sweepResults(i).knirps_traces_ra = NaN;
          sweepResults(i).reactivation_time_vec = NaN;
          sweepResults(i).tf_dependent_curves_ra = NaN; 
          
      end
  end