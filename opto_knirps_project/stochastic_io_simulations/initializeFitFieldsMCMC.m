function sweepInfo= initializeFitFieldsMCMC(sweepInfo)
  
%   sweepResults = struct;
  
  % set list of parameters to sample   
  sweepInfo.paramList = {'HC','KD','ks','ka','kon','koff'};
  sweepInfo.fitFlags = [1 1 1 1 1 1];
  sweepInfo.trueVals = [7,6e5,0,0,sweepInfo.R2_orig(2,1),sweepInfo.R2_orig(1,2)]; 
  
  if contains(sweepInfo.simType,'2')      
      ka_index = strcmp(sweepInfo.paramList,'ka');
      ks_index = strcmp(sweepInfo.paramList,'ks');
      sweepInfo.fitFlags([ks_index ka_index]) = 0; 
%       sweepInfo.trueVals([ks_index ka_index]) = 0; 
      if contains(sweepInfo.simType,'koff')
          sweepInfo.fitFlags(strcmp(sweepInfo.paramList,'kon')) = 0;
      elseif contains(sweepInfo.simType,'kon')
          sweepInfo.fitFlags(strcmp(sweepInfo.paramList,'koff')) = 0;
      end
  elseif strcmp(sweepInfo.simType,'match_exp')
      sweepInfo.fitFlags(1:end) = false;      
  else
      % koff is always fixed for 3 state models
      sweepInfo.fitFlags(strcmp(sweepInfo.paramList,'koff')) = 0;          
  end
  sweepInfo.nFit = sum(sweepInfo.fitFlags);
  if strcmp(sweepInfo.simType,'match_exp')
      sweepInfo.nIterations = 1;
      sweepInfo.nParamIncrement = 1;
  elseif ~isfield(sweepInfo, 'nIterations')
      sweepInfo.nIterations = sweepInfo.nParamIncrement^sweepInfo.nFit;
  end

  % set mcmc sampling hyperparameters
%   sweepInfo.param_values = NaN(2,sweepInfo.nParamIncrement);
  % hill coefficient
  sweepInfo.param_bounds(:,1) = linspace(0.5, 20, sweepInfo.nParamIncrement);
  % KD
  sweepInfo.param_bounds(:,2) = linspace(1, 20, sweepInfo.nParamIncrement); 
  % constrain rates to same range
  sweepInfo.param_bounds(:,3:6) = repmat(logspace(-3,0,sweepInfo.nParamIncrement),4,1)';%[1e-2 60]';    
% 
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