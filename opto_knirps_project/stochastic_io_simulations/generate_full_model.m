% Function wrapper for stochastic trace simulations in non-steady-state
% conditions

function sweepInfo = generate_full_model(sweepInfo,varargin)

% If we're looking at a model where koff is the concentration-dependent parameter,
% we need to ensure that the effective ON rate is consistent with cpHMM estimates
% if strcmp(sweepInfo.simType,'koff_only_2')     
%     R2_orig = sweepInfo.R2_orig;
%     kon_true = R2_orig(2,1)*(sweepInfo.ka + sweepInfo.ks)/sweepInfo.ka;
%     R2 = R2_orig;
%     R2(2,1) = kon_true;
%     R2(1,1) = -kon_true;
%     sweepInfo.R2 = R2;
% end    

% Build 2 state network
sweepInfo.RateMatrix = zeros(3,3);
sweepInfo.RateMatrix(2,3) = sweepInfo.koff;
sweepInfo.RateMatrix(3,2) = sweepInfo.kon;

% Add third, silent state
sweepInfo.RateMatrix(:,1) = [0 ; sweepInfo.ka ; 0];
sweepInfo.RateMatrix(1,:) = [0 sweepInfo.ks 0];

% normalize
diag_flags = eye(size(sweepInfo.RateMatrix,1))==1;
sweepInfo.RateMatrix(diag_flags) = 0;
sweepInfo.RateMatrix(diag_flags) = -sum(sweepInfo.RateMatrix);

% update other parameters
sweepInfo.r_emission = [0 sweepInfo.r2]; % loading rate for each state
sweepInfo.pi0 = [0 sweepInfo.pi0];

% indicate which rate is tf-dependent (assume only one possible for now)
sweepInfo.tf_dependent_flags = false(3,3);
if strcmp(sweepInfo.simType,'match_exp') 
    % do nothing
elseif contains(sweepInfo.simType,'out')      
    sweepInfo.tf_dependent_flags(1,2) = true;
elseif contains(sweepInfo.simType,'in')      
    sweepInfo.HC = -sweepInfo.HC;
    sweepInfo.tf_dependent_flags(2,1) = true;
elseif contains(sweepInfo.simType,'kon')      
    sweepInfo.HC = -sweepInfo.HC;
    sweepInfo.tf_dependent_flags(3,2) = true;
elseif contains(sweepInfo.simType,'koff')          
    sweepInfo.tf_dependent_flags(2,3) = true;    
end