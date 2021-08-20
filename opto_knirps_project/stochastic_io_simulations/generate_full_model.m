% Function wrapper for stochastic trace simulations in non-steady-state
% conditions

function sweepInfo = generate_full_model(sweepInfo,varargin)


%% %%%%%%%%%%%%%%%%%%%%%% Check for optional inputs %%%%%%%%%%%%%%%%%%%%%%%
for i = 1:numel(varargin)
   if ischar(varargin{i}) && i <  numel(varargin)
       eval(['simInfo.' varargin{i} ' = varargin{i+1};'])
   end
end

if contains(sweepInfo.simType,'in') % this signifies 3 state system where rate from SILENT -> OFF is regulated
    % Note that, in this scenario, we need to back-calculate the "true" 2
    % state parameters
    R2 = sweepInfo.R2;
    kon_true = R2(2,1)*(sweepInfo.ka + sweepInfo.ks)/sweepInfo.ka;
    R2(2,1) = kon_true;
    R2(1,1) = -kon_true;
    sweepInfo.R2 = R2;
elseif strcmp(sweepInfo.simType,'koff_only_2')
    sweepInfo.R2(1,2) = sweepInfo.rate_max;
end    

if contains(sweepInfo.simType,'on')
    sweepInfo.frac_init = 1e-3;
    sweepInfo.frac_final = 1-1e-3;
end

% Add third, silent state
sweepInfo.RateMatrix(2:3,2:3) = sweepInfo.R2;
if ~contains(sweepInfo.simType,'2') && ~strcmp(sweepInfo.simType,'match_exp') 
    sweepInfo.RateMatrix(:,1) = [0 ; sweepInfo.ka ; 0];
    sweepInfo.RateMatrix(1,:) = [0 sweepInfo.ks 0];
end    
diag_flags = eye(size(sweepInfo.RateMatrix,1))==1;
sweepInfo.RateMatrix(diag_flags) = 0;
sweepInfo.RateMatrix(diag_flags) = -sum(sweepInfo.RateMatrix);


sweepInfo.r_emission = [0 sweepInfo.r2]; % loading rate for each state
sweepInfo.noise = sweepInfo.noise;
sweepInfo.pi0 = [0 sweepInfo.pi0];
if contains(sweepInfo.simType,'_on') && ~contains(sweepInfo.simType,'2')
    sweepInfo.pi0(1) = 1;
end
% simInfo.noise = 1e4; 

% dictate which rate is tf-dependent (assume only one possible for now)
sweepInfo.tf_dependent_flags = false(size(sweepInfo.RateMatrix));
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
