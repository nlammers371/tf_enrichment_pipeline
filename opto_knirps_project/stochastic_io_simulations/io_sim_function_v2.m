% Function wrapper for stochastic trace simulations in non-steady-state
% conditions

function simInfo = io_sim_function_v2(simType,systemParams,KD,HC,K_out,K_in,tf_profile_array,n_traces,granularity,varargin)

% close all
addpath(genpath('../../utilities'));

% %%%%%%%%%%%%%%%%%%%%%% Specify initial settings %%%%%%%%%%%%%%%%%%%%%%%

% initialize data structure to store simulation info
simInfo = struct;
simInfo.simType = simType;%'out_only_off';

% define basic sim parameters
simInfo.seq_length = systemParams.seq_length;
simInfo.memory = systemParams.memory;
simInfo.deltaT = systemParams.deltaT;
simInfo.t_MS2 = systemParams.t_MS2;

simInfo.n_traces = n_traces;
simInfo.granularity = granularity;

% define response parameters
simInfo.KD = KD;
simInfo.HC = HC;
simInfo.K_out = K_out;
simInfo.K_in = K_in;

% specify network architecture
simInfo.systemParams = systemParams;

simInfo.frac_init = 1-1e-3;
simInfo.frac_final = 1e-3;

%% %%%%%%%%%%%%%%%%%%%%%% Check for optional inputs %%%%%%%%%%%%%%%%%%%%%%%
for i = 1:numel(varargin)
   if ischar(varargin{i}) && i <  numel(varargin)
       eval(['simInfo.' varargin{i} ' = varargin{i+1};'])
   end
end

if contains(simInfo.simType,'in')
    % Note that, in this scenario, we need to back-calculate the "true" 2
    % state parameters
    R2 = systemParams.R2;
    kon_true = R2(2,1)*(K_in + K_out)/K_in;
    R2(2,1) = kon_true;
    R2(1,1) = -kon_true;
    systemParams.R2 = R2;
end

if contains(simInfo.simType,'on')
    simInfo.frac_init = 1e-3;
    simInfo.frac_final = 1-1e-3;
end

% Add third, silent state
simInfo.RateMatrix(2:3,2:3) = systemParams.R2;
simInfo.RateMatrix(:,1) = [0 ; K_in ; 0];
simInfo.RateMatrix(1,:) = [0 K_out 0];
diag_flags = eye(size(simInfo.RateMatrix,1))==1;
simInfo.RateMatrix(diag_flags) = 0;
simInfo.RateMatrix(diag_flags) = -sum(simInfo.RateMatrix);


simInfo.r_emission = [0 systemParams.r2]; % loading rate for each state
simInfo.noise = systemParams.noise;
simInfo.pi0 = [0 systemParams.pi0];
if contains(simInfo.simType,'on')
    simInfo.pi0(1) = 1;
end
% simInfo.noise = 1e4; 

% dictate which rate is tf-dependent (assume only one possible for now)
simInfo.tf_dependent_flags = false(size(simInfo.RateMatrix));
if contains(simInfo.simType,'out')      
    simInfo.tf_dependent_flags(1,2) = true;
elseif contains(simInfo.simType,'in')      
    simInfo.HC = -simInfo.HC;
    simInfo.tf_dependent_flags(2,1) = true;
end

% simulate TF profiles
if isempty(tf_profile_array)
    simInfo.tf_profile_array = simulate_tf_profiles(simInfo.frac_init,simInfo.frac_final,simInfo);
else
    s_indices = randsample(1:size(tf_profile_array,3),simInfo.n_traces,true);
    simInfo.tf_profile_array = tf_profile_array(:,1,s_indices);
end
%% call simulation function
simInfo.gillespie = synthetic_rate_gillespie_io_v2(simInfo);

