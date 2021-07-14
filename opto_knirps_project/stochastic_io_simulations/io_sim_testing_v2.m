% Script to scan through parameter space for selected biological variables
clear
close all

addpath(genpath('../utilities'));

% %%%%%%%%%%%%%%%%%%%%%% Specify initial settings %%%%%%%%%%%%%%%%%%%%%%%

% initialize data structure to store simulation info
simInfo = struct;
simInfo.simType = 'out_only_off';

% define basic sim parameters
simInfo.seq_length = 120;
simInfo.memory = 7;
simInfo.deltaT = 20;
simInfo.t_MS2 = 1.4;
simInfo.n_traces = 100;
simInfo.granularity = 1;

% define response parameters
simInfo.KD = 10;
simInfo.HC = 6;

% specify network architecture
simInfo.RateMatrix = [-1 60 0; 
                      1 -62 1; 
                      0  2 -1]/60;
simInfo.r_emission = [0 0 1]; % loading rate for each state
simInfo.pi0 = [0 0.5 0.5];
simInfo.noise = 1; 

simInfo.frac_init = 1-1e-3;
simInfo.frac_final = 1e-3;

%% %%%%%%%%%%%%%%%%%%%%%% Check for optional inputs %%%%%%%%%%%%%%%%%%%%%%%
% for i = 1:numel(varargin)
%    if ischar(varargin{i}) && i <  numel(varargin)
%        eval(['simInfo.' varargin{i} ' = varargin{i+1};'])
%    end
% end

% dictate which rate is tf-dependent (assume only one possible for now)
simInfo.tf_dependent_flags = false(size(simInfo.RateMatrix));
if strcmp(simInfo.simType,'out_only_off')      
    simInfo.tf_dependent_flags(1,2) = true;
elseif strcmp(simInfo.simType,'in_only_off')      
    simInfo.tf_dependent_flags(2,1) = true;
end

% simulate TF profiles
simInfo.tf_profile_array = simulate_tf_profiles(simInfo.frac_init,simInfo.frac_final,simInfo);

%%

tic
gillespie = synthetic_rate_gillespie_io_v2(simInfo);
toc

% calculate last shut-off times
detection_thresh = 10;
on_off_array = ones(size(gillespie.fluo_ms2_array));
for n = 1:n_traces
    last_i = find(gillespie.fluo_ms2_array(:,n)>detection_thresh,1,'last');
    on_off_array(last_i+1:end,n) = 0;
end