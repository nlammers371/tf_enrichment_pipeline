% Script to scan through parameter space for selected biological variables
clear
close all

addpath(genpath('../utilities'));

% designate simulation type
simType = 'out_only_on';%'out_only_off';

% specify 2 state network architecture (eventually this will be drawn from
% actual fits)
systemParams = struct;
systemParams.R2 = [-2  1; 
                    2 -1]/60;
systemParams.r2 = [0 1]*1e4; % loading rate for each state
systemParams.pi0 = [0.5 0.5];
if contains(simType,'on')
    systemParams.pi0 = [0 0];
end
systemParams.noise = 5e3; 

% designate key variables to sweep
KD = 3e4;
HC = 5;
K_out = 1/60;%linspace(1,60,10);
if contains(simType,'out')
    K_out = 60;
end
K_in = 1/60;
F_min = 1e4;

%% call stochastic simulation function
tic
simInfo = io_sim_function(simType,systemParams,KD,HC,K_out(end),K_in,'n_traces',10,'granularity',1);
toc

%% use output to generate predicted cumulative OFF (or ON) curve(s)
simInfo.F_min = logspace(3,log(5e4));
tic
simInfo = calculate_cumulative_dist(simInfo);
toc