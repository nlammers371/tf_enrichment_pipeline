% Script to scan through parameter space for selected biological variables
clear
close all

addpath(genpath('../utilities'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%% Load experimental data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

resultsRoot = 'S:\Nick\Dropbox\ProcessedEnrichmentData\parameterSweeps\';
if ~isfolder(resultsRoot)
  resultsRoot = 'C:\Users\nlamm\Dropbox (Personal)\ProcessedEnrichmentData\parameterSweeps\';
end

% load experimental reference data
load([resultsRoot 'io_ref_ra.mat'])
load([resultsRoot 'io_ref_wt.mat'])


% designate simulation type
simTypeCell = {'koff_only_2','kon_only_2','out_only','in_only'};
tfDependentParamCell = {'koff','kon','ks','ka'};

% initialize structure to store results
master_struct = struct;

% load sweep results
for s = 1:length(simTypeCell)
  
    simType = simTypeCell{s};
    
    % load full sweep results
    load([resultsRoot 'sweepInfo_' simType '.mat'],'sweepInfo');    
    master_struct(s).sweepInfo = sweepInfo;
    
    % load best-performing parameter sets   
    load([resultsRoot 'sweepInfoBest_' simType '.mat'],'sweepInfoBest');
    master_struct(s).sweepInfoBest = sweepInfoBest;
    
end


%% Retroactively calculate combined fit scores

