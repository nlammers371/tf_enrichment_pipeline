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
paramIncVec = [5 5 5 5];
% keep_prediction_flag = false;
% designate simulation type
simTypeCell = {'koff_only_2','kon_only_2','out_only','in_only'};
tfDependentParamCell = {'koff','kon','ks','ka'};
for s = 1%:length(simTypeCell)
    simType = simTypeCell{s};
    nParamIncrement = paramIncVec(s);
    % call parameter sweep algorithm to conduct initial search of space
    sweepInfo = io_sweep_wrapper(resultsRoot,nParamIncrement,simType,[],false,'granularity',1e-1);
end

%% now identify best perfomers wrpt each metric

% RA
[~,best_ra_i] = nanmax(sweepInfo.ra_fit);
[~,worst_ra_i] = nanmin(sweepInfo.ra_fit);
param_vec_ra = sweepInfo.param_val_vec([best_ra_i,worst_ra_i],:);
sweepInfoRA = io_sweep_wrapper(resultsRoot,nParamIncrement,simType,param_vec_ra,true);

% MF
[~,best_ft_i] = nanmax(sweepInfo.fluo_time_fit);
[~,worst_ft_i] = nanmin(sweepInfo.fluo_time_fit);
param_vec_ft = sweepInfo.param_val_vec([best_ft_i,worst_ft_i],:);
sweepInfoFT = io_sweep_wrapper(resultsRoot,nParamIncrement,simType,param_vec_ft,true);

%%

% identify best thousand or best 1% (whichever is less) and run
% sweeps that save key info
sweepInfoBest = sweepInfo;
% save simulation type
sweepInfoBest.simType = simType;    
% load markov system parameter info
sweepInfoBest = getMarkovSystemInfo(sweepInfoBest);    
% generate ground truth reference curves       
sweepInfoBest = addGroundTruthFields(sweepInfoBest, io_ref_ra, io_ref_wt);

% generate array of weights to apply to different metrics. This will
% allow us to extract different versions of "optimal" networks
sweepInfoBest.fit_fields_to_use = {'ra_full_fit_R2','mean_fluo_fit_R2','off_time_fit_R2'};
n_fit_fields = length(sweepInfoBest.fit_fields_to_use);
sweepInfoBest.n_keep = min([size(sweepInfo.param_val_vec,1) max([100 min([round(0.01*size(sweepInfo.param_val_vec,1)), 1e3])])]);

sweepInfoBest.resultsRoot = resultsRoot;

% now iterate throug and find the best-scoring parameter set for each set
% of weights
sweepInfoBest = sweepBestPerformers(sweepInfoBest,sweepInfo);   

save([resultsRoot 'sweepInfoBest_' simType '.mat'],'sweepInfoBest', '-v7.3');