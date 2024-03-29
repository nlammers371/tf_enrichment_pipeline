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
paramIncVec = [25 25 16 16];
% keep_prediction_flag = false;
% designate simulation type
simTypeCell = {'koff_only_2','kon_only_2','out_only','in_only'};
tfDependentParamCell = {'koff','kon','ks','ka'};
master_struct = struct;

for s = 1:length(simTypeCell)
    simType = simTypeCell{s};
    nParamIncrement = paramIncVec(s);
    % call parameter sweep algorithm to conduct initial search of space
    sweepInfo = io_sweep_wrapper(resultsRoot,nParamIncrement,simType,[],false,'granularity',1);
    master_struct(s).sweepInfo = sweepInfo;
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
s = 4;
sweepInfo = master_struct(s).sweepInfo;
total_score = -sqrt(sweepInfo.fluo_time_fit.^2 + sweepInfo.fluo_time_ON_fit.^2 + sweepInfo.ra_fit.^2 + sweepInfo.still_on_fit.^2);
% ra_filter = sweepInfo.ra_fit>-0.5;
% total_score(~ra_filter) = -Inf;
[~,best_combined] = nanmax(total_score);
[~,best_ft_ON_i] = nanmax(sweepInfo.fluo_time_ON_fit);
[~,best_ft_i] = nanmax(sweepInfo.fluo_time_fit);
[~,best_pon_i] = nanmax(sweepInfo.still_on_fit);
[~,best_ra_i] = nanmax(sweepInfo.ra_fit);
param_vec = sweepInfo.param_val_vec([best_combined,best_ft_i,best_ft_ON_i,best_ra_i,best_pon_i],:);
sweepInfoBest = io_sweep_wrapper(resultsRoot,nParamIncrement,sweepInfo.simType,param_vec,true,'n_traces',1e3);
%%
close all

figure;
plot(sweepInfoBest.ra_time_cdf_predicted')
hold on
plot(sweepInfoBest.reactivation_cdf,'-k')

figure;
plot(permute(sweepInfoBest.fluo_time_predicted,[1 3 2]))
hold on
plot(sweepInfoBest.fluo_time_mean,'-k')

figure;
plot(permute(sweepInfoBest.fluo_time_predicted_ON,[1 3 2]))
hold on
plot(sweepInfoBest.fluo_time_mean_ON,'-k')

figure;
plot(permute(sweepInfoBest.p_still_on_predicted,[2 3 1]))
hold on
plot(sweepInfoBest.fraction_still_on,'-k')

%% identify ranges of parameter space that coincide with optimal performance
% Let's sweep these regions in more detail
s = 3;
sweepInfo = master_struct(s).sweepInfo;
fit_flags = sweepInfo.fitFlags;

% let's take best 5% of each metric and overall
prctile_cutoff = 99.99;
% overall
total_score = -sqrt(sweepInfo.fluo_time_fit.^2 + sweepInfo.ra_fit.^2 + sweepInfo.still_on_fit.^2);
prctile_total = prctile(total_score,prctile_cutoff);
best_total_ids = find(total_score>=prctile_total);
% ra
prctile_ra = prctile(sweepInfo.ra_fit,prctile_cutoff);
best_ra_ids = find(sweepInfo.ra_fit>=prctile_ra);
% ra
prctile_ft = prctile(sweepInfo.fluo_time_fit,prctile_cutoff);
best_ft_ids = find(sweepInfo.fluo_time_fit>=prctile_ft);
% ra
prctile_pon = prctile(sweepInfo.still_on_fit,prctile_cutoff);
best_pon_ids = find(sweepInfo.still_on_fit>=prctile_pon);

% get unique list
best_ids = unique([best_total_ids]);

% calculate ranges
best_params = sweepInfo.param_val_vec(best_ids,:);
best_param_max = nanmax(best_params)
best_param_min = nanmin(best_params)


