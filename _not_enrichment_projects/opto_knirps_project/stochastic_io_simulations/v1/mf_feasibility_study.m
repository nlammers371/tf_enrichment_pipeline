% script to make basic performance figures
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
% mak figure directory  
date_string = '22-Sep-2021';
resultsPath = [resultsRoot date_string filesep];

try
    FigurePath = [liveProject.figurePath 'io_sim_results' filesep date_string filesep];
catch
    FigurePath = ['C:\Users\nlamm\Dropbox (Personal)\LocalEnrichmentFigures\PipelineOutput\paramSweeps\' filesep date_string filesep];
end
mkdir(FigurePath)

paramIncVec = [25 25 16 16];
simTypeCell = {'koff_only_2','kon_only_2','out_only','in_only'};
tfDependentParamCell = {'koff','kon','ks','ka'};
master_struct = struct;

% load sweep results
for s = 1:length(simTypeCell)
    simType = simTypeCell{s};
    nInc = paramIncVec(s);
    readPath = [resultsPath 'sweeps_n' num2str(nInc) filesep];
    % load sweep result
    load([readPath 'sweepInfo_' simType '.mat'],'sweepInfo');
    
    % store
    master_struct(s).sweepInfo = sweepInfo;
end

% load raw data
load([resultsRoot 'io_ref_ra.mat'],'io_ref_ra')
load([resultsRoot 'io_ref_wt.mat'],'io_ref_wt')
load([resultsRoot 'io_ref_ON.mat'],'io_ref_ON')

%% obtain predicted trends for best systems from sim type
n_traces = 250;
n_best = 10;

for s = 1:length(master_struct)
    % extract
    sweepInfo = master_struct(s).sweepInfo;
    simType = sweepInfo.simType;
    
    % calculate aggregate score
    total_score = -sqrt(sweepInfo.fluo_time_ON_fit.^2 + sweepInfo.fluo_time_fit.^2 + sweepInfo.ra_fit.^2 + sweepInfo.still_on_fit.^2);

    % find best overall performers
    [~,score_ids] = sort(total_score,'descend');
    best_i_list = score_ids(1:n_best);
    
    % now find best parameter-specific scores
    [~,best_ft_i] = nanmax(sweepInfo.fluo_time_fit);
    [~,best_ft_ON_i] = nanmax(sweepInfo.fluo_time_ON_fit);
    [~,best_pon_i] = nanmax(sweepInfo.still_on_fit);
    [~,best_ra_i] = nanmax(sweepInfo.ra_fit);
    
    % get corresponding parameters
    param_vec = sweepInfo.param_val_vec([best_i_list',best_ft_i,best_ra_i,best_pon_i,best_ft_ON_i],:);
    
    % run sweep for selected networks
    sweepInfoBest = io_sweep_wrapper(resultsRoot,2,sweepInfo.simType,param_vec,true,'n_traces',n_traces);
    
    % store
    master_struct(s).sweepInfoBest = sweepInfoBest;
    master_struct(s).best_i_list = best_i_list;
    master_struct(s).best_ft_i = best_ft_i;
    master_struct(s).best_ra_i = best_ra_i;
    master_struct(s).best_pon_i = best_pon_i;
end

%% 
%%%%%%%%%%%%%%%%%%
% koff
sim_i = 1;
sub_ind = 1;
sweepInfo = master_struct(sim_i).sweepInfoBest;
paramVals = sweepInfo.param_val_vec(sub_ind,:);
ms2_trace_array = sweepInfo.ms2_traces_true_ON(:,:,sub_ind);
tf_rate_true = sweepInfo.tf_dependent_curves_ON(:,:,sub_ind);
% calculate predicted average rate
rate_pd_array = sweepInfo.r2(2) * paramVals(end-1) ./ (tf_rate_true + paramVals(end-1));
% calculate predicted fluorescence
fluo_kernel = ms2_loading_coeff (sweepInfo.t_MS2, sweepInfo.memory);
fluo_pd_array = convn(fluo_kernel',rate_pd_array,'full');

koff_fig = figure;
hold on
plot(nanmean(fluo_pd_array,2))
plot(nanmean(ms2_trace_array,2))
ylabel('average MS2 signal')
xlabel('time steps')
legend('mf prediction (koff)','simulation (koff)')


%%%%%%%%%%%%%%%%%%
% kon
sim_i = 2;
sub_ind = 1;
sweepInfo = master_struct(sim_i).sweepInfoBest;
paramVals = sweepInfo.param_val_vec(sub_ind,:);
ms2_trace_array = sweepInfo.ms2_traces_true_wt(:,:,sub_ind);
ms2_trace_array_nan = sweepInfo.ms2_traces_observed_wt(:,:,sub_ind);
ms2_trace_array_nan(isnan(ms2_trace_array_nan)) = 0;
tf_rate_true = sweepInfo.tf_dependent_curves_wt(:,:,sub_ind);
% calculate predicted average rate
rate_pd_array = sweepInfo.r2(2) * tf_rate_true ./ (tf_rate_true + paramVals(end));
% calculate predicted fluorescence
fluo_pd_array = convn(fluo_kernel',rate_pd_array,'full');
fluo_pd_array = fluo_pd_array(1:end-sweepInfo.memory+1,:);
knirps_traces = sweepInfo.knirps_traces_ON(:,:,sub_ind);

% exract experimental profile
exp_raw = io_ref_wt.fluo_array;
ap_bounds = io_ref_wt.ap_limits_still_on;
mean_ap = io_ref_wt.mean_ap;
ap_filter = mean_ap>=ap_bounds(1)&mean_ap<=ap_bounds(2);
exp_ap = exp_raw(:,ap_filter);
exp_ap(isnan(exp_ap)) = 0;

kon_fig = figure;
hold on
plot(sweepInfo.time_axis_wt,nanmean(fluo_pd_array,2))
plot(sweepInfo.time_axis_wt,nanmean(ms2_trace_array,2))
plot(sweepInfo.time_axis_wt,nanmean(ms2_trace_array_nan,2))
plot(io_ref_wt.time_axis,nanmean(exp_ap,2))

ylabel('average MS2 signal')
xlabel('time steps')
legend('mf prediction (kon)','simulation (kon)','simulation with threshold (kon)')

%%%%%%%%%%%%%%%%%%
%% kon
sim_i = 2;
sub_ind = 1;
sweepInfo = master_struct(sim_i).sweepInfoBest;
paramVals = sweepInfo.param_val_vec(sub_ind,:);
ms2_trace_array = sweepInfo.ms2_traces_true_ra(:,:,sub_ind);
ms2_trace_array_nan = sweepInfo.ms2_traces_observed_ra(:,:,sub_ind);
ms2_trace_array_nan(isnan(ms2_trace_array_nan)) = 0;
tf_rate_true = sweepInfo.tf_dependent_curves_ra(:,:,sub_ind);
% calculate predicted average rate
rate_pd_array = sweepInfo.r2(2) * tf_rate_true ./ (tf_rate_true + paramVals(end));
% calculate predicted fluorescence
fluo_pd_array = convn(fluo_kernel',rate_pd_array,'full');
fluo_pd_array = fluo_pd_array(1:end-sweepInfo.memory+1,:);
knirps_traces = sweepInfo.knirps_traces_ra(:,:,sub_ind);

% exract experimental profile
exp_ap = io_ref_ra.fluo_array;
% ap_bounds = io_ref_ON.ap_limits_still_on;
% mean_ap = io_ref_ON.mean_ap;
% ap_filter = mean_ap>=ap_bounds(1)&mean_ap<=ap_bounds(2);
% exp_ap = exp_raw(:,ap_filter);
exp_ap(isnan(exp_ap)) = 0;

kon_fig = figure;
hold on
plot(sweepInfo.time_axis_ra,nanmean(fluo_pd_array,2))
plot(sweepInfo.time_axis_ra,nanmean(ms2_trace_array,2))
plot(sweepInfo.time_axis_ra,nanmean(ms2_trace_array_nan,2))
plot(io_ref_ra.time_vec,nanmean(exp_ap,2))

ylabel('average MS2 signal')
xlabel('time steps')
legend('mf prediction (kon)','simulation (kon)','simulation with threshold (kon)')
%% make kon input-output plot
kni_vec = linspace(0, 10);

kon_fig = figure;
hold on
plot(kni_vec,60*paramVals(end-1)*paramVals(2)^paramVals(1)./(paramVals(2)^paramVals(1) + kni_vec.^paramVals(1)),'LineWidth',1.5)
plot(kni_vec,repelem(60*sweepInfo.R2_orig(2,1),length(kni_vec)),'--k','LineWidth',1.5)
xlabel('Knirps concentration (au)')
ylabel('k_{on} (events per minute)')


