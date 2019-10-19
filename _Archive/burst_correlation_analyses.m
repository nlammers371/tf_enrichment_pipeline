% Script to attempt a systematic dissection of various factors driving
% proteinxtranscription burst coincidence
clear
close all
addpath('utilities')
% define core ID variables
project = 'Dl-Ven_snaBAC-mCh';
% project = 'Dl-Ven_hbP2P-mCh';
dropboxFolder =  'E:\Nick\LivemRNA\Dropbox\';
dataPath = [dropboxFolder 'ProcessedEnrichmentData\' project '\'];
figPath = [dropboxFolder 'LocalEnrichmentFigures\' project '\_paper_figures\input_output02\'];
mkdir(figPath)
% load data
load([dataPath 'hmm_input_output_results.mat'])
w = 7;
K = 3;
load([dataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'],'hmm_input_output');
nBoots = 100;
% define size of window of interest
roi_window = 6; % in time steps
window_size = 15;
start = window_size + 2;

% extract roi_vectors from swap and locus arrays
spot_protein_vec = nanmean(results_struct.spot_array_dt(:,start:start + roi_window-1),2) - ...
    nanmean(results_struct.spot_array_dt(:,start-roi_window-1:start-2),2);
swap_protein_vec = nanmean(results_struct.swap_array_dt(:,start:start + roi_window-1),2) - ...
    nanmean(results_struct.swap_array_dt(:,start-roi_window-1:start-2),2);
virt_protein_vec = nanmean(results_struct.virtual_array_dt(:,start:start + roi_window-1),2) - ...
    nanmean(results_struct.virtual_array_dt(:,start-roi_window-1:start-2),2);
mf_protein_vec = nanmean(results_struct.mf_array(:,start-roi_window-1:start + roi_window-1),2);

% pull other trend vectors
feature_sign_vec = results_struct.feature_sign_vec';
lag_size_vec = results_struct.lag_size_vec';
lead_size_vec = results_struct.lead_size_vec';
lag_dur_vec = results_struct.lag_dur_vec';
lead_dur_vec = results_struct.lead_dur_vec';
tr_burst_size_vec = lag_dur_vec.*lag_size_vec;

% min size of preceding lag
min_prev_lag = 5;
% make rise filter
rise_ft = feature_sign_vec == 1;
analysis_ft = rise_ft & lead_dur_vec> min_prev_lag & ~isnan(spot_protein_vec) & ~isnan(lag_dur_vec);

%%% First make simple bivariate plots indicating trends 

% (1) protein burst size vs. average Dl concentration
mf_protein_analysis = mf_protein_vec(analysis_ft);
lag_dur_analysis = lag_dur_vec(analysis_ft);
swap_analysis = swap_protein_vec(analysis_ft);
spot_analysis = spot_protein_vec(analysis_ft);
virt_analysis = virt_protein_vec(analysis_ft);
% calculate bins
mf_bins = linspace(prctile(mf_protein_analysis,5),prctile(mf_protein_analysis,95),25);
mf_sigma = 1*median(diff(mf_bins));
% obtain average protein enrichment 
spot_mf_array = NaN(nBoots,numel(mf_bins));
swap_mf_array = NaN(nBoots,numel(mf_bins));
virt_mf_array = NaN(nBoots,numel(mf_bins));
sample_index = 1:sum(analysis_ft);
for n = 1:nBoots
    s_ids = randsample(sample_index,sum(analysis_ft),true);    
    spot_boot = spot_analysis(s_ids);
    swap_boot = swap_analysis(s_ids);
    virt_boot = virt_analysis(s_ids);
    mf_boot = mf_protein_analysis(s_ids);
    for m = 1:numel(mf_bins)
        mf_val = mf_bins(m);
        mf_weights = exp(-.5*((mf_boot-mf_val)/mf_sigma).^2);
        spot_mf_array(n,m) = nansum(spot_boot.*mf_weights) / nansum(mf_weights);
        swap_mf_array(n,m) = nansum(swap_boot.*mf_weights) / nansum(mf_weights);
        virt_mf_array(n,m) = nansum(virt_boot.*mf_weights) / nansum(mf_weights);
    end
end

% calculate bootstrap mean and standard error
spot_mf_mean = nanmean(spot_mf_array);
spot_mf_ste = nanstd(spot_mf_array);

swap_mf_mean = nanmean(swap_mf_array);
swap_mf_ste = nanstd(swap_mf_array);

virt_mf_mean = nanmean(virt_mf_array);
virt_mf_ste = nanstd(virt_mf_array);

% make figure showing spot trend
mf_dl_fig = figure;
hold on
cmap1 = brewermap([],'Set2');
e = errorbar(mf_bins, spot_mf_mean, spot_mf_ste,'Color','black','LineWidth',1);
e.CapSize = 0;
s = scatter(mf_bins, spot_mf_mean,'MarkerFaceColor',cmap1(2,:)/1.2,...
    'MarkerEdgeColor','black','LineWidth',.5);
% format axes
grid on
ylabel('average Dl enrichment (au)')
xlabel('nuclear Dl concentration (au)')
set(gca,'Fontsize',14)
axis([mf_bins(1) mf_bins(end) 0 65])
saveas(mf_dl_fig,[figPath 'spot_enrichment_vs_meanDl.tif'])
saveas(mf_dl_fig,[figPath 'spot_enrichment_vs_meanDl.pdf'])

% show controls as well
mf_dl_all_fig = figure;
hold on
% spot
e = errorbar(mf_bins, spot_mf_mean, spot_mf_ste,'Color','black','LineWidth',1);
e.CapSize = 0;
s1 = scatter(mf_bins, spot_mf_mean,'MarkerFaceColor',cmap1(2,:)/1.2,...
    'MarkerEdgeColor','black','LineWidth',.5);
% swap
e = errorbar(mf_bins, virt_mf_mean, virt_mf_ste,'Color','black','LineWidth',1);
e.CapSize = 0;
% virtual
s2 = scatter(mf_bins, virt_mf_mean,'MarkerFaceColor',cmap1(5,:)/1.2,...
    'MarkerEdgeColor','black','LineWidth',.5);
% swap
e = errorbar(mf_bins, swap_mf_mean, swap_mf_ste,'Color','black','LineWidth',1);
e.CapSize = 0;
s3 = scatter(mf_bins, swap_mf_mean,'MarkerFaceColor',cmap1(3,:)/1.2,...
    'MarkerEdgeColor','black','LineWidth',.5);
% format axes
grid on
ylabel('average Dl enrichment (au)')
xlabel('nuclear Dl concentration (au)')
legend([s1 s2 s3], 'locus', 'virtual spot', 'nn control','Location','northwest')
set(gca,'Fontsize',14)
axis([mf_bins(1) mf_bins(end) -10 65])
saveas(mf_dl_all_fig,[figPath 'all_enrichment_vs_meanDl.tif'])
saveas(mf_dl_all_fig,[figPath 'all_enrichment_vs_meanDl.pdf'])


% (2) enrichment vs burst duration
% calculate bins
dur_bins = linspace(0,12,25);
dur_sigma = .5*median(diff(dur_bins));
% obtain average protein enrichment 
spot_dur_array = NaN(nBoots,numel(dur_bins));
swap_dur_array = NaN(nBoots,numel(dur_bins));
virt_dur_array = NaN(nBoots,numel(dur_bins));
sample_index = 1:sum(analysis_ft);
for n = 1:nBoots
    s_ids = randsample(sample_index,sum(analysis_ft),true);    
    spot_boot = spot_analysis(s_ids);
    swap_boot = swap_analysis(s_ids);
    virt_boot = virt_analysis(s_ids);
    dur_boot = lag_dur_analysis(s_ids);
    for d = 1:numel(dur_bins)
        dur_val = dur_bins(d);
        dur_weights = exp(-.5*((dur_boot-dur_val)/dur_sigma).^2);
        spot_dur_array(n,d) = nansum(spot_boot.*dur_weights) / nansum(dur_weights);
        swap_dur_array(n,d) = nansum(swap_boot.*dur_weights) / nansum(dur_weights);
        virt_dur_array(n,d) = nansum(virt_boot.*dur_weights) / nansum(dur_weights);
    end
end

% calculate bootstrap mean and standard error
spot_dur_mean = nanmean(spot_dur_array);
spot_dur_ste = nanstd(spot_dur_array);

swap_dur_mean = nanmean(swap_dur_array);
swap_dur_ste = nanstd(swap_dur_array);

virt_dur_mean = nanmean(virt_dur_array);
virt_dur_ste = nanstd(virt_dur_array);

%%
% spto trend only
dur_axis = dur_bins/60*20;
dur_dl_fig = figure;
hold on
e = errorbar(dur_axis, spot_dur_mean, spot_dur_ste,'Color','black','LineWidth',1);
e.CapSize = 0;
s = scatter(dur_axis, spot_dur_mean,'MarkerFaceColor',cmap1(2,:)/1.2,...
    'MarkerEdgeColor','black','LineWidth',.5);
% format axes
grid on
ylabel('average Dl enrichment (au)')
xlabel('{\itsnail} burst duration (au)')
set(gca,'Fontsize',14)
axis([dur_bins(1)/60*20 dur_bins(end)/60*20 -10 55])
saveas(dur_dl_fig,[figPath 'spot_enrichment_vs_durTr.tif'])
saveas(dur_dl_fig,[figPath 'spot_enrichment_vs_durTr.pdf'])

% show controls as well
dur_dl_all_fig = figure;
hold on
% spot
e = errorbar(dur_axis, spot_dur_mean, spot_dur_ste,'Color','black','LineWidth',1);
e.CapSize = 0;
s1 = scatter(dur_axis, spot_dur_mean,'MarkerFaceColor',cmap1(2,:)/1.2,...
    'MarkerEdgeColor','black','LineWidth',.5);
% swap
e = errorbar(dur_axis, virt_dur_mean, virt_mf_ste,'Color','black','LineWidth',1);
e.CapSize = 0;
% virtual
s2 = scatter(dur_axis, virt_dur_mean,'MarkerFaceColor',cmap1(5,:)/1.2,...
    'MarkerEdgeColor','black','LineWidth',.5);
% swap
e = errorbar(dur_axis, swap_dur_mean, swap_dur_ste,'Color','black','LineWidth',1);
e.CapSize = 0;
s3 = scatter(dur_axis, swap_dur_mean,'MarkerFaceColor',cmap1(3,:)/1.2,...
    'MarkerEdgeColor','black','LineWidth',.5);
% format axes
grid on
ylabel('average Dl enrichment (au)')
xlabel('{\itsnail} burst duration')
legend([s1 s2 s3], 'locus', 'virtual spot', 'nn control','Location','northwest')
set(gca,'Fontsize',14)
axis([dur_axis(1) dur_axis(end) -10 55])
saveas(dur_dl_all_fig,[figPath 'all_enrichment_vs_durTr.tif'])
saveas(dur_dl_all_fig,[figPath 'all_enrichment_vs_durTr.pdf'])

%%% linear regression to see if burst duration is predictive of enrichment 
burst_dur_norm = (lag_dur_analysis);
mf_norm = (mf_protein_analysis);
% locus
lin_spot = fitlm([mf_norm burst_dur_norm], spot_analysis);
tri_spot = fitlm([mf_norm mf_norm.^2 mf_norm.^3 burst_dur_norm burst_dur_norm.^2 ...
    burst_dur_norm.^3], spot_analysis);
% controls
lin_swap = fitlm([mf_norm burst_dur_norm], swap_analysis);
lin_virt = fitlm([mf_norm burst_dur_norm], virt_analysis);

% generate predicted trend and compare
pd_enrichment_trend = dur_bins*lin_spot.Coefficients.Estimate(3) + ...
    spot_dur_mean(1);

pd_dur_dl_fig = figure;
hold on
s = scatter(dur_axis, spot_dur_mean,'MarkerFaceColor',cmap1(2,:)/1.2,...
    'MarkerEdgeColor','black','LineWidth',.5);
p = plot(dur_axis, pd_enrichment_trend,'Color','black','LineWidth',1.5);
% format axes
grid on
ylabel('average Dl enrichment (au)')
xlabel('{\itsnail} burst duration (au)')
set(gca,'Fontsize',14)
legend([s p], 'enrichment at locus', 'prediction (linear trend)','Location','northwest')
axis([dur_bins(1)/60*20 dur_bins(end)/60*20 0 45])
saveas(pd_dur_dl_fig,[figPath 'spot_enrichment_vs_durTr_trend.tif'])
saveas(pd_dur_dl_fig,[figPath 'spot_enrichment_vs_durTr_trend.pdf'])
