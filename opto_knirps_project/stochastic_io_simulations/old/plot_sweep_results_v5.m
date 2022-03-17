test(% Script to scan through parameter space for selected biological variables
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
try
    FigurePath = [liveProject.figurePath 'io_sim_results' filesep];
catch
    FigurePath = 'C:\Users\nlamm\Dropbox (Personal)\LocalEnrichmentFigures\PipelineOutput\paramSweeps\';
end
mkdir(FigurePath)

% load experimental reference data
load([resultsRoot 'io_ref_ra.mat'])
load([resultsRoot 'io_ref_wt.mat'])


% designate simulation type
simTypeCell = {'koff_only_2','kon_only_2','out_only','in_only'};
tfDependentParamCell = {'koff','kon','ks','ka'};
labelCell = {'k_{off}','k_{on}','k_{s}','k_{a}'};
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

%% Plot N best-performing networks for each model type
close all

n_plot = 25;
labelCell2 = [{'data'} labelCell];
ap_lims = [master_struct(s).sweepInfoBest.ap_axis_mean(1)-0.005 master_struct(s).sweepInfoBest.ap_axis_mean(end)+0.005];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Mean fluorescence
for s = 1:length(simTypeCell)
  
    mean_fluo_fig = figure;
    cmap = brewermap([],'Set2');
    hold on
    simType = simTypeCell{s};
    
    % plot experimental data trend
    errorbar(master_struct(1).sweepInfoBest.ap_axis_mean,master_struct(1).sweepInfoBest.mean_fluo_ap,...
                                master_struct(1).sweepInfoBest.mean_fluo_ap_ste,'Color','k','Capsize',0);
                              
    
    plots(1) = scatter(master_struct(1).sweepInfoBest.ap_axis_mean,master_struct(1).sweepInfoBest.mean_fluo_ap,75,'MarkerFaceColor',cmap(8,:),'MarkerEdgeColor','k');

    % calculate approximate feasible range for 
    plot(master_struct(s).sweepInfoBest.ap_axis_mean, master_struct(s).sweepInfoBest.mean_fluo_predicted(1:n_plot,:)','Color',[cmap(s,:) 0.75])
%     plots(end+1) = plot(master_struct(s).sweepInfoBest.ap_axis_mean, master_struct(s).sweepInfoBest.mean_fluo_predicted(1,:)','Color',[cmap(s,:) 1],'LineWidth',2);  
    
    xlabel('relative AP position');
    ylabel('average spot fluorescence (au)')
%     legend(plots,labelCell2{:},'Orientation','horizontal','Location','south');

    set(gca,'FontSize',14)

    set(gca,'Color',[228,221,209]/255) 
    grid on
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.XAxis(1).Color = 'k';
    xlim(ap_lims)
    ylim([0 1.8e5])
    mean_fluo_fig.InvertHardcopy = 'off';
    set(gcf,'color','w');
    title(labelCell{s});
    
    saveas(mean_fluo_fig,[FigurePath simType '_mean_fluo_fits.png'])
end  
                          

%%
mean_fluo_fig = figure;
cmap = brewermap([],'Set2');
hold on

% plot experimental data trend
errorbar(master_struct(1).sweepInfoBest.ap_axis_mean,master_struct(1).sweepInfoBest.mean_fluo_ap,...
                            master_struct(1).sweepInfoBest.mean_fluo_ap_ste,'Color','k','Capsize',0);
plots = [];
plots(1) = scatter(master_struct(1).sweepInfoBest.ap_axis_mean,master_struct(1).sweepInfoBest.mean_fluo_ap,75,'MarkerFaceColor',cmap(8,:),'MarkerEdgeColor','k');

for s = 1:length(simTypeCell)
%     plot(master_struct(s).sweepInfoBest.ap_axis_mean, master_struct(s).sweepInfoBest.mean_fluo_predicted(1:n_plot,:)','Color',[cmap(s,:) 0.25])
    plots(end+1) = plot(master_struct(s).sweepInfoBest.ap_axis_mean, master_struct(s).sweepInfoBest.mean_fluo_predicted(1,:)','Color',[cmap(s,:) 1],'LineWidth',2);  
end  
                          

xlabel('relative AP position');
ylabel('average spot fluorescence (au)')
legend(plots,labelCell2{:},'Orientation','horizontal','Location','south');

set(gca,'FontSize',14)

set(gca,'Color',[228,221,209]/255) 
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
xlim(ap_lims)
ylim([0 1.8e5])
mean_fluo_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(mean_fluo_fig,[FigurePath 'mean_fluo_fits.png'])

%% Average off time

for s = 1:length(simTypeCell)
  
    off_time_fig = figure;
    cmap = brewermap([],'Set2');
    hold on
    simType = simTypeCell{s};
    
    % plot experimental data trend
    errorbar(master_struct(1).sweepInfoBest.ap_axis_mean,master_struct(1).sweepInfoBest.off_time_ap/60,...
                                master_struct(1).sweepInfoBest.off_time_ap_ste/60,'Color','k','Capsize',0);
                              
    
    plots(1) = scatter(master_struct(1).sweepInfoBest.ap_axis_mean,master_struct(1).sweepInfoBest.off_time_ap/60,75,'MarkerFaceColor',cmap(8,:),'MarkerEdgeColor','k');

    % calculate approximate feasible range for 
    plot(master_struct(s).sweepInfoBest.ap_axis_mean, (master_struct(s).sweepInfoBest.off_time_predicted(1:n_plot,:)/60)','Color',[cmap(s,:) 0.75])
%     plots(end+1) = plot(master_struct(s).sweepInfoBest.ap_axis_mean, master_struct(s).sweepInfoBest.mean_fluo_predicted(1,:)','Color',[cmap(s,:) 1],'LineWidth',2);  
    
    xlabel('relative AP position');
    ylabel('average off time (minutes)')
%     legend(plots,labelCell2{:},'Orientation','horizontal','Location','south');

    set(gca,'FontSize',14)
    ylim([10 40])
    set(gca,'Color',[228,221,209]/255) 
    grid on
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.XAxis(1).Color = 'k';
    xlim(ap_lims)
%     ylim([0 1.8e5])
    off_time_fig.InvertHardcopy = 'off';
    set(gcf,'color','w');
    title(labelCell{s});
    
    saveas(off_time_fig,[FigurePath simType '_off_time_fits.png'])
end  

%%
off_time_fig = figure;
cmap = brewermap([],'Set2');
hold on

% plot experimental data trend
errorbar(master_struct(1).sweepInfoBest.ap_axis_mean,master_struct(1).sweepInfoBest.off_time_ap/60,...
                            master_struct(1).sweepInfoBest.off_time_ap_ste/60,'Color','k','Capsize',0);
plots = [];
plots(1) = scatter(master_struct(1).sweepInfoBest.ap_axis_mean,master_struct(1).sweepInfoBest.off_time_ap/60,75,'MarkerFaceColor',cmap(8,:),'MarkerEdgeColor','k');

for s = 1:length(simTypeCell)
%     plot(master_struct(s).sweepInfoBest.ap_axis_mean, master_struct(s).sweepInfoBest.off_time_predicted'/60,'Color',[cmap(s,:) 0.5])
    plots(end+1) = plot(master_struct(s).sweepInfoBest.ap_axis_mean, master_struct(s).sweepInfoBest.off_time_predicted(1,:)/60,'Color',[cmap(s,:) 1],'LineWidth',2);  
end  
                          

xlabel('relative AP position');
ylabel('average off time (minutes)')
legend(plots,labelCell2{:},'Orientation','horizontal','Location','south');

set(gca,'FontSize',14)

set(gca,'Color',[228,221,209]/255) 
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
xlim(ap_lims)
ylim([10 40])
off_time_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(off_time_fig,[FigurePath 'off_time_fits.png'])

%%

for s = 1:length(simTypeCell)
  
    ra_fig = figure;
    cmap = brewermap([],'Set2');
    hold on
    simType = simTypeCell{s};
    
    % plot experimental data trend
    errorbar(master_struct(1).sweepInfoBest.reactivation_time/60,master_struct(1).sweepInfoBest.reactivation_cdf,...
                                master_struct(1).sweepInfoBest.reactivation_cdf_ste,'Color','k','Capsize',0);
                                  
    plots(1) = scatter(master_struct(1).sweepInfoBest.reactivation_time/60,master_struct(1).sweepInfoBest.reactivation_cdf,75,'MarkerFaceColor',cmap(8,:),'MarkerEdgeColor','k');

    % calculate approximate feasible range for 
    plot(master_struct(1).sweepInfoBest.reactivation_time/60, (master_struct(s).sweepInfoBest.ra_time_cdf_predicted(1:n_plot,:))','Color',[cmap(s,:) 0.75])
%     plots(end+1) = plot(master_struct(s).sweepInfoBest.ap_axis_mean, master_struct(s).sweepInfoBest.mean_fluo_predicted(1,:)','Color',[cmap(s,:) 1],'LineWidth',2);  
    
    xlabel('relative AP position');
    ylabel('average off time (minutes)')
%     legend(plots,labelCell2{:},'Orientation','horizontal','Location','south');

    set(gca,'FontSize',14)
%     ylim([10 40])
    set(gca,'Color',[228,221,209]/255) 
    grid on
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.XAxis(1).Color = 'k';
%     xlim(ap_lims)
%     ylim([0 1.8e5])
    ra_fig.InvertHardcopy = 'off';
    set(gcf,'color','w');
    title(labelCell{s});
    
    saveas(ra_fig,[FigurePath simType '_ra_fits.png'])
end  

%%
ra_fig = figure;
cmap = brewermap([],'Set2');
hold on

% plot experimental data trend
errorbar(master_struct(1).sweepInfoBest.reactivation_time/60,master_struct(1).sweepInfoBest.reactivation_cdf,...
                            sweepInfoBest.reactivation_cdf_ste,'Color','k','Capsize',0);
plots = [];
plots(1) = scatter(master_struct(1).sweepInfoBest.reactivation_time/60,master_struct(1).sweepInfoBest.reactivation_cdf,75,'MarkerFaceColor',cmap(8,:),'MarkerEdgeColor','k');

for s = 1:length(simTypeCell)
%     plot(master_struct(s).sweepInfoBest.reactivation_time/60, master_struct(s).sweepInfoBest.ra_time_cdf_predicted(1:n_plot,:)','Color',[cmap(s,:) 0.5])
    plots(end+1) = plot(master_struct(s).sweepInfoBest.reactivation_time/60, master_struct(s).sweepInfoBest.ra_time_cdf_predicted(1,:),'Color',[cmap(s,:) 1],'LineWidth',2);  
end  
                          

ylabel('cumulative fraction reactivated');
xlabel('reactivation time (minutes)')
legend(plots,labelCell2{:},'Orientation','horizontal','Location','south');

set(gca,'FontSize',14)

set(gca,'Color',[228,221,209]/255) 
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
% xlim(ap_lims)
% ylim([10 40])
ra_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(ra_fig,[FigurePath 'ra_time_fits.png'])

%%
still_on_fig = figure;
cmap = brewermap([],'Set2');
hold on

% plot experimental data trend
plot(master_struct(1).sweepInfoBest.knirps_axis_still_on,master_struct(1).sweepInfoBest.fraction_still_on,...
                            'Color','k');
plots = [];
plots(1) = scatter(master_struct(1).sweepInfoBest.knirps_axis_still_on,master_struct(1).sweepInfoBest.fraction_still_on,75,'MarkerFaceColor',cmap(8,:),'MarkerEdgeColor','k');

for s = 1:length(simTypeCell)
    plot(master_struct(s).sweepInfoBest.knirps_axis_still_on, master_struct(s).sweepInfoBest.p_still_on_predicted(1:n_plot,:)','Color',[cmap(s,:) 0.5])
    plots(end+1) = plot(master_struct(s).sweepInfoBest.knirps_axis_still_on(4:end), master_struct(s).sweepInfoBest.p_still_on_predicted(1,4:end),'Color',[cmap(s,:) 1],'LineWidth',2);  
end  
                          

xlabel('knirps concentration');
ylabel('fraction still on')
legend(plots,labelCell2{:},'Orientation','vertical','Location','northeast');

set(gca,'FontSize',14)

set(gca,'Color',[228,221,209]/255) 
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
% xlim(ap_lims)
% ylim([10 40])
still_on_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(still_on_fig,[FigurePath 'fraction_still_on_fits.png'])

%% Make parameter fit scatters 
n_plot = 100;
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RA vs mean rate
symbolCell = {'o','s','d','^'};
xb = [0 1.5];
yb = [0 0.4];
% calculate experimental envelope 
mf_envelope = 3*mean((master_struct(1).sweepInfoBest.mean_fluo_ap_ste./nanmean(master_struct(1).sweepInfoBest.mean_fluo_ap)));
ra_envelope = 3*mean((master_struct(1).sweepInfoBest.reactivation_cdf_ste./nanmean(master_struct(1).sweepInfoBest.reactivation_cdf)));

ra_vs_mr = figure;
hold on
p = [];
for s = length(simTypeCell):-1:1
    p(end+1) = scatter(sqrt(master_struct(s).sweepInfoBest.ra_fit_R2),sqrt(master_struct(s).sweepInfoBest.off_time_fit_R2),25,symbolCell{s},...
      'MarkerFaceColor',cmap(s,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',0.75);
end  
p = fliplr(p);

plot(linspace(xb(1),xb(2)),repelem(mf_envelope,100),'--k','LineWidth',2)
plot(repelem(ra_envelope,100),linspace(yb(1),yb(2)),'--k','LineWidth',2)

xlim(xb)
ylim(yb)

xlabel('reactivation time fit error');
ylabel('mean fluorescence fit error')

legend(p,labelCell{:},'Orientation','horizontal','Location','south');

set(gca,'FontSize',14)

set(gca,'Color',[228,221,209]/255) 
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

grid on
ra_vs_mr.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(ra_vs_mr,[FigurePath 'ra_vs_mr_scatter.png'])

% set(gca,'xscale','log');
% set(gca,'yscale','log');

%%
ot_envelope = 3*mean((master_struct(1).sweepInfoBest.off_time_ap_ste./nanmean(master_struct(1).sweepInfoBest.off_time_ap)));

ra_vs_ot = figure;
yb = [0 0.4];
hold on
p = [];
for s = length(simTypeCell):-1:1
    p(end+1) = scatter(sqrt(master_struct(s).sweepInfoBest.ra_fit_R2),sqrt(master_struct(s).sweepInfoBest.mean_fluo_fit_R2),25,symbolCell{s},...
      'MarkerFaceColor',cmap(s,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.75);
end  
p = fliplr(p);

plot(linspace(xb(1),xb(2)),repelem(ot_envelope,100),'--k','LineWidth',2)
plot(repelem(ra_envelope,100),linspace(yb(1),yb(2)),'--k','LineWidth',2)

xlim(xb)
% ylim(yb)

xlabel('reactivation time fit error');
ylabel('mean off time fit error')

legend(p,labelCell{:},'Orientation','horizontal','Location','south');

set(gca,'FontSize',14)

set(gca,'Color',[228,221,209]/255) 
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

grid on
ra_vs_ot.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(ra_vs_ot,[FigurePath 'ra_vs_ot_scatter.png'])


%% Re-plot showing only networks that made the cut for fluorescence

mf_vs_ot = figure;
xb = [0 0.4];
yb = [0 0.4];
hold on
p = [];
for s = length(simTypeCell):-1:1    
    p(end+1) = scatter(sqrt(master_struct(s).sweepInfoBest.off_time_fit_R2),sqrt(master_struct(s).sweepInfoBest.mean_fluo_fit_R2),25,symbolCell{s},...
      'MarkerFaceColor',cmap(s,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.75);
end  
p = fliplr(p);

plot(linspace(xb(1),xb(2)),repelem(ot_envelope,100),'--k','LineWidth',2)
plot(repelem(mf_envelope,100),linspace(yb(1),yb(2)),'--k','LineWidth',2)

xlim(xb)
% ylim(yb)

xlabel('mean fluorescence fit error')
ylabel('mean off time fit error')

legend(p,labelCell{:},'Orientation','horizontal','Location','south');

set(gca,'FontSize',14)

set(gca,'Color',[228,221,209]/255) 
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

grid on
mf_vs_ot.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(mf_vs_ot,[FigurePath 'mf_vs_ot_scatter.png'])

%% What happens if we take the best overall score?
total_score_cell = cell(1,length(simTypeCell));
best_score_array = NaN(50,length(simTypeCell));
best_score_ids = NaN(50,length(simTypeCell));

for s = 1:length(simTypeCell)
    total_score_cell{s} = sqrt(master_struct(s).sweepInfoBest.off_time_fit_R2 + master_struct(s).sweepInfoBest.mean_fluo_fit_R2 + master_struct(s).sweepInfoBest.ra_fit_R2);
    [scores_sorted,score_ids] = sort(total_score_cell{s},'ascend');
    best_score_array(:,s) = scores_sorted(1:50);
    best_score_ids(:,s) = score_ids(1:50);
end  

% make scatter plot of results

score_fig = figure;
hold on
p = [];
for s = 1:length(simTypeCell)
    score_vec = total_score_cell{s};
    jitter = rand(size(score_vec))*0.25;
    p(end+1) = scatter(s+jitter,score_vec,25,symbolCell{s},...
      'MarkerFaceColor',cmap(s,:),'MarkerEdgeColor','k','MarkerFaceAlpha',0.75);
end  
% p = fliplr(p);


xlim([0.5,4.5])
% ylim([1e-1 3])

% xlabel('mean fluorescence fit error')
ylabel('aggregate error')
set(gca,'xtick',1:4,'xticklabels',labelCell)
xtickangle(30);
legend(p,labelCell{:},'Orientation','horizontal','Location','north');

set(gca,'FontSize',14)

set(gca,'Color',[228,221,209]/255) 
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
% set(gca,'yscale','log');

saveas(score_fig,[FigurePath 'aggregate_scores.png'])

%% identify 10 best networks and plot
id_array = repmat(1:4,50,1);
id_vec = id_array(:);
best_id_vec = best_score_ids(:);


[scores_ranked,si] = sort(best_score_array(:),'ascend');
ids_ranked = best_id_vec(si);
sim_ids_ranked = id_vec(si);

n_take = 10;
bestGlobalFits = struct;
for n = 1:n_take
  
    % basic info
    bestGlobalFits(n).sim_id = sim_ids_ranked(n);
    bestGlobalFits(n).sim_sub_id = ids_ranked(n);
    bestGlobalFits(n).score = scores_ranked(n);
    
    % fit profiles
    bestGlobalFits(n).param_val_vec = master_struct(sim_ids_ranked(n)).sweepInfoBest.param_val_vec(ids_ranked(n),:);
    bestGlobalFits(n).mean_fluo_predicted = master_struct(sim_ids_ranked(n)).sweepInfoBest.mean_fluo_predicted(ids_ranked(n),:);
    bestGlobalFits(n).off_time_predicted = master_struct(sim_ids_ranked(n)).sweepInfoBest.off_time_predicted(ids_ranked(n),:);
    bestGlobalFits(n).ra_time_cdf = master_struct(sim_ids_ranked(n)).sweepInfoBest.ra_time_cdf_predicted(ids_ranked(n),:);
    bestGlobalFits(n).p_still_on_predicted = master_struct(sim_ids_ranked(n)).sweepInfoBest.p_still_on_predicted(ids_ranked(n),:);
    
    % other info
    bestGlobalFits(n).ms2_traces_observed_wt = master_struct(sim_ids_ranked(n)).sweepInfoBest.ms2_traces_observed_wt(:,:,ids_ranked(n));
    bestGlobalFits(n).knirps_traces_wt = master_struct(sim_ids_ranked(n)).sweepInfoBest.knirps_traces_wt(:,:,ids_ranked(n));
    bestGlobalFits(n).knirps_traces_ra = master_struct(sim_ids_ranked(n)).sweepInfoBest.knirps_traces_ra(:,:,ids_ranked(n));
    bestGlobalFits(n).tf_dependent_curves_wt = master_struct(sim_ids_ranked(n)).sweepInfoBest.tf_dependent_curves_wt(:,:,ids_ranked(n));
    bestGlobalFits(n).tf_dependent_curves_ra = master_struct(sim_ids_ranked(n)).sweepInfoBest.tf_dependent_curves_ra(:,:,ids_ranked(n));
    bestGlobalFits(n).ms2_traces_observed_ra = master_struct(sim_ids_ranked(n)).sweepInfoBest.ms2_traces_observed_ra(:,:,ids_ranked(n));
    bestGlobalFits(n).trace_ap_vec_wt = master_struct(sim_ids_ranked(n)).sweepInfoBest.trace_ap_vec_wt(ids_ranked(n),:);
end

%% for the 10 best fits, let's re-calculate the "still on" response curve
ap_limits = master_struct(1).sweepInfo.ap_limits_still_on;

for n = 1:n_take
    % take subset of traces from relevant region
    trace_ap_vec = bestGlobalFits(n).trace_ap_vec_wt;
    ap_indices = find(trace_ap_vec>=ap_limits(1) & trace_ap_vec<=ap_limits(2));
    fluo_array_filtered = bestGlobalFits(n).ms2_traces_observed_wt(:,ap_indices);
    knirps_array = bestGlobalFits(n).knirps_traces_wt(:,ap_indices);
    
    % generate indexing vectors
    index_vec = (1:size(fluo_array_filtered,1))';
    
    % calculate first and last active frames for each trace
    active_indices = 1*(fluo_array_filtered>0) .* index_vec;    
    active_indices2 = active_indices;
    active_indices2(active_indices2==0) = Inf;
    first_i_vec = min(active_indices2);
    last_i_vec = max(active_indices);        
    
    still_on_array = ones(size(fluo_array_filtered));
    for i = 1:size(fluo_array_filtered,2)
        if ~isinf(first_i_vec(i))
%             still_on_array(first_i_vec(i):end,i) = 1;
            still_on_array(last_i_vec(i):end,i) = 0;
        end
    end         
    still_on_vec = still_on_array(:);
%     knirps_array = permute(tf_profile_array(:,:,ap_indices),[1 3 2]);
    knirps_groups = discretize(knirps_array(:),master_struct(1).sweepInfo.knirps_bins_still_on);
    first_group = nanmin(knirps_groups);
    predicted_fraction_still_on = NaN(size(master_struct(1).sweepInfo.knirps_axis_still_on));
    for k = 1:length(predicted_fraction_still_on)
        val = nanmean(still_on_vec(knirps_groups==k));
        if sum(knirps_groups==k)>=10
            predicted_fraction_still_on(k) = val;
        end
    end 
    predicted_fraction_still_on(1:first_group-1) = 1;
    bestGlobalFits(n).predicted_fraction_still_on = predicted_fraction_still_on;
end
%%
still_on_fig = figure;
cmap = brewermap([],'Set2');
hold on

% plot experimental data trend
plot(master_struct(1).sweepInfoBest.knirps_axis_still_on,master_struct(1).sweepInfoBest.fraction_still_on,...
                            'Color','k');
plots = [];
scatter(master_struct(1).sweepInfoBest.knirps_axis_still_on,master_struct(1).sweepInfoBest.fraction_still_on,75,'MarkerFaceColor',cmap(8,:),'MarkerEdgeColor','k');
sim_id_vec = [];

xlabel('knirps concentration');
ylabel('fraction still on')


set(gca,'FontSize',14)

set(gca,'Color',[228,221,209]/255) 
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
% xlim(ap_lims)
% ylim([10 40])
still_on_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(still_on_fig,[FigurePath 'fraction_still_on_data_only.png'])

for s = 1:length(bestGlobalFits)
    sim_id_vec(end+1) = bestGlobalFits(s).sim_id;
%     plot(master_struct(s).sweepInfoBest.knirps_axis_still_on, master_struct(s).sweepInfoBest.p_still_on_predicted(1:n_plot,:)','Color',[cmap(s,:) 0.5])
    plots(end+1) = plot(master_struct(1).sweepInfoBest.knirps_axis_still_on, bestGlobalFits(s).predicted_fraction_still_on,'Color',[cmap(bestGlobalFits(s).sim_id,:) 1],'LineWidth',2);  
end  
[unique_ids,ia] =  unique(sim_id_vec);                         
plots = plots(ia);
legend(plots,labelCell{unique_ids},'Orientation','vertical','Location','southwest','Color','w');

saveas(still_on_fig,[FigurePath 'fraction_still_on_best_fits.png'])

%%
miss_fig = figure;
plot(master_struct(1).sweepInfoBest.fluo_ref_curve, master_struct(1).sweepInfoBest.p_miss_ref_vec,'-k','LineWidth',2)

xlabel('true spot fluorescence (au)');
ylabel('likelihood of missed detection')

set(gca,'FontSize',14)

set(gca,'Color',[228,221,209]/255) 
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
xlim([0 5e5])
% ylim([10 40])
miss_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(miss_fig,[FigurePath 'miss_prob_plot.png'])

%% Retroactively calculate combined fit scores

bins = linspace(0,2,25);
figure;
hold on
histogram(master_struct(1).sweepInfoBest.bestScore,bins,'Normalization','probability')
histogram(master_struct(2).sweepInfoBest.bestScore,bins,'Normalization','probability')
histogram(master_struct(4).sweepInfoBest.bestScore,bins,'Normalization','probability')
