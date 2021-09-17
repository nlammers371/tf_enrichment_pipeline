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
try
    FigurePath = [liveProject.figurePath 'io_sim_results' filesep];
catch
    FigurePath = 'C:\Users\nlamm\Dropbox (Personal)\LocalEnrichmentFigures\PipelineOutput\paramSweeps\';
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
    readPath = [resultsRoot 'sweeps_n' num2str(nInc) filesep];
    % load sweep result
    load([readPath 'sweepInfo_' simType '.mat'],'sweepInfo');
    
    % store
    master_struct(s).sweepInfo = sweepInfo;
end

% load raw data
load([resultsRoot 'io_ref_ra.mat'],'io_ref_ra')
load([resultsRoot 'io_ref_wt.mat'],'io_ref_wt')

%% obtain predicted trends for best systems from sim type
n_traces = 250;
n_best = 10;

for s = 2%1:length(master_struct)
    % extract
    sweepInfo = master_struct(s).sweepInfo;
    simType = sweepInfo.simType;
    
    % calculate aggregate score
    total_score = -sqrt(sweepInfo.fluo_time_fit.^2 + sweepInfo.ra_fit.^2 + sweepInfo.still_on_fit.^2);

    % find best overall performers
    [~,score_ids] = sort(total_score,'descend');
    best_i_list = score_ids(1:n_best);
    
    % now find best parameter-specific scores
    [~,best_ft_i] = nanmax(sweepInfo.fluo_time_fit);
    [~,best_pon_i] = nanmax(sweepInfo.still_on_fit);
    [~,best_ra_i] = nanmax(sweepInfo.ra_fit);
    
    % get corresponding parameters
    param_vec = sweepInfo.param_val_vec([best_i_list',best_ft_i,best_ra_i,best_pon_i],:);
    
    % run sweep for selected networks
    sweepInfoBest = io_sweep_wrapper(resultsRoot,2,sweepInfo.simType,param_vec,true,'n_traces',n_traces);
    
    % store
    master_struct(s).sweepInfoBest = sweepInfoBest;
    master_struct(s).best_i_list = best_i_list;
    master_struct(s).best_ft_i = best_ft_i;
    master_struct(s).best_ra_i = best_ra_i;
    master_struct(s).best_pon_i = best_pon_i;
end
%% Make plots
metricPath = [FigurePath 'metricPlots' filesep];
mkdir(metricPath)
labelCell = {'k_{off}','k_{on}','k_{s}','k_{a}'};
labelCell2 = {'global optima','mean fluo','reactivation time','silencing'};
labelCell3 = [{'data'} labelCell2];

for s = 1:length(simTypeCell)
  
    % Silencing fig
    still_on_fig = figure;
    cmap = brewermap([],'Set2');
    hold on
    simType = simTypeCell{s};
    plots = [];
    
    % plot experimental data trend
    errorbar(master_struct(s).sweepInfoBest.knirps_axis_still_on,master_struct(s).sweepInfoBest.fraction_still_on,...
                                master_struct(s).sweepInfoBest.fraction_still_on_ste,'Color','k','Capsize',0);    
                              
    plots(1) = scatter(master_struct(s).sweepInfoBest.knirps_axis_still_on,master_struct(s).sweepInfoBest.fraction_still_on,...
                        75,'MarkerFaceColor',cmap(8,:),'MarkerEdgeColor','k');
    
    xlabel('knirps concentration (au)');
    ylabel('fraction still on')
    title(labelCell{s})  

    set(gca,'FontSize',14)

    set(gca,'Color',[228,221,209]/255) 
    grid on
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.XAxis(1).Color = 'k';
%     xlim(ap_lims)
    ylim([0 1.05])
    still_on_fig.InvertHardcopy = 'off';
    set(gcf,'color','w');
    
    if s == 1
        saveas(still_on_fig,[metricPath 'still_on_dadta_only.png'])
    end
    title(labelCell{s});
    
    

    % plot best N performers
    plot(master_struct(s).sweepInfoBest.knirps_axis_still_on, master_struct(s).sweepInfoBest.p_still_on_predicted(1:n_best,:)','Color',[cmap(2,:) 0.5])
    
    % plot best rep for each metric
    plots(2) = plot(master_struct(s).sweepInfoBest.knirps_axis_still_on, master_struct(s).sweepInfoBest.p_still_on_predicted(1,:)','Color',cmap(2,:),'LineWidth',2);
    plots(3) = plot(master_struct(s).sweepInfoBest.knirps_axis_still_on, master_struct(s).sweepInfoBest.p_still_on_predicted(end-2,:)','--','Color',cmap(3,:),'LineWidth',2);
    plots(4) = plot(master_struct(s).sweepInfoBest.knirps_axis_still_on, master_struct(s).sweepInfoBest.p_still_on_predicted(end-1,:)','--','Color',cmap(4,:),'LineWidth',2);
    plots(5) = plot(master_struct(s).sweepInfoBest.knirps_axis_still_on, master_struct(s).sweepInfoBest.p_still_on_predicted(end,:)','-','Color',cmap(5,:),'LineWidth',2);
        
    legend(plots,labelCell3{:},'Location','southwest','Color','w');
    
    saveas(still_on_fig,[metricPath simType '_still_on_fits.png'])
    
end  

%%
close all

% Mean fluo vs time
for s = 1:length(simTypeCell)
  
    % Silencing fig
    fluo_time_fig = figure;
    cmap = brewermap([],'Set2');
    hold on
    simType = simTypeCell{s};
    plots = [];
    
    % plot experimental data trendsweepInfoBest.time_axis_mf,sweepInfoBest.fluo_time_mean,'-k')
    errorbar(master_struct(s).sweepInfoBest.time_axis_mf,master_struct(s).sweepInfoBest.fluo_time_mean,...
                                master_struct(s).sweepInfoBest.fluo_time_ste,'Color','k','Capsize',0);                              
    
    plots(1) = scatter(master_struct(s).sweepInfoBest.time_axis_mf,master_struct(s).sweepInfoBest.fluo_time_mean,...
                        75,'MarkerFaceColor',cmap(8,:),'MarkerEdgeColor','k');

    
    xlabel('time into nc14 (minutes)');
    ylabel('average active spot fluorescence (au)')
    title(labelCell{s})    

    set(gca,'FontSize',14)

    set(gca,'Color',[228,221,209]/255) 
    grid on
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.XAxis(1).Color = 'k';
    xlim([8 25])
    ylim([0 3e5])
    fluo_time_fig.InvertHardcopy = 'off';
    set(gcf,'color','w');
    
    if s == 1
        saveas(fluo_time_fig,[metricPath 'fluo_time_data_only.png'])
    end
    
    % plot best N performers
    plot(master_struct(s).sweepInfoBest.time_axis_mf', permute(master_struct(s).sweepInfoBest.fluo_time_predicted(:,1,1:n_best),[1 3 2])','Color',[cmap(2,:) 0.5])
    
    % plot best rep for each metric
    plots(2) = plot(master_struct(s).sweepInfoBest.time_axis_mf, master_struct(s).sweepInfoBest.fluo_time_predicted(:,:,1),'Color',cmap(2,:),'LineWidth',2);
    plots(3) = plot(master_struct(s).sweepInfoBest.time_axis_mf, master_struct(s).sweepInfoBest.fluo_time_predicted(:,:,end-2),'-','Color',cmap(3,:),'LineWidth',2);
    plots(4) = plot(master_struct(s).sweepInfoBest.time_axis_mf, master_struct(s).sweepInfoBest.fluo_time_predicted(:,:,end-1),'--','Color',cmap(4,:),'LineWidth',2);
    plots(5) = plot(master_struct(s).sweepInfoBest.time_axis_mf, master_struct(s).sweepInfoBest.fluo_time_predicted(:,:,end),'--','Color',cmap(5,:),'LineWidth',2);

    legend(plots,labelCell3{:},'Location','northeast','Color','w');
    
    title(labelCell{s});
    
    saveas(fluo_time_fig,[metricPath simType '_mean_fluo_time.png'])
    
end  
%%

% Reactivation time CDF
for s = 1:length(simTypeCell)
    ra_time_axis = master_struct(s).sweepInfoBest.reactivation_time/60;
    
    ra_fig = figure;
    cmap = brewermap([],'Set2');
    hold on
    simType = simTypeCell{s};
    plots = [];
    
    % plot experimental data trendsweepInfoBest.time_axis_mf,sweepInfoBest.fluo_time_mean,'-k')
    errorbar(ra_time_axis,master_struct(s).sweepInfoBest.reactivation_cdf,...
                                master_struct(s).sweepInfoBest.reactivation_cdf_ste,'Color','k','Capsize',0);                              
    
    plots(1) = scatter(ra_time_axis,master_struct(s).sweepInfoBest.reactivation_cdf,...
                            75,'MarkerFaceColor',cmap(8,:),'MarkerEdgeColor','k');   

    
    xlabel('time since perturbation');
    ylabel('average active spot fluorescence (au)')
    title(labelCell{s})    

    set(gca,'FontSize',14)

    set(gca,'Color',[228,221,209]/255) 
    grid on
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.XAxis(1).Color = 'k';
    xlim([ra_time_axis(1) ra_time_axis(end)])
    ylim([0 1.05])
    ra_fig.InvertHardcopy = 'off';
    set(gcf,'color','w');
    
    if s == 1
        saveas(ra_fig,[metricPath 'ra_cdf_data_only.png'])
    end
    
    % plot best N performers
    plot(ra_time_axis, master_struct(s).sweepInfoBest.ra_time_cdf_predicted(1:n_best,:),'Color',[cmap(2,:) 0.5])
    
    % plot best rep for each metric
    plots(2) = plot(ra_time_axis, master_struct(s).sweepInfoBest.ra_time_cdf_predicted(1,:),'Color',cmap(2,:),'LineWidth',2);
    plots(3) = plot(ra_time_axis, master_struct(s).sweepInfoBest.ra_time_cdf_predicted(end-2,:),'--','Color',cmap(3,:),'LineWidth',2);
    plots(4) = plot(ra_time_axis, master_struct(s).sweepInfoBest.ra_time_cdf_predicted(end-1,:),'-','Color',cmap(4,:),'LineWidth',2);
    plots(5) = plot(ra_time_axis, master_struct(s).sweepInfoBest.ra_time_cdf_predicted(end,:),'--','Color',cmap(5,:),'LineWidth',2);
    
    title(labelCell{s});
    
    legend(plots,labelCell3{:},'Location','southeast','Color','w');
    
    saveas(ra_fig,[metricPath simType '_reactivation_time.png'])
    
end  

%% Make figure illustrating probabilistic threshold

thresh_fig = figure;
hold on
plot(sweepInfoBest.ms2_traces_true_wt(:,1,1),'-o','LineWidth',2,'Color',cmap(8,:))
plot(sweepInfoBest.ms2_traces_observed_wt(:,1,1),'-','LineWidth',2,'Color','k')

xlabel('time');
ylabel('spot fluorescence (au)')
title(labelCell{s})
legend('ground truth','observed','Location','northeast');

set(gca,'FontSize',14)

set(gca,'Color',[228,221,209]/255) 
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
thresh_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(thresh_fig,[metricPath simType '_thresh_fig.png'])