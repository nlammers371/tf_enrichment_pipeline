% Script to scan through parameter space for selected biological variables
clear
close all

addpath(genpath('../utilities'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%% Load sweep results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

projectName = 'optokni_eve4+6_ON'; 
try
    liveProject = LiveEnrichmentProject(projectName);
    resultsRoot = [liveProject.dataPath filesep];
catch
    resultsRoot = 'C:\Users\nlamm\Dropbox (Personal)\ProcessedEnrichmentData\optokni_eve4+6_ON\';%[liveProject.dataPath filesep];
end


simTypeCell = {'out_only','in_only','koff_only_2','kon_only_2'};
labelCell = {'k_{s}','k_{a}','k_{off}','k_{on}'};
% simTypeCell = {'koff_only_2','kon_only_2'};
% labelCell = {'k_{off}','k_{on}'};
master_struct = struct;
for s = 1:length(simTypeCell)
    % load and store sweep result sets
    simType = simTypeCell{s};
    load([resultsRoot 'sweepInfo_' simType '.mat'],'sweepInfo');
    master_struct(s).sweepInfo = sweepInfo;
    clear sweepInfo
end    

% load basic data structure containing reference metrics
load([resultsRoot 'io_ref_struct.mat']);

% mak figure directory  
try
    FigurePath = [liveProject.figurePath 'io_sim_results' filesep];
catch
    FigurePath = 'C:\Users\nlamm\Dropbox (Personal)\LocalEnrichmentFigures\PipelineOutput\optokni_eve4+6_ON\';
end
mkdir(FigurePath)

%% Make basic heatmap to illustrate knirps repression
t_vec = io_ref_struct.time_vec;

hm_fig = figure;
cmap = flipud(brewermap([],'Spectral'));
colormap(cmap);
imagesc(flipud(io_ref_struct.fluo_array));

h = colorbar;

xlabel('trace ID');
ylabel('time since ON perturbation')
ylabel(h,'spot fluorescence (au)')

set(gca,'ytick',1:10:length(t_vec),'yticklabels',round(t_vec(1:10:length(t_vec))/60))
% grid on
set(gca,'FontSize',14)

saveas(hm_fig,[FigurePath 'illustrative_heatmap.png'])

%% Make fluo detection threshold figure
close all

gauss_fit_parameters = io_ref_struct.fit_parameters;
bins_fit = io_ref_struct.min_fluo_bins;
gauss_fun = io_ref_struct.gauss_fun;
min_fluo_counts = io_ref_struct.min_fluo_counts;

% get predicted profile
profile_pd = gauss_fun(gauss_fit_parameters);

detection_fig = figure;
cmap = brewermap([],'Set2');
hold on

b = bar(bins_fit,min_fluo_counts,1,'FaceColor',cmap(5,:));
p = plot(bins_fit,profile_pd,'Color','k','LineWidth',2);

xlabel('minimum detected spot intensity (au)');
ylabel('frequency')
% grid on
set(gca,'FontSize',14)

legend([b,p],'data',['Gaussian fit (\mu=' num2str(round(gauss_fit_parameters(1),-2)) ')'])
set(gca,'Color',[228,221,209]/255) 
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

saveas(detection_fig,[FigurePath 'detection_threshold_fit.png'])

%% Next, look at reactivation waiting times
labelCell2 = [{'data'} labelCell];

reactivation_time_axis = io_ref_struct.reactivation_time_axis;
wt_pdf_true = histcounts(io_ref_struct.reactivation_time_vec*60,[reactivation_time_axis reactivation_time_axis(end)+20],'Normalization','probability');
wt_cdf_true = io_ref_struct.reactivation_time_cdf;
wt_cdf_true(end) = 1;

% find best fits using CDF
ra_cdf_times = 0:20:50*60;
n_fits = 1e3;
best_ra_fit_ids = NaN(n_fits,length(master_struct));
best_ra_fit_scores = NaN(n_fits,length(master_struct));
best_ra_cdf_fits = NaN(length(reactivation_time_axis),n_fits,length(master_struct));
best_ra_pdf_fits = NaN(length(reactivation_time_axis),n_fits,length(master_struct));

for s = 1:length(master_struct)
    % extract and filter
    ra_cdf_array = master_struct(s).sweepInfo.reactivation_time_cdf_array;    
    time_filter = ismember(ra_cdf_times,reactivation_time_axis);
    ra_cdf_array_ft = ra_cdf_array(time_filter,:);
    % calculate msq diff
    mean_squared_differences = sqrt(mean((ra_cdf_array_ft-wt_cdf_true').^2));
    master_struct(s).mean_squared_differences = mean_squared_differences;
    % sort results according to fit
    [~ ,fit_ranking] = sort(master_struct(s).mean_squared_differences);
    % decide how many to take, take n_fits if there are enough that meet
    % cirterio
    n_quality = nansum(mean_squared_differences<=2*nanmin(mean_squared_differences));
    n_keep = min([n_quality n_fits]);
    % take top N IDs
    best_ra_fit_ids(1:n_keep,s) = fit_ranking(1:n_keep);
    best_ra_fit_scores(1:n_keep,s) = mean_squared_differences(fit_ranking(1:n_keep));
    % And top N profiles for CDF
    best_ra_cdf_fits(:,1:n_keep,s) = ra_cdf_array_ft(:,fit_ranking(1:n_keep));      
end
 
close all
% make plots for each individual set
for s = 1:length(master_struct)
    cdf_fit_figure = figure;
    hold on
    % plot fits
    plot(reactivation_time_axis/60,best_ra_cdf_fits(:,:,s),'Color',[cmap(s,:) .3])
    p1 = plot(reactivation_time_axis/60,best_ra_cdf_fits(:,1,s),'-.','Color',cmap(s,:),'LineWidth',3);
    
    % plot ground truth    
    p2 = plot(reactivation_time_axis/60,wt_cdf_true,'Color','k','LineWidth',3);
    
    xlabel('minutes since ON perturbation');
    ylabel('cumulative fraction reactivated')
    
    % grid on
    set(gca,'FontSize',14)

    legend([p2 p1],'data',['best fit (' labelCell{s} ')'],'Location','southeast')
    set(gca,'Color',[228,221,209]/255) 
    grid on
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.XAxis(1).Color = 'k';
    xlim([0 8])
    
    cdf_fit_figure.InvertHardcopy = 'off';
    set(gcf,'color','w');

    saveas(cdf_fit_figure,[FigurePath 'cdf_fit_fig_' simTypeCell{s} '.png'])
    
end  

% make pdf plot 
pdf_fig = figure;
bar(reactivation_time_axis/60,wt_pdf_true,1,'FaceColor',cmap(8,:))

xlabel('minutes since ON perturbation');
ylabel('fraction reactivated')

set(gca,'FontSize',14)

set(gca,'Color',[228,221,209]/255) 
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
xlim([0 8])

pdf_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(pdf_fig,[FigurePath 'pdf_fig_data_only.png'])

% plot them all together
p = [];
cdf_fig = figure;

hold on

% plot ground truth    
p(1) = plot(reactivation_time_axis/60,wt_cdf_true,'Color','k','LineWidth',3);

xlabel('minutes since ON perturbation');
ylabel('cumulative fraction reactivated')

set(gca,'FontSize',14)

set(gca,'Color',[228,221,209]/255) 
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
xlim([0 8])

cdf_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
saveas(cdf_fig,[FigurePath 'cdf_fit_fig_data_only.png'])

for s = 1:length(master_struct)
%     plot(reactivation_time_axis/60,best_ra_cdf_fits(:,:,s),'Color',[cmap(s,:) .3])
    p(s+1) = plot(reactivation_time_axis/60,best_ra_cdf_fits(:,1,s),'-.','Color',cmap(s,:),'LineWidth',2);
end

legend(labelCell2{:},'Location','southeast','Color','w')
saveas(cdf_fig,[FigurePath 'cdf_fit_fig_all.png'])

%% Generate example traces for best-fitting parameter set in each model scenario
swap_fields = {'n_traces','granularity','systemParams','simType','paramList','fitFlags','tf_profile_array'};
rng(123)
close all
n_plot = 10;
% pt_time = 41*20 / 60;

for s = 1:length(master_struct)
    best_id = best_ra_fit_ids(1,s);
    % generate temporay structure
    simTemp = struct;
    simTemp.param_fit_array = master_struct(s).sweepInfo.param_fit_array(best_id,:);
    simTemp.step = 1;
    for f = 1:length(swap_fields)
        simTemp.(swap_fields{f}) = master_struct(s).sweepInfo.(swap_fields{f});
    end
    
    % run simulation
    simInfoPD = io_prediction_wrapper_v3(simTemp);
    F_min = master_struct(s).sweepInfo.param_fit_array(best_id,3)';
    % store
    master_struct(s).simExample = simInfoPD;
    master_struct(s).sweepInfoTrunc = simTemp;
    
    % plot results
    % look for traces that turned off
    off_flags = ~isnan(simInfoPD.reactivation_time_vec);
    
    plot_ids = randsample(1:simTemp.n_traces,n_plot,true,off_flags); 
    
    % generate time axis
    time_axis = simInfoPD.gillespie.t_ref;
    time_axis = (time_axis - time_axis(41))/60;
    time_filter_plot = time_axis >= -10 & time_axis <= 10; 
    
    % generate mean knirps profile
    mean_rate_vec = nanmean(simInfoPD.gillespie.rate_curve_in,3);
    mean_knirps_vec = nanmean(simInfoPD.gillespie.tf_ref_in,3);
    
    trace_fig_knirps = figure;
%     colormap(cmap)
    hold on
    plot(time_axis(time_filter_plot),simInfoPD.gillespie.fluo_ms2_array(time_filter_plot,plot_ids))    
    xlabel('minutes since ON perturbation');
    ylabel('simulated fluorescence (au)')

    yyaxis right        
    plot(time_axis(time_filter_plot),mean_knirps_vec(time_filter_plot),'Color',cmap(5,:),'LineWidth',2)
    set(gca,'FontSize',14)
    
    ylabel('knirps concentration (au)')
    
    set(gca,'Color',[228,221,209]/255) 
    grid on
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.XAxis(1).Color = 'k';
    ax.YAxis(2).Color = cmap(5,:);
    xlim([-10 10])
    
    trace_fig_knirps.InvertHardcopy = 'off';
    set(gcf,'color','w');
    set(gca,'Clipping','on')
    saveas(trace_fig_knirps,[FigurePath simTypeCell{s} '_sim_w_knirps'  '.png'])
    
    yyaxis left
    plot(time_axis(time_filter_plot),repelem(F_min,sum(time_filter_plot)),'--k','LineWidth',1.5)
    saveas(trace_fig_knirps,[FigurePath simTypeCell{s} '_sim_w_knirps'  '_thresh.png'])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % now with rate instead
    trace_fig_rate = figure;
%     colormap(cmap)
    hold on
    plot(time_axis(time_filter_plot),simInfoPD.gillespie.fluo_ms2_array(time_filter_plot,plot_ids))    
    xlabel('minutes since ON perturbation');
    ylabel('simulated fluorescence (au)')

    yyaxis right        
    plot(time_axis(time_filter_plot),mean_rate_vec(time_filter_plot),'Color',cmap(s,:),'LineWidth',2)
    set(gca,'FontSize',14)
    
    ylabel([labelCell{s} ' (s^{-1})'])
    
    set(gca,'Color',[228,221,209]/255) 
    grid on
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.XAxis(1).Color = 'k';
    ax.YAxis(2).Color = cmap(s,:);
    xlim([-10 10])
    
    trace_fig_rate.InvertHardcopy = 'off';
    set(gcf,'color','w');
    set(gca,'Clipping','on')
    saveas(trace_fig_rate,[FigurePath  simTypeCell{s} '_sim_w_rate' '.png'])
    
    yyaxis left
    plot(time_axis(time_filter_plot),repelem(F_min,sum(time_filter_plot)),'--k','LineWidth',1.5)
    saveas(trace_fig_rate,[FigurePath  simTypeCell{s} '_sim_w_rate' '_thresh.png'])
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%% Now attempt to calculate global best fit %%%%%%%%%%%%%%%%%%
pon_true_profile = master_struct(1).sweepInfo.p_on_true;
ff_true_profile = master_struct(1).sweepInfo.fluo_true;

time_vec = master_struct(1).sweepInfo.time_vec/60;
t_bounds_fit = [-5 5]; % demarcates time poiunts that we care about
t_filter = time_vec <= t_bounds_fit(2) & time_vec >= t_bounds_fit(1);

n_keep = 25;

overall_best_fit_array = NaN(n_keep,length(master_struct));

% calculate mean-squared error for "PON" on "Fluo Full" metrics
% Need a way to normalize for comparison with RA cdf ift

for s = 1:length(master_struct)
    %%%%%%%%%%%%%%%%%%%%
    % PON    
    pon_fit_array = master_struct(s).sweepInfo.p_on_fit_array;
    
    pon2_array = (pon_fit_array-pon_true_profile).^2;
    master_struct(s).pon_objective_vec_alt = mean(sqrt(pon2_array(t_filter,:)),1);
    
    pon_mean = nanmean(master_struct(s).pon_objective_vec_alt);
    pon_std = nanstd(master_struct(s).pon_objective_vec_alt);
    
    ob_pon_norm = (master_struct(s).pon_objective_vec_alt)/pon_std;
    
    
    %%%%%%%%%%%%%%%%%%%%
    % Fluo full
    ff_fit_array = master_struct(s).sweepInfo.fluo_fit_array;    
    ff2_array = (ff_fit_array-ff_true_profile).^2;
    master_struct(s).fluo_full_objective_vec_alt = mean(sqrt(ff2_array(t_filter,:)),1);
    
    ff_mean = nanmean(master_struct(s).fluo_full_objective_vec_alt);
    ff_std = nanstd(master_struct(s).fluo_full_objective_vec_alt);
    
    ob_ff_norm = (master_struct(s).fluo_full_objective_vec_alt)/ff_std;
    
    % normalize ra
    cdf_mean = nanmean(master_struct(s).mean_squared_differences);
    cdf_std = nanstd(master_struct(s).mean_squared_differences);
    ob_ra_norm = (master_struct(s).mean_squared_differences)/cdf_std;
    
    % alright...let's compare 
    objective_full = sqrt(ob_ra_norm.^2 + ob_pon_norm.^2 + ob_ff_norm.^2);
    objective_full(isnan(objective_full)) = Inf;
    
    % sort
    [~,best_ids_full] = sort(objective_full);
    overall_best_fit_array(:,s) = best_ids_full(1:n_keep);
    
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%% Make basic fit figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PON

pon_fig = figure;
cmap = brewermap([],'Set2');
hold on
p = [];
% plot the true profile
p(1) = plot(time_vec,pon_true_profile,'Color','k','LineWidth',2.5);

ylabel('instantaneous active fraction')

% yyaxis right
% plot(time_vec,mean_tf_vec,'-','Color',cmap(s+1,:),'LineWidth',2)

% ylim([0 2])
xlim(t_bounds_fit)

xlabel('minutes from light ON event');
% ylabel('average knirps concentration (au)')
% grid on
set(gca,'FontSize',14)

set(gca,'Color',[228,221,209]/255) 
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
% ax.YAxis(2).Color = cmap(s+1,:);
ax.XAxis(1).Color = 'k';

pon_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
ylim([0 1])
saveas(pon_fig,[FigurePath 'pon_true_profile.png'])


% iterate through sim types
for s = 1:length(master_struct)
    
    best_overall_ids = overall_best_fit_array(:,s);
    best_overall_ids = best_overall_ids(~isnan(best_overall_ids));
    pon_fit_profiles = master_struct(s).sweepInfo.p_on_fit_array(:,best_overall_ids);            
    
    % plot%   
    p(s+1) = plot(time_vec,pon_fit_profiles(:,1),'-.','Color',cmap(s,:),'LineWidth',2);
    plot(time_vec,pon_fit_profiles,'-','Color',[cmap(s,:) .2],'LineWidth',1);
%
end

legend(p,labelCell2{:},'Location','southeast','Color','w')

saveas(pon_fig,[FigurePath 'best_p_on_fits_full.png'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fluo Full

fluo_fig = figure;
cmap = brewermap([],'Set2');
hold on
p = [];
% plot the true profile
p(1) = plot(time_vec,ff_true_profile,'Color','k','LineWidth',2.5);

ylabel('average spot fluorescence (all loci) (au)')

% yyaxis right
% plot(time_vec,mean_tf_vec,'-','Color',cmap(s+1,:),'LineWidth',2)

% ylim([0 2])
xlim(t_bounds_fit)

xlabel('minutes from light ON event');
% ylabel('average knirps concentration (au)')
% grid on
set(gca,'FontSize',14)

set(gca,'Color',[228,221,209]/255) 
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
% ax.YAxis(2).Color = cmap(s+1,:);
ax.XAxis(1).Color = 'k';

fluo_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
ylim([0 2.5e5])
saveas(fluo_fig,[FigurePath 'fluo_full_true.png'])

% iterate through sim types
for s = 1:length(master_struct)
    
    best_overall_ids = overall_best_fit_array(:,s);
    best_overall_ids = best_overall_ids(~isnan(best_overall_ids));
    ff_fit_profiles = master_struct(s).sweepInfo.fluo_fit_array(:,best_overall_ids);            
    
    % plot%   
    p(s+1) = plot(time_vec,ff_fit_profiles(:,1),'-.','Color',cmap(s,:),'LineWidth',2);
    plot(time_vec,ff_fit_profiles,'-','Color',[cmap(s,:) .2],'LineWidth',1);
%
end

legend(p,labelCell2{:},'Location','northwest','Color','w')
saveas(fluo_fig,[FigurePath 'best_fluo_fits_full.png'])

% Reactivation time distribution
%% make plots for each individual set
cdf_fit_figure = figure;
hold on

p = [];
p(1) = plot(reactivation_time_axis/60,wt_cdf_true,'Color','k','LineWidth',3); 

ylabel('cumulative fraction reactivated')

% yyaxis right
% plot(time_vec,mean_tf_vec,'-','Color',cmap(s+1,:),'LineWidth',2)

% ylim([0 2])
xlim([0 8])
ylim([0 1.05])
xlabel('minutes from light ON event');
% ylabel('average knirps concentration (au)')
% grid on
set(gca,'FontSize',14)
% legend(p,labelCell2{:},'Location','southeast')
set(gca,'Color',[228,221,209]/255) 
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
% ax.YAxis(2).Color = cmap(s+1,:);
ax.XAxis(1).Color = 'k';

cdf_fit_figure.InvertHardcopy = 'off';
set(gcf,'color','w');
saveas(cdf_fit_figure,[FigurePath 'cdf_true.png'])

for s = 1:length(master_struct)
  
    best_overall_ids = overall_best_fit_array(:,s);
    best_overall_ids = best_overall_ids(~isnan(best_overall_ids));
    ra_fit_profiles = master_struct(s).sweepInfo.reactivation_time_cdf_array(:,best_overall_ids);  
    
    time_filter_cdf = ismember(ra_cdf_times,reactivation_time_axis);
    ra_fit_profiles = ra_fit_profiles(time_filter_cdf,:);
    
    % plot fits
    plot(reactivation_time_axis/60,ra_fit_profiles,'Color',[cmap(s,:) .2])
    p(s+1) = plot(reactivation_time_axis/60,ra_fit_profiles(:,1),'-.','Color',cmap(s,:),'LineWidth',2);
            
end  

legend(p,labelCell2{:},'Location','northwest','Color','w')
saveas(cdf_fit_figure,[FigurePath 'best_cdf_fits_full.png'])