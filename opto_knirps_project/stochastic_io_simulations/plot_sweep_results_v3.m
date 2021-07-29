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
    FigurePath = 'C:\Users\nlamm\Dropbox (Personal)\LocalEnrichmentFigures\PipelineOutput\optokni_eve4+6_ON';
end
mkdir(FigurePath)


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
p = plot(bins_fit,profile_pd,'Color','k');

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
    mean_squared_differences = sqrt(sum((ra_cdf_array_ft-wt_cdf_true').^2));
    % sort results according to fit
    [~ ,fit_ranking] = sort(mean_squared_differences);
    % decide how many to take, take n_fits if there are enough that meet
    % cirterio
    n_quality = nansum(mean_squared_differences<=2*nanmin(mean_squared_differences));
    n_take = min([n_quality n_fits]);
    % take top N IDs
    best_ra_fit_ids(1:n_take,s) = fit_ranking(1:n_take);
    best_ra_fit_scores(1:n_take,s) = mean_squared_differences(fit_ranking(1:n_take));
    % And top N profiles for CDF
    best_ra_cdf_fits(:,1:n_take,s) = ra_cdf_array_ft(:,fit_ranking(1:n_take));      
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

saveas(cdf_fig,[FigurePath 'cdf_fit_fig_data_only.png'])

for s = 1:length(master_struct)
%     plot(reactivation_time_axis/60,best_ra_cdf_fits(:,:,s),'Color',[cmap(s,:) .3])
    p(s+1) = plot(reactivation_time_axis/60,best_ra_cdf_fits(:,1,s),'-.','Color',cmap(s,:),'LineWidth',2);
end

legend(labelCell2{:},'Location','southeast','Color','w')
saveas(cdf_fig,[FigurePath 'cdf_fit_fig_all.png'])

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%% Make basic fit figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%

pon_true_profile = master_struct(1).sweepInfo.p_on_true;
time_vec = master_struct(1).sweepInfo.time_vec/60;
t_filter = master_struct(1).sweepInfo.t_filter;
mean_tf_vec = nanmean(master_struct(1).sweepInfo.tf_profile_array,2);

% generate an alternative filter that is more focused on the immediate
% vicinity of the event
t_bounds_fit = [-8 6];
t_filter_alt = time_vec <= t_bounds_fit(2) & time_vec >= t_bounds_fit(1);


close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PON
pon_fig = figure;
cmap = brewermap([],'Set2');
hold on
p = [];
% plot the true profile
p(1) = plot(time_vec,pon_true_profile,'Color','k','LineWidth',2.5);


% iterate through sim types
for s = 1:length(master_struct)
    % calculate errors for each predicted profile
    pon_fit_array = master_struct(s).sweepInfo.p_on_fit_array;
    diff2_array = (pon_fit_array-pon_true_profile).^2;
    master_struct(s).pon_objective_vec_alt = sum(diff2_array(t_filter_alt,:),1);
    
    % find best fit
    [~,mi_pon] = min(master_struct(s).pon_objective_vec_alt);
    master_struct(s).mi_pon = mi_pon;
    best_ids = best_ra_fit_ids(:,s);
    best_ids = best_ids(~isnan(best_ids));
    pon_fit_profiles = master_struct(s).sweepInfo.p_on_fit_array(:,best_ids);    
    master_struct(s).best_pon_fit_profile = pon_fit_profiles;
    
    % plot
%     plot(time_vec,pon_fit_profiles,'-','Color',[cmap(s,:) .2],'LineWidth',1)
    fit_spread = std(pon_fit_profiles,[],2)';
    fit_mean = mean(pon_fit_profiles,2)';
    fit_ub = fit_mean + fit_spread;
    fit_lb = fit_mean - fit_spread;
%     p(s+1) = plot(time_vec,fit_mean,'-.','Color',cmap(s,:),'LineWidth',1.5);
%     fill([time_vec' fliplr(time_vec')],[fit_ub fliplr(fit_lb)],cmap(s,:),'FaceAlpha',0.5,'EdgeColor','k')
    p(s+1) = errorbar(time_vec+rand(size(time_vec))*0.1,fit_mean,fit_spread,'Color',cmap(s,:),'LineWidth',1,'CapSize',3);
end

ylabel('instantaneous active fraction')

% yyaxis right
% plot(time_vec,mean_tf_vec,'-','Color',cmap(s+1,:),'LineWidth',2)

% ylim([0 2])
xlim(t_bounds_fit)

xlabel('minutes from light ON event');
% ylabel('average knirps concentration (au)')
% grid on
set(gca,'FontSize',14)
legend(p,labelCell2{:},'Location','southeast')
set(gca,'Color',[228,221,209]/255) 
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = cmap(s+1,:);
ax.XAxis(1).Color = 'k';

pon_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

% saveas(pon_fig,[FigurePath 'best_p_on_fits.png'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fluo full

fluo_full_true_profile = master_struct(1).sweepInfo.fluo_true;
p = [];
fluo_full_fig = figure;
cmap = brewermap([],'Set2');
hold on

% plot the true profile
p(1) = plot(time_vec,fluo_full_true_profile,'Color','k','LineWidth',2.5);


% iterate through sim types
for s = 1:length(master_struct)
    % calculate errors for each predicted profile
    fluo_full_fit_array = master_struct(s).sweepInfo.fluo_fit_array;
    diff2_array = (fluo_full_fit_array-fluo_full_true_profile).^2;
    master_struct(s).fluo_full_objective_vec_alt = sum(diff2_array(t_filter_alt,:),1);
    
    % find best fit
    [~,mi_fluo] = min(master_struct(s).fluo_full_objective_vec_alt);
    master_struct(s).mi_fluo= mi_fluo;
    
    best_ids = best_ra_fit_ids(:,s);
    best_ids = best_ids(~isnan(best_ids));
    
    fluo_full_fit_profiles = master_struct(s).sweepInfo.fluo_fit_array(:,best_ids);    
    master_struct(s).fluo_full_fit_profile2 = fluo_full_fit_profiles;
    
    fit_spread = std(fluo_full_fit_profiles,[],2)';
    fit_mean = mean(fluo_full_fit_profiles,2)';
    fit_ub = fit_mean + fit_spread;
    fit_lb = fit_mean - fit_spread;
    p(s+1) = errorbar(time_vec+rand(size(time_vec))*0.1,fit_mean,fit_spread,'Color',cmap(s,:),'LineWidth',1,'CapSize',3);
end

ylabel('average spot fluorescence (all loci) (au)')

% yyaxis right
% plot(time_vec,mean_tf_vec,'-','Color',cmap(s+1,:),'LineWidth',2)

% ylim([0 2])
xlim(t_bounds_fit)

xlabel('minutes from light ON event');
% ylabel('average knirps concentration (au)')
% grid on
set(gca,'FontSize',14)
legend(p,labelCell2{:},'Location','southwest')
set(gca,'Color',[228,221,209]/255) 
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = cmap(s+1,:);
ax.XAxis(1).Color = 'k';

fluo_full_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(fluo_full_fig,[FigurePath 'best_fluo_full_fits.png'])

%%
close all
fluo_raw_true_profile = master_struct(1).sweepInfo.fluo_true_raw;

fluo_raw_fig = figure;
cmap = brewermap([],'Set2');
hold on
p = [];
% plot the true profile
p(1) = plot(time_vec,fluo_raw_true_profile,'Color','k','LineWidth',2.5);


% iterate through sim types
for s = 1%1:length(master_struct)
    % calculate errors for each predicted profile
%     fluo_raw_fit_array = master_struct(s).sweepInfo.fluo_raw_fit_array;
%     diff2_array = (fluo_raw_fit_array-fluo_raw_true_profile).^2;
%     nan_flags = isnan(diff2_array(t_filter_alt,:));
%     master_struct(s).fluo_raw_objective_vec_alt = nanmean(diff2_array(t_filter_alt,:),1);
%     master_struct(s).fluo_raw_objective_vec_alt(mean(nan_flags)>0.25) = NaN;
    
    % find best fit
%     [~,mi_fluo_raw] = min(master_struct(s).fluo_raw_objective_vec_alt);
%     master_struct(s).mi_fluo_raw = mi_fluo_raw;
    
    best_ids = best_ra_fit_ids(:,s);
    best_ids = best_ids(~isnan(best_ids));
    
    fluo_raw_fit_profiles = master_struct(s).sweepInfo.fluo_raw_fit_array(:,best_ids);    
    master_struct(s).fluo_raw_fit_profiles = fluo_raw_fit_profiles;
    
    fit_spread = nanstd(fluo_raw_fit_profiles,[],2)';
    fit_mean = nanmean(fluo_raw_fit_profiles,2)';
    fit_ub = fit_mean + fit_spread;
    fit_lb = fit_mean - fit_spread;
    p(s+1) = errorbar(time_vec+rand(size(time_vec))*0.1,fit_mean,fit_spread,'Color',cmap(s,:),'LineWidth',1,'CapSize',3);
end

ylabel('average spot fluorescence (active loci only) (au)')
ylim([0 2.5e5]);

% yyaxis right
% plot(time_vec,mean_tf_vec,'-','Color',cmap(s+1,:),'LineWidth',2)

% ylim([0 2])
xlim(t_bounds_fit)

xlabel('minutes from light ON event');
% ylabel('average knirps concentration (au)')
% grid on
set(gca,'FontSize',14)
legend(p,labelCell2{:},'Location','southwest')
set(gca,'Color',[228,221,209]/255) 
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
% ax.YAxis(2).Color = cmap(s+1,:);
ax.XAxis(1).Color = 'k';

fluo_raw_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(fluo_raw_fig,[FigurePath 'best_fluo_raw_fits.png'])
%%
close all
figure;
hold on
for s = 1:length(master_struct)
    scatter(1./master_struct(s).pon_objective_vec_alt,1./master_struct(s).fluo_raw_objective_vec_alt)
end
legend(labelCell{:},'Location','northeast')
%%
% find best fit
[~,mi_pon] = min(sweepInfo.objective_val_p_on);
bestParamVals = sweepInfo.param_fit_array(mi_pon,:);
time_axis = sweepInfo.time_vec - sweepInfo.time_vec(pert_ind);

% examine goodness of fit
p_on_fit_fig = figure;
hold on
plot(time_axis,sweepInfo.p_on_fit_array(:,mi_pon)','--','LineWidth',2)
plot(time_axis,sweepInfo.p_on_true,'Color','k','LineWidth',2)

yyaxis right
plot(time_axis,nanmean(sweepInfo.tf_profile_array_true,2),'-g','LineWidth',1.5)
ylabel('knirps concentration (au)')
ax = gca;
ax.YAxis(2).Color = 'g';
grid on

set(gca,'Fontsize',14)
xlabel('minutes from perturbation')
ylabel('fraction active')

xlim([-11 11])

legend('best fit','experimental data','Location','southwest')
saveas(p_on_fit_fig,[FigurePath 'p_on_fit_' simTypeCell '.png'])

% now for fluorescence
[~,mi_fluo] = min(sweepInfo.objective_val_fluo);
fluo_fit_fig = figure;
hold on
plot(time_axis,sweepInfo.fluo_fit_array(:,mi_fluo)','--','LineWidth',2)
plot(time_axis,sweepInfo.fluo_true,'Color','k','LineWidth',2)
ylabel('average fluorescence (au)')

yyaxis right
plot(time_axis,nanmean(sweepInfo.tf_profile_array_true,2),'-g','LineWidth',1.5)
ylabel('knirps concentration (au)')

ax = gca;
ax.YAxis(2).Color = 'g';
grid on

set(gca,'Fontsize',14)
xlabel('minutes from perturbation')


xlim([-11 11])

legend('best fit','experimental data','Location','southwest')
saveas(fluo_fit_fig,[FigurePath 'fluo_fit_' simTypeCell '.png'])

%%
kd_fig = figure;
hold on
scatter(sweepInfo.param_fit_array(:,2),sweepInfo.objective_val_p_on)
set(gca,'Fontsize',14)
xlabel(['K_D (' rate_str ')'])
ylabel('residual error')
set(gca,'yscale','log');
grid on
saveas(kd_fig,[FigurePath 'KD_fit_scatter_' simTypeCell '.png'])

thresh_fig = figure;
hold on
scatter(sweepInfo.param_fit_array(:,3),sweepInfo.objective_val_p_on)
set(gca,'Fontsize',14)
xlabel('fluorescence detection threshold')
ylabel('residual error')
set(gca,'yscale','log');
grid on
saveas(thresh_fig,[FigurePath 'thresh_fit_scatter_' simTypeCell '.png'])

Kout_fig = figure;
hold on
scatter(sweepInfo.param_fit_array(:,4),sweepInfo.objective_val_p_on)
set(gca,'xscale','log');
set(gca,'yscale','log');
set(gca,'Fontsize',14)
xlabel('K_{out} (s^{-1})')
ylabel('residual error')
grid on
saveas(Kout_fig,[FigurePath 'kout_fit_scatter_' simTypeCell '.png'])

Kin_fig = figure;
hold on
scatter(sweepInfo.param_fit_array(:,5),sweepInfo.objective_val_p_on)
set(gca,'xscale','log');
set(gca,'yscale','log');
set(gca,'Fontsize',14)
xlabel('K_{in} (s^{-1})')
ylabel('residual error')
grid on
saveas(Kin_fig,[FigurePath 'kin_fit_scatter_' simTypeCell '.png'])

%%
ex_fig = figure;
plot(sweepInfo.gillespie{mi_pon}.t_ref,sweepInfo.gillespie{mi_pon}.fluo_ms2_array(:,10))
ylabel('simulated fluorescence (au)')
yyaxis right
plot(sweepInfo.gillespie{mi_pon}.t_ref,sweepInfo.tf_profile_array(:,10))
ylabel('measured single-nucles [knirps]')

set(gca,'Fontsize',14)
xlabel('minutes from perturbation')
saveas(ex_fig,[FigurePath 'io_trace_' simTypeCell '.png'])

io_fig = figure;
plot(sweepInfo.gillespie{mi_pon}.t_ref,sweepInfo.gillespie{mi_pon}.io_ref_in(:,1,10))
ylabel(rate_str)
yyaxis right
plot(sweepInfo.gillespie{mi_pon}.t_ref,sweepInfo.tf_profile_array(:,10))
ylabel('measured single-nucles [knirps]')

set(gca,'Fontsize',14)
xlabel('minutes from perturbation')
saveas(ex_fig,[FigurePath 'io_rate_' simTypeCell '.png'])

