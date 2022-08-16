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
master_struct = struct;
for s = 1:length(simTypeCell)
    % load and store sweep result sets
    simType = simTypeCell{s};
    load([resultsRoot 'sweepInfo_' simType '.mat'],'sweepInfo');
    master_struct(s).sweepInfo = sweepInfo;
    clear sweepInfo
end    

% mak figure directory  
try
    FigurePath = [liveProject.figurePath 'io_sim_results' filesep];
catch
    FigurePath = 'C:\Users\nlamm\Dropbox (Personal)\LocalEnrichmentFigures\PipelineOutput\optokni_eve4+6_ON';
end
mkdir(FigurePath)

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


labelCell2 = [{'data'} labelCell];
%%
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PON
pon_fig = figure;
cmap = brewermap([],'Set2');
hold on

% plot the true profile
plot(time_vec,pon_true_profile,'Color','k','LineWidth',2.5)


% iterate through sim types
for s = 1:length(master_struct)
    % calculate errors for each predicted profile
    pon_fit_array = master_struct(s).sweepInfo.p_on_fit_array;
    diff2_array = (pon_fit_array-pon_true_profile).^2;
    master_struct(s).pon_objective_vec_alt = sum(diff2_array(t_filter_alt,:),1);
    
    % find best fit
    [~,mi_pon] = min(master_struct(s).pon_objective_vec_alt);
    master_struct(s).mi_pon = mi_pon;
    pon_fit_profile = master_struct(s).sweepInfo.p_on_fit_array(:,master_struct(s).mi_pon);    
    master_struct(s).best_pon_fit_profile = pon_fit_profile;
    
    % plot
    plot(time_vec,pon_fit_profile,'-.','Color',cmap(s,:),'LineWidth',2)
end

ylabel('instantaneous active fraction')

yyaxis right
plot(time_vec,mean_tf_vec,'-','Color',cmap(s+1,:),'LineWidth',2)

% ylim([0 2])
xlim(t_bounds_fit)

xlabel('minutes from light ON event');
ylabel('average knirps concentration (au)')
% grid on
set(gca,'FontSize',14)
legend(labelCell2{:},'Location','southwest')
set(gca,'Color',[228,221,209]/255) 
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = cmap(s+1,:);
ax.XAxis(1).Color = 'k';

pon_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(pon_fig,[FigurePath 'best_p_on_fits.png'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fluo full

fluo_full_true_profile = master_struct(1).sweepInfo.fluo_true;

fluo_full_fig = figure;
cmap = brewermap([],'Set2');
hold on

% plot the true profile
plot(time_vec,fluo_full_true_profile,'Color','k','LineWidth',2.5)


% iterate through sim types
for s = 1:length(master_struct)
    % calculate errors for each predicted profile
    fluo_full_fit_array = master_struct(s).sweepInfo.fluo_fit_array;
    diff2_array = (fluo_full_fit_array-fluo_full_true_profile).^2;
    master_struct(s).fluo_full_objective_vec_alt = sum(diff2_array(t_filter_alt,:),1);
    
    % find best fit
    [~,mi_fluo] = min(master_struct(s).fluo_full_objective_vec_alt);
    master_struct(s).mi_fluo= mi_fluo;
    fluo_full_fit_profile = master_struct(s).sweepInfo.fluo_fit_array(:,master_struct(s).mi_pon);    
    master_struct(s).fluo_full_fit_profile = fluo_full_fit_profile;
    
    % plot
    plot(time_vec,fluo_full_fit_profile,'-.','Color',cmap(s,:),'LineWidth',2)
end

ylabel('average spot fluorescence (all loci) (au)')

yyaxis right
plot(time_vec,mean_tf_vec,'-','Color',cmap(s+1,:),'LineWidth',2)

% ylim([0 2])
xlim(t_bounds_fit)

xlabel('minutes from light ON event');
ylabel('average knirps concentration (au)')
% grid on
set(gca,'FontSize',14)
legend(labelCell2{:},'Location','southwest')
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

fluo_raw_true_profile = master_struct(1).sweepInfo.fluo_true_raw;

fluo_raw_fig = figure;
cmap = brewermap([],'Set2');
hold on

% plot the true profile
plot(time_vec,fluo_raw_true_profile,'Color','k','LineWidth',2.5)


% iterate through sim types
for s = 1:length(master_struct)
    % calculate errors for each predicted profile
    fluo_raw_fit_array = master_struct(s).sweepInfo.fluo_raw_fit_array;
    diff2_array = (fluo_raw_fit_array-fluo_raw_true_profile).^2;
    nan_flags = isnan(diff2_array(t_filter_alt,:));
    master_struct(s).fluo_raw_objective_vec_alt = nanmean(diff2_array(t_filter_alt,:),1);
    master_struct(s).fluo_raw_objective_vec_alt(mean(nan_flags)>0.25) = NaN;
    
    % find best fit
    [~,mi_fluo_raw] = min(master_struct(s).fluo_raw_objective_vec_alt);
    master_struct(s).mi_fluo_raw = mi_fluo_raw;
    fluo_raw_fit_profile = master_struct(s).sweepInfo.fluo_raw_fit_array(:,master_struct(s).mi_fluo_raw);    
    master_struct(s).fluo_raw_fit_profile = fluo_raw_fit_profile;
    
    % plot
    plot(time_vec,fluo_raw_fit_profile,'-.','Color',cmap(s,:),'LineWidth',2)
end

ylabel('average spot fluorescence (active loci only) (au)')
ylim([0 2.5e5]);

yyaxis right
plot(time_vec,mean_tf_vec,'-','Color',cmap(s+1,:),'LineWidth',2)

% ylim([0 2])
xlim(t_bounds_fit)

xlabel('minutes from light ON event');
ylabel('average knirps concentration (au)')
% grid on
set(gca,'FontSize',14)
legend(labelCell2{:},'Location','southwest')
set(gca,'Color',[228,221,209]/255) 
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = cmap(s+1,:);
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

