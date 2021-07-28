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

simTypeCell = {'koff_only_2','kon_only_2'};
if contains(simTypeCell,'in_only')
  rate_str = 'k_{in}';
elseif contains(simTypeCell,'out_only')
  rate_str = 'k_{out}';
end
  
load([resultsRoot 'sweepInfo_' simTypeCell '.mat'])
load([resultsRoot 'gillespie_' simTypeCell '.mat'],'gillespie_struct');

try
    FigurePath = [liveProject.figurePath 'io_sim_results' filesep];
catch
    FigurePath = 'C:\Users\nlamm\Dropbox (Personal)\LocalEnrichmentFigures\PipelineOutput\optokni_eve4+6_WT';
end
mkdir(FigurePath)

close all
paramList = sweepInfo.paramList;
pert_ind = 34;

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

