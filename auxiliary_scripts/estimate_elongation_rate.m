% Script to estimate elongation time using average auto-correlation of MS2
% traces
clear
close all
addpath('../utilities')
% define core ID variables
project = 'Dl-Ven_snaBAC-mCh_v3';
% project = 'Dl-Ven_hbP2P-mCh_v2';
DropboxFolder =  'E:\Nick\LivemRNA\Dropbox (Personal)\';
[~, DataPath, FigureRoot] =   header_function(DropboxFolder, project);

FigPath = [FigureRoot '\elongation_rate\' project '\'];
mkdir(FigPath)
% w = 7;
% K = 3;
n_lags = 20;
n_boots = 100;
load([DataPath 'nucleus_struct.mat'])
%%
qc_indices = find([nucleus_struct.qc_flag]==1);
%%
% load([DataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'],'hmm_input_output')
% compile traces into single array
time_grid = unique([nucleus_struct.time_interp]);
time_grid = time_grid(~isnan(time_grid));
trace_array = zeros(numel(time_grid),numel(qc_indices));
for i = 1:numel(qc_indices)
    qc = qc_indices(i);
    time = nucleus_struct(qc).time_interp;
    fluo = nucleus_struct(qc).fluo_interp;
    trace_array(ismember(time_grid,time),i) = fluo;
end
%%
[wt_autocorr, a_boot_errors, wt_dd, dd_boot_errors, wt_ddd, ddd_boot_errors] = ...
    weighted_autocorrelation(trace_array, n_lags, 1,n_boots,ones(size(trace_array)));
%% make figure
close all
x_axis1 = (0:n_lags);
x_axis2 = 1:n_lags-1;
dd_elongation_fig = figure;
hold on
e = errorbar(x_axis2,wt_dd,dd_boot_errors,'-o','color','black','LineWidth',1.5);
% scatter(x_axis,wt_dd,20,'filled');
e.CapSize = 0;
grid on 
xlabel('time delay (time steps)')
ylabel('second derivative of ACF')
set(gca,'Fontsize',12)
saveas(dd_elongation_fig, [FigPath 'acf_dd_plot.tif'])
saveas(dd_elongation_fig, [FigPath 'acf_dd_plot.pdf'])

elongation_fig = figure;
hold on
e = errorbar(x_axis1,wt_autocorr, a_boot_errors,'-o','color','black','LineWidth',1.5);
e.CapSize = 0;
grid on 
xlabel('time delay (time steps)')
ylabel('ACF')
set(gca,'Fontsize',12)
saveas(elongation_fig, [FigPath 'acf_plot.tif'])
saveas(elongation_fig, [FigPath 'acf_plot.pdf'])



ddd_elongation_fig = figure;
hold on
e = errorbar(x_axis2(1:end-1),wt_ddd,ddd_boot_errors,'-o','color','black','LineWidth',1.5);
% scatter(x_axis,wt_dd,20,'filled');
e.CapSize = 0;
grid on 
xlabel('time delay (time steps)')
ylabel('third derivative of ACF')
set(gca,'Fontsize',12)
ylim([-.01 .01])
saveas(ddd_elongation_fig, [FigPath 'acf_ddd_plot.tif'])
saveas(ddd_elongation_fig, [FigPath 'acf_ddd_plot.pdf'])
