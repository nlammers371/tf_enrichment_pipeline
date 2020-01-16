% Script to estimate elongation time using average auto-correlation of MS2
% traces
clear
close all
addpath('../utilities')
% define core ID variables
project = 'Dl-Ven_snaBAC-mCh_v3';
% project = 'Dl-Ven_hbP2P-mCh';
DropboxFolder =  'E:\Nick\LivemRNA\Dropbox (Personal)\';
[~, DataPath, FigureRoot] =   header_function(DropboxFolder, project);

FigPath = [FigureRoot '\elongation_rate\' project '\'];
mkdir(FigPath)
w = 7;
K = 3;
n_lags = 20;
n_boots = 100;
load([DataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'],'hmm_input_output')
% compile traces into single array
time_grid = unique([hmm_input_output.time]);
trace_array = zeros(numel(time_grid),numel(hmm_input_output));
for i = 1:numel(hmm_input_output)
    time = hmm_input_output(i).time;
    fluo = hmm_input_output(i).fluo;
    trace_array(ismember(time_grid,time),i) = fluo;
end

[wt_autocorr, a_boot_errors, wt_dd, dd_boot_errors] = ...
    weighted_autocorrelation(trace_array, n_lags, 1,n_boots,ones(size(trace_array)));
%%% make figure
x_axis = (0:n_lags)*20/60;
dd_elongation_fig = figure;
hold on
e = errorbar(x_axis(2:end-1),wt_dd,dd_boot_errors,'color','black','LineWidth',1.5);
% scatter(x_axis,wt_dd,20,'filled');
e.CapSize = 0;
grid on 
xlabel('time delay (min)')
ylabel('second derivative of ACF')
set(gca,'Fontsize',12)
saveas(dd_elongation_fig, [FigPath 'acf_dd_plot.tif'])
saveas(dd_elongation_fig, [FigPath 'acf_dd_plot.pdf'])

elongation_fig = figure;
hold on
e = errorbar(x_axis,wt_autocorr, a_boot_errors,'color','black','LineWidth',1.5);
% scatter(x_axis,wt_dd,20,'filled');
e.CapSize = 0;
grid on 
xlabel('time delay (min)')
ylabel('ACF')
set(gca,'Fontsize',12)
saveas(elongation_fig, [FigPath 'acf_plot.tif'])
saveas(elongation_fig, [FigPath 'acf_plot.pdf'])