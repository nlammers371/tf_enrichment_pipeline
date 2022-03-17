clear all
close all
clc


%% initialization

projectName = 'optokni_eve4+6_WT'; 

liveProject = LiveEnrichmentProject(projectName);
resultsRoot = [liveProject.dataPath filesep];    

FigurePath = [liveProject.figurePath 'bursting' filesep];
mkdir(FigurePath)

% Define some other colors  
yw = [234 194 100]/255; % yellow
bl = [115 143 193]/255; % blue
gr = [191 213 151]/255; % green

%% Figure: Plot k_on, k_off vs ap

% load inference results
%WT_bursting = load('P:\Jake\Pipeline\tf_enrichment_pipeline\opto_knirps_project\bursting_data\compiledResults_w7_K2_p0_ap7_t1_f2D_temp_WT.mat');
%WT_bursting = load('P:\Jake\Pipeline\tf_enrichment_pipeline\opto_knirps_project\bursting_data\compiledResults_w7_K2_p0_ap9_t1_f2D_WT.mat');
WT_bursting = load('P:\Jake\Pipeline\tf_enrichment_pipeline\opto_knirps_project\bursting_data\compiledResults_w7_K2_p0_ap11_t1_f2D_v1_4500.mat');


burst_axis = (WT_bursting.compiledResults.apBins(1:end-1)+ WT_bursting.compiledResults.apBins(2:end))/2;
ap_axis = burst_axis*100;
% data from WT
burst_freq_WT = WT_bursting.compiledResults.freq_vec_mean;
burst_freq_ste_WT = WT_bursting.compiledResults.freq_vec_ste;
burst_dur_WT = WT_bursting.compiledResults.dur_vec_mean;
burst_dur_ste_WT = WT_bursting.compiledResults.dur_vec_ste;
burst_rate_WT= WT_bursting.compiledResults.init_vec_mean*1e-5;
burst_rate_ste_WT = WT_bursting.compiledResults.init_vec_ste*1e-5;


burst_dur_center = burst_dur_WT(4);
burst_rate_center = burst_rate_WT(4);
burst_freq_center = burst_freq_WT(4);

burst_k_off_center = 1/burst_dur_center;
burst_k_on_center = burst_freq_center;

burst_k_off_mean_WT = 1./WT_bursting.compiledResults.dur_vec_mean;
burst_k_off_ste_WT = 1./(WT_bursting.compiledResults.dur_vec_mean.^2).*WT_bursting.compiledResults.dur_vec_ste;
burst_k_on_mean_WT = WT_bursting.compiledResults.freq_vec_mean;
burst_k_on_ste_WT = WT_bursting.compiledResults.freq_vec_ste;

% k_off vs ap
burst_k_off_fig = figure;
hold on

set(gca,'FontSize',14)


errorbar(burst_axis*100,burst_k_off_mean_WT/burst_k_off_center,burst_k_off_ste_WT/burst_k_off_center,'Color','k','CapSize',0)
plot(burst_axis*100,burst_k_off_mean_WT/burst_k_off_center,'-k')
scatter(burst_axis*100,burst_k_off_mean_WT/burst_k_off_center,50,'MarkerFaceColor',yw,'MarkerEdgeColor','k')
%set(gca,'YColor',bl);
ylabel(['k_{off} (relative to WT center)'])
ylim([0.25 1.75])


xlabel('AP position (% embryo length)');
xlim([ap_axis(1) ap_axis(end)])
%mRNA_fig2.InvertHardcopy = 'off';
set(gcf,'color','w'); 
pbaspect([3 2 1])
%legend('','','','','','')


% k_on vs ap
burst_k_on_fig = figure;
hold on

set(gca,'FontSize',14)

errorbar(burst_axis*100,burst_k_on_mean_WT/burst_k_on_center,burst_k_on_ste_WT/burst_k_on_center,'Color','k','CapSize',0)
plot(burst_axis*100,burst_k_on_mean_WT/burst_k_on_center,'-k')
scatter(burst_axis*100,burst_k_on_mean_WT/burst_k_on_center,50,'MarkerFaceColor',yw,'MarkerEdgeColor','k')

%set(gca,'YColor',bl);
ylabel(['k_{on} (relative to WT center)'])
ylim([0.25 1.75])


xlabel('AP position (% embryo length)');
xlim([ap_axis(1) ap_axis(end)])
%mRNA_fig2.InvertHardcopy = 'off';
set(gcf,'color','w'); 
pbaspect([3 2 1])
legend('','','','','','')


% loading rate vs ap
burst_loading_rate_ap_fig = figure;
hold on

set(gca,'FontSize',14)

errorbar(burst_axis*100,burst_rate_WT/burst_rate_center,burst_rate_ste_WT/burst_rate_center,'Color','k','CapSize',0)
plot(burst_axis*100,burst_rate_WT/burst_rate_center,'-k')
scatter(burst_axis*100,burst_rate_WT/burst_rate_center,50,'MarkerFaceColor',yw,'MarkerEdgeColor','k')
%set(gca,'YColor',gr);
ylabel(['mRNA loading rate (relative to WT center)'])
ylim([0 2])


xlabel('AP position (% embryo length)');
xlim([ap_axis(1) ap_axis(end)])
%mRNA_fig2.InvertHardcopy = 'off';
set(gcf,'color','w'); 
pbaspect([3 2 1])


saveas(burst_k_on_fig,[FigurePath 'figure_k_on_vs_ap_WT.pdf'])
saveas(burst_k_off_fig,[FigurePath 'figure_k_off_vs_ap_WT.pdf'])
saveas(burst_loading_rate_ap_fig,[FigurePath 'figure_loading_rate_vs_ap.pdf'])

%% plot bursting rates + loading rate

% loading rate vs ap
burst_loading_rate_ap_fig = figure;
hold on

set(gca,'FontSize',14)
%set(gca, 'YScale', 'log')

errorbar(burst_axis*100,burst_rate_WT,burst_rate_ste_WT,'Color','k','CapSize',0)
plot(burst_axis*100,burst_rate_WT,'-k')
scatter(burst_axis*100,burst_rate_WT,50,'MarkerFaceColor',bl,'MarkerEdgeColor','k')
%set(gca,'YColor',gr);
ylabel(['mRNA loading rate (AU/min)'])
ylim([0 4])

xlabel('AP position (% embryo length)');
xlim([ap_axis(1) ap_axis(end)])
%mRNA_fig2.InvertHardcopy = 'off';
set(gcf,'color','w'); 
pbaspect([2 1 1])



% transition rate vs ap
burst_rate_fig = figure;
hold on

set(gca,'FontSize',14)
set(gca, 'YScale', 'log')

errorbar(burst_axis*100,burst_k_on_mean_WT,burst_k_on_ste_WT,'Color','k','CapSize',0)
plot(burst_axis*100,burst_k_on_mean_WT,'-k')
scatter(burst_axis*100,burst_k_on_mean_WT,50,'MarkerFaceColor',bl,'MarkerEdgeColor','k')

errorbar(burst_axis*100,burst_k_off_mean_WT,burst_k_off_ste_WT,'Color','k','CapSize',0)
plot(burst_axis*100,burst_k_off_mean_WT,'-k')
scatter(burst_axis*100,burst_k_off_mean_WT,50,'MarkerFaceColor',bl,'MarkerEdgeColor','k')

%set(gca,'YColor',bl);
ylabel(['transition rates (1/min)'])
ylim([0 10])


xlabel('AP position (% embryo length)');
xlim([ap_axis(1) ap_axis(end)])
%mRNA_fig2.InvertHardcopy = 'off';
set(gcf,'color','w'); 
pbaspect([2 1 1])
%legend('','','','','','')


%saveas(burst_rate_fig,[FigurePath 'figure_burst_rate_vs_ap.pdf'])
%saveas(burst_loading_rate_ap_fig,[FigurePath 'figure_loading_rate_vs_ap.pdf'])


%% Figure: Plot k_on, k_off vs time

% load inference results
WT_temp_bursting = load('P:\Jake\Pipeline\tf_enrichment_pipeline\opto_knirps_project\bursting_data\compiledResults_w7_K2_p0_ap1_t11_f2D_WT.mat');

burst_axis = linspace(7.5,32.5,11);

% data from WT
burst_freq_WT = WT_temp_bursting.compiledResults.freq_vec_mean;
burst_freq_ste_WT = WT_temp_bursting.compiledResults.freq_vec_ste;
burst_dur_WT = WT_temp_bursting.compiledResults.dur_vec_mean;
burst_dur_ste_WT = WT_temp_bursting.compiledResults.dur_vec_ste;
burst_rate_WT= WT_temp_bursting.compiledResults.init_vec_mean*1e-5;
burst_rate_ste_WT = WT_temp_bursting.compiledResults.init_vec_ste*1e-5;


burst_dur_first = burst_dur_WT(1);
burst_rate_first = burst_rate_WT(1);
burst_freq_first = burst_freq_WT(1);

burst_k_off_center = 1/burst_dur_first;
burst_k_on_center = burst_freq_first;
burst_rate_center = burst_rate_first;

burst_k_off_mean_WT = 1./WT_temp_bursting.compiledResults.dur_vec_mean;
burst_k_off_ste_WT = 1./(WT_temp_bursting.compiledResults.dur_vec_mean.^2).*WT_temp_bursting.compiledResults.dur_vec_ste;
burst_k_on_mean_WT = WT_temp_bursting.compiledResults.freq_vec_mean;
burst_k_on_ste_WT = WT_temp_bursting.compiledResults.freq_vec_ste;


% k_off vs time
burst_k_off_fig = figure;
hold on

set(gca,'FontSize',14)


errorbar(burst_axis,burst_k_off_mean_WT/burst_k_off_center,burst_k_off_ste_WT/burst_k_off_center,'Color','k','CapSize',0)
plot(burst_axis,burst_k_off_mean_WT/burst_k_off_center,'-k')
scatter(burst_axis,burst_k_off_mean_WT/burst_k_off_center,50,'MarkerFaceColor',yw,'MarkerEdgeColor','k')
%set(gca,'YColor',bl);
ylabel(['k_{off}'])
ylim([0.5 1.5])


xlabel('time (min) into nc14');
xlim([5 30])
%mRNA_fig2.InvertHardcopy = 'off';
set(gcf,'color','w'); 
pbaspect([2 1 1])


% k_on vs time
burst_k_on_fig = figure;
hold on

set(gca,'FontSize',14)

errorbar(burst_axis,burst_k_on_mean_WT/burst_k_on_center,burst_k_on_ste_WT/burst_k_on_center,'Color','k','CapSize',0)
plot(burst_axis,burst_k_on_mean_WT/burst_k_on_center,'-k')
scatter(burst_axis,burst_k_on_mean_WT/burst_k_on_center,50,'MarkerFaceColor',yw,'MarkerEdgeColor','k')

%set(gca,'YColor',bl);
ylabel(['k_{on}'])
%ylim([0 2])
ylim([0.5 1.25])


xlabel('time (min) into nc14');
xlim([5 30])
%mRNA_fig2.InvertHardcopy = 'off';
set(gcf,'color','w'); 
pbaspect([2 1 1])
legend('','','','','','')


% loading rate vs time
burst_loading_rate_ap_fig = figure;
hold on

set(gca,'FontSize',14)

errorbar(burst_axis,burst_rate_WT/burst_rate_center,burst_rate_ste_WT/burst_rate_center,'Color','k','CapSize',0)
plot(burst_axis,burst_rate_WT/burst_rate_center,'-k')
scatter(burst_axis,burst_rate_WT/burst_rate_center,50,'MarkerFaceColor',yw,'MarkerEdgeColor','k')
%set(gca,'YColor',gr);
ylabel(['mRNA loading rate'])
ylim([0 2])


xlabel('time (min) into nc14');
xlim([5 30])
%mRNA_fig2.InvertHardcopy = 'off';
set(gcf,'color','w'); 
pbaspect([2 1 1])


%saveas(burst_k_on_fig,[FigurePath 'figure_k_on_vs_ap_WT_CONST.pdf'])
%saveas(burst_k_off_fig,[FigurePath 'figure_k_off_vs_ap_WT_CONST.pdf'])
%saveas(burst_loading_rate_ap_fig,[FigurePath 'figure_loading_rate_vs_ap.pdf'])

%% plot temporal trends on the same figure

% k_on vs time
burst_k_on_off_fig = figure;
hold on

set(gca,'FontSize',14)
%set(gca, 'YScale', 'log')

errorbar(burst_axis,burst_k_on_mean_WT,burst_k_on_ste_WT,'Color','k','CapSize',0)
plot(burst_axis,burst_k_on_mean_WT,'-k')
scatter(burst_axis,burst_k_on_mean_WT,50,'MarkerFaceColor',yw,'MarkerEdgeColor','k')

errorbar(burst_axis,burst_k_off_mean_WT,burst_k_on_ste_WT,'Color','k','CapSize',0)
plot(burst_axis,burst_k_off_mean_WT,'-k')
scatter(burst_axis,burst_k_off_mean_WT,50,'MarkerFaceColor',yw,'MarkerEdgeColor','k')

%set(gca,'YColor',bl);
ylabel(['transition rates (1/min)'])
ylim([0 7])


xlabel('time (min) into nc14');
xlim([5 30])
%mRNA_fig2.InvertHardcopy = 'off';
set(gcf,'color','w'); 
pbaspect([2 1 1])

%saveas(burst_k_on_off_fig,[FigurePath 'figure_k_on_off_vs_time_WT_CONST.pdf'])
