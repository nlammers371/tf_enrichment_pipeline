clear 
close all

load("S:\Gabriella\Dropbox\ProcessedEnrichmentData\hbBAC-MS2-25C\cpHMM_results\compiledResults_w7_K3_p0_ap14_t8_f2D.mat")
FigurePath = 'S:\Gabriella\Dropbox\ProcessedEnrichmentData\hbBAC-MS2-25C\cpHMM_results\figures\w7_K3_p0_ap14_t8\';
mkdir(FigurePath);


time_group_vec = compiledResults.timeGroupVec;
time_group_index = unique(time_group_vec);
ap_group_vec = compiledResults.apGroupVec;
ap_group_index = unique(ap_group_vec);
ap_axis = compiledResults.apBins(1:end-1) + diff(compiledResults.apBins);
time_axis = [];
for t = 1:length(time_group_index)
   time_axis(t) = mean(compiledResults.timeBins{t});
end

%% burst duration

dur_fig = figure;
cmap = flipud(brewermap(length(time_axis),'Spectral'));
colormap(cmap);
hold on
for t = 1:length(time_group_index)
    time_filter = time_group_vec==t;
    ap_ids = ap_group_vec(time_filter);
    dur_vec_mean = compiledResults.dur_vec_mean(time_filter);
    dur_vec_ste = compiledResults.dur_vec_ste(time_filter);
    
    errorbar(ap_axis(ap_ids),dur_vec_mean,dur_vec_ste,'Color','k','Capsize',0)
    scatter(ap_axis(ap_ids),dur_vec_mean,'MarkerFaceColor',cmap(t,:),'MarkerEdgeColor','k')
end

h = colorbar;
% h.Ticks = 1:8;
h.TickLabels = time_axis/60;
ylabel(h,'time cohort (minutes into nc14)')
grid on
xlabel('AP position')
ylabel('burst duration (minutes)')
set(gca,'Fontsize',14)
xlim([20 47.5])

saveas(dur_fig,[FigurePath 'burst_duration.png'])


dur_fig = figure;
cmap = flipud(brewermap(length(time_axis),'Spectral'));
colormap(cmap);
hold on
for t = 1:length(time_group_index)
    time_filter = time_group_vec==t;
    ap_ids = ap_group_vec(time_filter);
    dur_vec_mean = compiledResults.dur_vec_mean(time_filter);
    dur_vec_ste = compiledResults.dur_vec_ste(time_filter);
    
    errorbar(ap_axis(ap_ids),dur_vec_mean,dur_vec_ste,'Color','k','Capsize',0)
    scatter(ap_axis(ap_ids),dur_vec_mean,'MarkerFaceColor',cmap(t,:),'MarkerEdgeColor','k')
end

h = colorbar;
% h.Ticks = 1:8;
h.TickLabels = time_axis/60;
ylabel(h,'time cohort (minutes into nc14)')
grid on
xlabel('AP position')
ylabel('burst duration (minutes)')
set(gca,'Fontsize',14)
xlim([20 47.5])

saveas(dur_fig,[FigurePath 'burst_duration.png'])

%% Frequency
freq_fig = figure;
cmap = flipud(brewermap(length(time_axis),'Spectral'));
colormap(cmap);
hold on
for t = 1:length(time_group_index)
    time_filter = time_group_vec==t;
    ap_ids = ap_group_vec(time_filter);
    freq_vec_mean = compiledResults.freq_vec_mean(time_filter);
    freq_vec_ste = compiledResults.freq_vec_ste(time_filter);
    
    errorbar(ap_axis(ap_ids),freq_vec_mean,freq_vec_ste,'Color','k','Capsize',0)
    scatter(ap_axis(ap_ids),freq_vec_mean,'MarkerFaceColor',cmap(t,:),'MarkerEdgeColor','k')
end

h = colorbar;
% h.Ticks = 1:8;
h.TickLabels = time_axis/60;
ylabel(h,'time cohort (minutes into nc14)')
grid on
xlabel('AP position')
ylabel('burst frequency (1/min)')
set(gca,'Fontsize',14)
xlim([20 47.5])

saveas(freq_fig,[FigurePath 'burst_frequency.png'])


%% Initiation
init_fig = figure;
colormap(cmap);
hold on
for t = 1:length(time_group_index)
    time_filter = time_group_vec==t;
    ap_ids = ap_group_vec(time_filter);
    init_vec_mean = compiledResults.init_vec_mean(time_filter);
    init_vec_ste = compiledResults.init_vec_ste(time_filter);
    
    errorbar(ap_axis(ap_ids),init_vec_mean,init_vec_ste,'Color','k','Capsize',0)
    scatter(ap_axis(ap_ids),init_vec_mean,'MarkerFaceColor',cmap(t,:),'MarkerEdgeColor','k')
end

h = colorbar;
% h.Ticks = 1:8;
h.TickLabels = time_axis/60;
ylabel(h,'time cohort (minutes into nc14)')
grid on
xlabel('AP position')
ylabel('initiation rate (au/min)')
set(gca,'Fontsize',14)
xlim([20 47.5])

saveas(init_fig,[FigurePath 'burst_initiation.png'])