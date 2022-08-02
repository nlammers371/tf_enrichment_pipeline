clear
close all

projectName = '2xDl-Ven_twiPEe-mCh';
nStates = 2;
% load data
if isfolder('C:\Users\nlamm\Dropbox (Personal)\')
    resultsRoot = 'C:\Users\nlamm\Dropbox (Personal)\ProcessedEnrichmentData';
else
   resultsRoot = 'S:\Nick\Dropbox (Personal)\ProcessedEnrichmentData';
end
load([resultsRoot filesep projectName filesep 'cpHMM_results' filesep 'compiledResults_w7_K' num2str(nStates) '_p1_ap1_t1_f2D.mat'])
figDir = [resultsRoot filesep projectName filesep 'cpHMM_results' filesep 'fig_K' num2str(nStates) filesep];
mkdir(figDir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Dorsal-dependent trends
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
cmap = brewermap([],'set2');

freq_trend = figure;
hold on

errorbar(compiledResults.protein_intensity_vec,compiledResults.freq_vec_mean,compiledResults.freq_vec_ste,'o','CapSize',0,'Color','k')
scatter(compiledResults.protein_intensity_vec,compiledResults.freq_vec_mean,75,'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k')
    
grid on

xlabel('[Dorsal] (au) ') 
ylabel('bursts per minute')

title('burst frequncy (k_{on}) vs. [Dorsal]')
set(gca,'Fontsize',14)

box on
saveas(freq_trend,[figDir, 'freq_vs_dl.png'])
saveas(freq_trend,[figDir, 'freq_vs_dl.pdf'])


dur_trend = figure;
hold on

errorbar(compiledResults.protein_intensity_vec,compiledResults.dur_vec_mean,compiledResults.dur_vec_ste,'o','CapSize',0,'Color','k')
scatter(compiledResults.protein_intensity_vec,compiledResults.dur_vec_mean,75,'MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k')
    
grid on

xlabel('[Dorsal] (au) ') 
ylabel('burst duration (minutes)')

title('burst duration (1/k_{off}) vs. [Dorsal]')
set(gca,'Fontsize',14)

box on
saveas(dur_trend,[figDir, 'dur_vs_dl.png'])
saveas(dur_trend,[figDir, 'dur_vs_dl.pdf'])


r_trend = figure;
hold on

errorbar(compiledResults.protein_intensity_vec,compiledResults.init_vec_mean,compiledResults.init_vec_ste,'o','CapSize',0,'Color','k')
scatter(compiledResults.protein_intensity_vec,compiledResults.init_vec_mean,75,'MarkerFaceColor',cmap(4,:),'MarkerEdgeColor','k')
    
grid on

xlabel('[Dorsal] (au) ') 
ylabel('au per minute')

title('burst amplitude (r) vs. [Dorsal]')
set(gca,'Fontsize',14)

box on
saveas(r_trend,[figDir, 'r_vs_dl.png'])
saveas(r_trend,[figDir, 'r_vs_dl.pdf'])


f_trend = figure;
hold on

errorbar(compiledResults.protein_intensity_vec,compiledResults.fluo_mean,compiledResults.fluo_ste,'o','CapSize',0,'Color','k')
scatter(compiledResults.protein_intensity_vec,compiledResults.fluo_mean,75,'MarkerFaceColor',cmap(5,:),'MarkerEdgeColor','k')
    
grid on

xlabel('[Dorsal] (au) ') 
ylabel('spot intensity (au)')

title('average transcription rate vs. [Dorsal]')
set(gca,'Fontsize',14)

box on
saveas(f_trend,[figDir, 'fluo_vs_dl.png'])
saveas(f_trend,[figDir, 'fluo_vs_dl.pdf'])


f_dur_trend = figure;
hold on

errorbar(compiledResults.protein_intensity_vec,compiledResults.fluo_mean,compiledResults.fluo_ste,'o','CapSize',0,'Color','k')
scatter(compiledResults.protein_intensity_vec,compiledResults.fluo_mean,75,'MarkerFaceColor',cmap(5,:),'MarkerEdgeColor','k')
xlabel('[Dorsal] (au) ') 
ylabel('spot intensity (au)')
ylim([30 55])
yyaxis right

errorbar(compiledResults.protein_intensity_vec,compiledResults.dur_vec_mean,compiledResults.dur_vec_ste,'o','CapSize',0,'Color','k')
scatter(compiledResults.protein_intensity_vec,compiledResults.dur_vec_mean,75,'MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k')

% grid on
ylabel('bursts duration (minutes)')

ax = gca;
ax.YAxis(1).Color = cmap(5,:);
ax.YAxis(2).Color = cmap(3,:);
set(gca,'Fontsize',14)
% ylim([0.10 .5])

box on
saveas(f_dur_trend,[figDir, 'fluo_dur_vs_dl.png'])
saveas(f_dur_trend,[figDir, 'fluo_dur_vs_dl.pdf'])


freq_trend = figure;
hold on

errorbar(compiledResults.fluo_mean,compiledResults.freq_vec_mean,compiledResults.freq_vec_ste,'o','CapSize',0,'Color','k')
scatter(compiledResults.fluo_mean,compiledResults.freq_vec_mean,75,'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k')
    
grid on

xlabel('mean spot intensity (au)') 
ylabel('bursts per minute')

title('burst frequncy (k_{on}) vs. spot intensity')
set(gca,'Fontsize',14)

box on
saveas(freq_trend,[figDir, 'freq_vs_fluo.png'])
saveas(freq_trend,[figDir, 'freq_vs_fluo.pdf'])


dur_trend = figure;
hold on

errorbar(compiledResults.fluo_mean,compiledResults.dur_vec_mean,compiledResults.dur_vec_ste,'o','CapSize',0,'Color','k')
scatter(compiledResults.fluo_mean,compiledResults.dur_vec_mean,75,'MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k')
    
grid on

xlabel('mean spot intensity (au)') 
ylabel('burst duration (minutes)')

title('burst duration (1/k_{off}) vs. spot intensity')
set(gca,'Fontsize',14)

box on
saveas(dur_trend,[figDir, 'dur_vs_fluo.png'])
saveas(dur_trend,[figDir, 'dur_vs_fluo.pdf'])


r_trend = figure;
hold on

errorbar(compiledResults.fluo_mean,compiledResults.init_vec_mean,compiledResults.init_vec_ste,'o','CapSize',0,'Color','k')
scatter(compiledResults.fluo_mean,compiledResults.init_vec_mean,75,'MarkerFaceColor',cmap(4,:),'MarkerEdgeColor','k')
    
grid on

xlabel('mean spot intensity (au)')  
ylabel('au per minute')

title('burst amplitude (r) vs. mean spot intensity')
set(gca,'Fontsize',14)
ylim([50 75])
box on
saveas(r_trend,[figDir, 'r_vs_fluo.png'])
saveas(r_trend,[figDir, 'r_vs_fluo.pdf'])

%%
pd_fluo = compiledResults.init_vec_mean.*(compiledResults.freq_vec_mean ./ ...
          (compiledResults.freq_vec_mean + 1./compiledResults.dur_vec_mean));
        
consistency_fig = figure;
scatter(compiledResults.fluo_mean,pd_fluo)
