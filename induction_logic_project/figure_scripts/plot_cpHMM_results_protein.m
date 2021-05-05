% script to plot results from cpHMM inference
clear
close all
addpath(genpath('utilities'))

projectName = '2xDl-Ven_snaBAC-mCh';
inferenceName = 'w7_K3_p1_ap1_t1_f2D';
venusFlag = contains(projectName,'Ven');
proteinString = 'Dl';
proteinStringLong = 'Dorsal';

if contains(projectName,'Bcd')
    proteinString = 'Bcd';
    proteinStringLong = 'Bicoid';
end

MarkerSize = 50;
blue = [115 143 193]/256;
purple = [171 133 172]/256;
red = [213 108 85]/256;
yellow = [234 195 100]/256;

% set axes
dur_lims = [0 3.5];
freq_lims = [0 3];
init_lims = [0 18]*1e4;

% get path to results
liveProject = LiveEnrichmentProject(projectName);
resultsRoot = [liveProject.dataPath filesep];
resultsDir = [resultsRoot 'cpHMM_results' filesep];

% make figure directory
figureDir = [resultsDir 'figures' filesep];
mkdir(figureDir);

% get list of projects
inferenceFilePath = [resultsDir 'compiledResults_' inferenceName '.mat'];            

% load data
load(inferenceFilePath);

% get index of x axis
load([resultsRoot 'proteinSamplingInfo.mat'])
load("C:\Users\nlamm\Dropbox (Personal)\ProcessedEnrichmentData\absolute_calibration\calibration_info.mat")

% calculate the volume of the Gaussian integration kernel
samplingXY_um = proteinSamplingInfo.xy_sigma_um;
samplingZ_um = proteinSamplingInfo.z_sigma_um;
samplingKernelVolume_um = ((2*pi)^1.5*samplingXY_um^2*samplingZ_um);
samplingKernelVolume_px = samplingKernelVolume_um./calibration_info.VoxelSize;

if venusFlag
    calFactorAbs = calibration_info.venus_au_per_molecule;
    calFactorNM = calibration_info.venus_au_per_nM;
else
    calFactorAbs = calibration_info.gfp_au_per_molecule;
    calFactorNM = calibration_info.gfp_au_per_nM;
end

protein_axis = compiledResults.protein_intensity_vec/calFactorNM;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make figures
close all
pt_axis_long = linspace(0, protein_axis(end)*1.5);


r_mdl = fitlm(protein_axis,compiledResults.init_vec_mean);
r_pd = r_mdl.Coefficients.Estimate(2)*pt_axis_long + r_mdl.Coefficients.Estimate(1);

r_trend = figure;

hold on
plot(pt_axis_long,r_pd,'--k','LineWidth',1.5)
e = errorbar(protein_axis,compiledResults.init_vec_mean,compiledResults.init_vec_ste, 'o','Color','black','LineWidth',1);          
e.CapSize = 0;
scatter(protein_axis,compiledResults.init_vec_mean,MarkerSize,'o','MarkerFaceColor',yellow,'MarkerEdgeColor','black');   

% grid on

xlabel(['nuclear [' proteinString '] (nM)'])
ylabel('burst amplitude (au/min)')
grid on
% title(['Burst Amplitude (r): ' projectName])
set(gca,'Fontsize',14)
ylim([50 75]);
xlim([protein_axis(1)-2.5 protein_axis(end)+2.5])
%         StandardFigure([],gca)
ax = gca;
ax.YColor = 'black';%cmap1(2,:);
set(gca,'Fontsize',14)

set(gca,    'Box','off',...
            'Color',[228,221,209]/255,...            
            'TickLength',[0.02,0.05])    
r_trend.Color = 'white';        
r_trend.InvertHardcopy = 'off';

saveas(r_trend,[figureDir, 'burst_amp_vs_' proteinString '.tif'])
saveas(r_trend,[figureDir, 'burst_amp_vs_' proteinString '.pdf'])

%% Burst frequency
f_mdl = fitlm(protein_axis,compiledResults.freq_vec_mean);
f_pd = f_mdl.Coefficients.Estimate(2)*pt_axis_long + f_mdl.Coefficients.Estimate(1);

f_trend = figure;

hold on
plot(pt_axis_long,f_pd,'--k','LineWidth',1.5)
e = errorbar(protein_axis,compiledResults.freq_vec_mean,compiledResults.freq_vec_ste, 'o','Color','black','LineWidth',1);          
e.CapSize = 0;
scatter(protein_axis,compiledResults.freq_vec_mean,MarkerSize,'o','MarkerFaceColor',purple,'MarkerEdgeColor','black');   

% grid on

xlabel(['nuclear [' proteinString '] (nM)'])
ylabel('burst frequency (1/min)')
grid on
% title(['Burst Amplitude (r): ' projectName])
set(gca,'Fontsize',14)
% ylim([50 75]);
xlim([protein_axis(1)-2.5 protein_axis(end)+2.5])
%         StandardFigure([],gca)
ax = gca;
ax.YColor = 'black';%cmap1(2,:);
set(gca,'Fontsize',14)

set(gca,    'Box','off',...
            'Color',[228,221,209]/255,...            
            'TickLength',[0.02,0.05])    
f_trend.Color = 'white';        
f_trend.InvertHardcopy = 'off';

saveas(f_trend,[figureDir, 'burst_freq_vs_' proteinString '.tif'])
saveas(f_trend,[figureDir, 'burst_freq_vs_' proteinString '.pdf'])

%% Burst duration
d_mdl = fitlm(protein_axis,compiledResults.dur_vec_mean);
d_pd = d_mdl.Coefficients.Estimate(2)*pt_axis_long + d_mdl.Coefficients.Estimate(1);

d_trend = figure;

hold on
plot(pt_axis_long,d_pd,'--k','LineWidth',1.5)
e = errorbar(protein_axis,compiledResults.dur_vec_mean,compiledResults.dur_vec_ste, 'o','Color','black','LineWidth',1);          
e.CapSize = 0;
scatter(protein_axis,compiledResults.dur_vec_mean,MarkerSize,'o','MarkerFaceColor',blue,'MarkerEdgeColor','black');   

% grid on

xlabel(['nuclear [' proteinString '] (nM)'])
ylabel('burst duration (min)')
grid on
% title(['Burst Amplitude (r): ' projectName])
set(gca,'Fontsize',14)
ylim([.5 2.7]);
xlim([protein_axis(1)-2.5 protein_axis(end)+2.5])
%         StandardFigure([],gca)
ax = gca;
ax.YColor = 'black';%cmap1(2,:);
set(gca,'Fontsize',14)

set(gca,    'Box','off',...
            'Color',[228,221,209]/255,...            
            'TickLength',[0.02,0.05])    
d_trend.Color = 'white';        
d_trend.InvertHardcopy = 'off';

saveas(d_trend,[figureDir, 'burst_dur_vs_' proteinString '.tif'])
saveas(d_trend,[figureDir, 'dur_vs_' proteinString '.pdf'])