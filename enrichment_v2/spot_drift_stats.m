% script to make x-y-z drift histograms
clear
close all

projectName = 'Bcd-GFP_hbMS2-mCh_AiryscanTest_';

liveProject = LiveEnrichmentProject(projectName);
resultsRoot = [liveProject.dataPath filesep];
FigurePath = liveProject.figurePath;

% load data
load([resultsRoot 'spot_struct.mat'])
load([resultsRoot 'proteinSamplingInfo.mat'])

%% calculate frame-over-frame differences
close all

x_drift = [];
y_drift = [];
z_drift = [];
t_drift = [];

for i = 1:length(spot_struct)
    dt = diff(spot_struct(i).time);
    
    dx = diff(spot_struct(i).xPosParticle);
    x_drift = [x_drift dx(~isnan(dx))./dt(~isnan(dx))];
    
    dy = diff(spot_struct(i).yPosParticle);
    y_drift = [y_drift dy(~isnan(dy))./dt(~isnan(dy))];
    
    dz = diff(spot_struct(i).zPosParticle);
    z_drift = [z_drift dz(~isnan(dz))./dt(~isnan(dz))];
        
    t_drift = [t_drift dt(~isnan(dx))];
end    

% get pixel size
currExperiment = liveProject.includedExperiments{1};
PixelSize = currExperiment.pixelSize_nm/1e3;
zSize = 0.5;
p_dt = median(t_drift)/2;

% make figures
r_drift = sqrt(x_drift.^2+y_drift.^2)*PixelSize;
r_fig = figure;
cmap1 = brewermap([],'Set2');
histogram(r_drift,'Normalization','probability','FaceColor',cmap1(2,:))

xlabel('radial drift (\mu m per second)')
ylabel('probability')
set(gca,'Fontsize',14)
xlim([0 0.15])
grid on
saveas(r_fig,[FigurePath 'r_drift_hist.png'])

z_fig = figure;
histogram(z_drift*zSize,'Normalization','probability','FaceColor',cmap1(5,:))

xlabel('axial drift (\mu m per second)')
ylabel('probability')
set(gca,'Fontsize',14)
xlim([-0.15 0.15])
grid on
saveas(z_fig,[FigurePath 'z_drift_hist.png'])

