% script to make x-y-z drift histograms
clear
close all

projectName = 'Bcd-GFP_hbMS2-mCh_Airy_fast';

liveProject = LiveEnrichmentProject(projectName);
resultsRoot = [liveProject.dataPath filesep];
FigurePath = liveProject.figurePath;


projectNameCell = {'Bcd-GFP_hbMS2-mCh_Airy_fast','Bcd-GFP_hbMS2-mCh_NoAiry_02'};
master_struct = struct;

for p = 1:length(projectNameCell)
    projectName = projectNameCell{p};
    liveProject = LiveEnrichmentProject(projectName);
    resultsRoot = [liveProject.dataPath filesep];
    
    % load data
    load([resultsRoot 'spot_struct.mat'])
    load([resultsRoot 'proteinSamplingInfo.mat'])
    load([resultsRoot 'snip_data.mat'])
    master_struct(p).spot_struct = spot_struct;
    master_struct(p).pixelSize = liveProject.includedExperiments{1}.pixelSize_um;
    master_struct(p).snip_data = snip_data;
    master_struct(p).projectName = projectName;
    clear spot_struct
end

%% 
close all
frame_ind = 60;
slice_ind_list = [3:5];

for s = 1:length(slice_ind_list)
    a_fig = figure;
    airy_snip = master_struct(1).snip_data(frame_ind).spot_protein_snips(:,:,slice_ind_list(s));
    imagesc(airy_snip)
    c = colorbar;
    ylabel('position (\mum)','FontSize',14)
    xlabel('position (\mum)','FontSize',14)
    ylabel(c,'Bcd-GFP concentration (AU)','FontSize',14)
    xy_lim = (size(airy_snip,1) + 1) / 2 * master_struct(1).pixelSize; % um, calculated with center of center pixel as 0 
    set(gca,'xtick',([-1.0 0 1.0] + xy_lim) ./ master_struct(1).pixelSize, ...
            'xticklabel',[-1.0, 0, 1.0])
    set(gca,'ytick',([-1.0 0 1.0] + xy_lim) ./ master_struct(1).pixelSize, ...
            'yticklabel',[-1.0, 0, 1.0])
    saveas(a_fig,[FigurePath 'airy_980_slice_' num2str(s) '.png'])
%     
    l_fig = figure;
    conf_snip = master_struct(2).snip_data(frame_ind).spot_protein_snips(:,:,slice_ind_list(s));
    imagesc(conf_snip)
    c = colorbar;
    ylabel('position (\mum)','FontSize',14)
    xlabel('position (\mum)','FontSize',14)
    ylabel(c,'Bcd-GFP concentration (AU)','FontSize',14)
    xy_lim = (size(airy_snip,1) + 1) / 2 * master_struct(2).pixelSize; % um, calculated with center of center pixel as 0 
    set(gca,'xtick',([-1.0 0 1.0] + xy_lim) ./ master_struct(2).pixelSize, ...
            'xticklabel',[-1.0, 0, 1.0])
    set(gca,'ytick',([-1.0 0 1.0] + xy_lim) ./ master_struct(2).pixelSize, ...
            'yticklabel',[-1.0, 0, 1.0])
    saveas(l_fig,[FigurePath 'conf_980_slice_' num2str(s) '.png'])

    l_fig = figure;
    conf_snip = imgaussfilt(master_struct(2).snip_data(frame_ind).spot_protein_snips(:,:,slice_ind_list(s)),1);
    imagesc(conf_snip)
    c = colorbar;
    ylabel('position (\mum)','FontSize',14)
    xlabel('position (\mum)','FontSize',14)
    ylabel(c,'Bcd-GFP concentration (AU)','FontSize',14)
    xy_lim = (size(airy_snip,1) + 1) / 2 * master_struct(2).pixelSize; % um, calculated with center of center pixel as 0 
    set(gca,'xtick',([-1.0 0 1.0] + xy_lim) ./ master_struct(2).pixelSize, ...
            'xticklabel',[-1.0, 0, 1.0])
    set(gca,'ytick',([-1.0 0 1.0] + xy_lim) ./ master_struct(2).pixelSize, ...
            'yticklabel',[-1.0, 0, 1.0])
    saveas(l_fig,[FigurePath 'conf_980_slice_' num2str(s) '_sm.png'])
end

%%
ind = 101;

xy_fig = figure;
xvec = master_struct(1).spot_struct(ind).xPosParticle;
yvec = master_struct(1).spot_struct(ind).yPosParticle;

plot(xvec,yvec,'-o')
xlabel('x position (pixels)')
xlabel('y position (pixels)')
grid on
saveas(xy_fig,[FigurePath 'xy_particle_path.png'])

time = master_struct(1).spot_struct(ind).time;
t_interp = time(1:end-1) + diff(time)/2;
x_interp = interp1(time,xvec,t_interp,'linear',0);
y_interp = interp1(time,yvec,t_interp,'linear',0);

xy_fig = figure;
hold on
plot(xvec,yvec,'-o')
scatter(x_interp,y_interp,20,'filled')
xlabel('x position (pixels)')
xlabel('y position (pixels)')
grid on
saveas(xy_fig,[FigurePath 'xy_particle_path_int.png'])

