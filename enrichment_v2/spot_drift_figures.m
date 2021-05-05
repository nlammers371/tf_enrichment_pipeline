% script to make x-y-z drift histograms
clear
close all

projectName = 'Bcd-GFP_hbMS2-mCh_AiryscanTest_';

liveProject = LiveEnrichmentProject(projectName);
resultsRoot = [liveProject.dataPath filesep];
FigurePath = liveProject.figurePath;


projectNameCell = {'Bcd-GFP_hbMS2-mCh_AiryscanTest_','Bcd-GFP_hbP2P-mCh'};
master_struct = struct;

for p = 1:length(projectNameCell)
    projectName = projectNameCell{p};
    liveProject = LiveEnrichmentProject(projectName);
    resultsRoot = [liveProject.dataPath filesep];

    % load data
    load([resultsRoot 'spot_struct.mat'])
    load([resultsRoot 'snip_data.mat'])
    master_struct(p).spot_struct = spot_struct;
    master_struct(p).snip_data = snip_data;
    master_struct(p).projectName = projectName;
    clear spot_struct
end

%% 
close all
frame_ind = 42;
slice_ind_list = [12:15];

for s = 1:length(slice_ind_list)
    a_fig = figure;
    imagesc(master_struct(1).snip_data(frame_ind).spot_protein_snips(:,:,slice_ind_list(s)))
    saveas(a_fig,[FigurePath 'airy_slice_' num2str(s) '.png'])
    
    l_fig = figure;
    imagesc(master_struct(2).snip_data(frame_ind).spot_protein_snips(:,:,slice_ind_list(s)))
    saveas(l_fig,[FigurePath 'leica_slice_' num2str(s) '.png'])
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

