clear
close all

projectName = 'Bcd-GFP_hbMS2-mCh_Airy_fast_int';

liveProject = LiveEnrichmentProject(projectName);
resultsRoot = [liveProject.dataPath filesep];
FigurePath = [liveProject.figurePath filesep 'nucleus_movies'];
RefPath = [liveProject.dataPath 'refFrames' filesep];

load([resultsRoot 'nucleusMovieData.mat'],'movieData') 

%% first, we need to reconfigure stuff to be in particle-centric

nucleus_snip_radius_um = 3; % in microns
PixelSize = liveProject.includedExperiments{1}.pixelSize_um;
nucleus_snip_radius_px = round(nucleus_snip_radius_um/PixelSize);
snip_size = 2*nucleus_snip_radius_px+1;
int_radius = 0.25 / PixelSize;
nMovies = 100; % number of nucleus movies

% arrangements
longSliceArray = cat(3,movieData.sliceArray);
longSpotXPosArray = [movieData.x_pos_vec];
longSpotYPosArray = [movieData.y_pos_vec];
longAPPosArray = [movieData.ap_pos_vec];
longTimeArray = [movieData.time_vec];
longFrameArray = [movieData.frame_vec];
longFluoArray = [movieData.fluo_vec];
longParticleArray = [movieData.particleID_vec];
longNucXArray = [movieData.x_nucleus];
longNucYArray = [movieData.y_nucleus];

% iterate through particles
particle_index = unique(longParticleArray);

% calculate fluorescence scale
fluo_max = prctile(longFluoArray,90);
markerSize = 75;
n_avg_frames = 7;
gray = [0.5 0.5 0.5];

% set frame and particle to image
pID = 20;
frames_plot = [28 30];

% make subdirectories
prefix = ['particle_' num2str(particle_index(pID))];
fig_dir = [FigurePath filesep prefix '_surf'];
mkdir(fig_dir);

% apply particle filter
particle_filter = longParticleArray==particle_index(pID);
particle_slice = longSliceArray(:,:,particle_filter);
rel_spot_x = longSpotXPosArray(particle_filter) - longNucXArray(particle_filter) + nucleus_snip_radius_px + 1;
rel_spot_y = longSpotYPosArray(particle_filter) - longNucYArray(particle_filter) + nucleus_snip_radius_px + 1;
spot_fluo = longFluoArray(particle_filter);
frame_vec = longFrameArray(particle_filter);
time_vec = longTimeArray(particle_filter);                        
     
% set fluorescence scale
pt_min = prctile(particle_slice(:),50);
pt_max = prctile(particle_slice(:),99.9);
    
%%%%%%%%%%%%%%%%%%
% surface plot
close all
for frame = frames_plot
    surf_fig = figure;
    cm1 = brewermap(9,'Greens');
    cm2 = flipud(brewermap(1e3,'RdYlBu'));
    colormap(cm2)              

    hold on
    % plot protein heatmap
    bcd_slice = particle_slice(:,:,frame);
    sf = surf(bcd_slice,'EdgeAlpha',0.1);
    view(-45,55)
    % plot spot position
    % get spot protein coordinate
    bcd_val = bcd_slice(round(rel_spot_y(frame)),round(rel_spot_x(frame)));
    scatter3(rel_spot_x(frame),rel_spot_y(frame),bcd_val,10 + sz_val/fluo_max*60,'MarkerFaceColor',cm1(6,:),'MarkerEdgeColor','k','MarkerFaceAlpha',1,'LineWidth',1)


    % make colorbar
    caxis([pt_min,pt_max])
%     h = colorbar;

    % formatng stuff
    set(gca,'Fontsize',14);
    xlim([1 snip_size])
    ylim([1 snip_size])

    grid on
    % labels
    tick_vec = round([1/PixelSize 3/PixelSize 5/PixelSize]);
    tick_labels = [-2 0 2];
    set(gca,'xtick',tick_vec,'xticklabels',tick_labels)
    set(gca,'ytick',tick_vec,'yticklabels',tick_labels)

    xlabel('\mum')
    ylabel('\mum')
    zlabel('[Bcd] (au)')
%     ylabel(h,'Bcd concentration (au)')

%     text(2,4,[num2str(round(time_vec(frame)/60,2)) ' minutes'],'Fontsize',14,'Color','w')

    saveas(surf_fig,[fig_dir filesep prefix '_contour_' sprintf('%03d',frame) '.tif'])
end      