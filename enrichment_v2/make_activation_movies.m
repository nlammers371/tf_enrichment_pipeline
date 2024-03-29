clear
close all

projectName = 'Bcd-GFP_hbMS2-mCh_Airy_fast_int';

liveProject = LiveEnrichmentProject(projectName);
resultsRoot = [liveProject.dataPath filesep];
FigurePath = [liveProject.figurePath filesep 'activation_movies'];
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
try
  parpool(18)
catch
end
for p = 60%:length(particle_index)
  
    % make subdirectories
    prefix = ['particle_' num2str(particle_index(p))];
    cf_dir = [FigurePath filesep prefix '_contour'];
    mkdir(cf_dir);
    cf_mean_dir = [FigurePath filesep prefix '_mean_contour'];
    mkdir(cf_mean_dir);
    hm_dir = [FigurePath filesep prefix '_heatmap'];
    mkdir(hm_dir);
    im_dir = [FigurePath filesep prefix '_raw'];
    mkdir(im_dir);
    
    % apply particle filter
    particle_filter = longParticleArray==particle_index(p);
    particle_slice = longSliceArray(:,:,particle_filter);
    rel_spot_x = longSpotXPosArray(particle_filter) - longNucXArray(particle_filter) + nucleus_snip_radius_px + 1;
    rel_spot_y = longSpotYPosArray(particle_filter) - longNucYArray(particle_filter) + nucleus_snip_radius_px + 1;
    spot_fluo = longFluoArray(particle_filter);
    frame_vec = longFrameArray(particle_filter);
    time_vec = longTimeArray(particle_filter);
    
    % use convolution to generate a lagged protein average indicator
    conv_kernel = ones(1,1,n_avg_frames);
    dummy_stack = ones(size(particle_slice));
    particle_slice_mean = convn(conv_kernel,particle_slice);
    norm_stack = convn(conv_kernel,dummy_stack);
    particle_slice_mean = particle_slice_mean ./ norm_stack;
    particle_slice_mean = particle_slice_mean(:,:,1:size(particle_slice_mean,3)-length(conv_kernel)+1);
    
    % set fluorescence scale
    lb = prctile(particle_slice(:),0.1);
    ub = prctile(particle_slice(:),99.9);
    pt_min = prctile(particle_slice(:),40);
    pt_max = prctile(particle_slice(:),99.9);
    
    pt_min_mean = prctile(particle_slice_mean(:),5);
    pt_max_mean = prctile(particle_slice_mean(:),99.9);

    for f = 1:length(frame_vec)
      
        % first write the raw tif snip to file
%         im_snip = mat2gray(particle_slice(:,:,f),[lb ub]);
%         imwrite(im_snip,[im_dir filesep prefix '_snip_' sprintf('%03d',frame_vec(f)) '.tif'])
        
        %%%%%%%%%%%%%%%%%%
        % contour plot
        snip_fig1 = figure('Visible','off');
        cm1 = brewermap(1e3,'Greens');
        cm2 = flipud(brewermap(1e3,'RdYlBu'));
        colormap(cm2)              
          
        % calculate color
        mc = gray/2;
        sz_val = 0;
        if ~isnan(spot_fluo(f))
            sz_val = max([1 spot_fluo(f)]);
            c_val = ceil(1e3*spot_fluo(f)/fluo_max);
            c_val = max([1 min(1e3,c_val)]);
            mc = cm1(c_val,:);
        end  
        
        hold on
        % plot protein heatmap
        cf = contourf(particle_slice(:,:,f),8);
        % plot spot position
        scatter(rel_spot_x(f),rel_spot_y(f),10 + (sz_val/fluo_max*70)^1.3,'MarkerFaceColor',mc,'MarkerEdgeColor','k','MarkerFaceAlpha',1,'LineWidth',1.5)
        % plot circle indicating integration radisu
        th = 0:pi/50:2*pi;
        xunit = int_radius * cos(th) + rel_spot_x(f);
        yunit = int_radius * sin(th) + rel_spot_y(f);
        plot(xunit, yunit,'Color','k');
        
        % make colorbar
        caxis([pt_min,pt_max])
        h = colorbar;
        
        % formatng stuff
        set(gca,'Fontsize',14);
        xlim([1 snip_size])
        ylim([1 snip_size])
        box on
        
        % labels
        tick_vec = round([1/PixelSize 3/PixelSize 5/PixelSize]);
        tick_labels = [-2 0 2];
        set(gca,'xtick',tick_vec,'xticklabels',tick_labels)
        set(gca,'ytick',tick_vec,'yticklabels',tick_labels)
        
        xlabel('\mum')
        ylabel('\mum')
        ylabel(h,'Bcd concentration (au)')
        
        text(2,4,[num2str(round(time_vec(f)/60,2)) ' minutes'],'Fontsize',14,'Color','w')
        
        saveas(snip_fig1,[cf_dir filesep prefix '_contour_' sprintf('%03d',frame_vec(f)) '.tif'])
        
        %%%%%%%%%%%%%%%%%%
        % time-averaged contour
        %%%%%%%%%%%%%%%%%%%%
        
        snip_fig2 = figure('Visible','off');
       
        colormap(cm2)
        % calculate color            
        lead_vec = max([1 f-n_avg_frames]):f; 
        
        sz_vec = spot_fluo(lead_vec);
        sz_vec(isnan(sz_vec))=0;
        sz_vec(sz_vec<0)=0;      
        
        hold on        
        contourf(particle_slice_mean(:,:,f),8)
        plot(rel_spot_x(lead_vec),rel_spot_y(lead_vec),'Color','k')
        scatter(rel_spot_x(lead_vec),rel_spot_y(lead_vec),10 + (sz_vec/fluo_max*70).^1.3,'MarkerFaceColor',gray,'MarkerEdgeAlpha',0.25,'MarkerFaceAlpha',0.25,'LineWidth',1,'MarkerEdgeColor','k')
        scatter(rel_spot_x(f),rel_spot_y(f),10 + (sz_val/fluo_max*70).^1.3,'MarkerFaceColor',mc,'MarkerEdgeColor','k','MarkerFaceAlpha',1,'LineWidth',1.5)
        % plot circle indicating integration radisu
%         viscircles([rel_spot_x(f),rel_spot_y(f)],int_radius,'Color',[0.5 0.5 0.5 0.5],'LineStyle','--')
        th = 0:pi/50:2*pi;
        xunit = int_radius * cos(th) + rel_spot_x(f);
        yunit = int_radius * sin(th) + rel_spot_y(f);
        plot(xunit, yunit,'Color','k');
        
        % make colorbar
        caxis([pt_min_mean,pt_max_mean])
        h = colorbar;
        
        % formatng stuff
        set(gca,'Fontsize',14);
        xlim([1 snip_size])
        ylim([1 snip_size])
        box on
        
        % labels
        tick_vec = round([1/PixelSize 3/PixelSize 5/PixelSize]);
        tick_labels = [-2 0 2];
        set(gca,'xtick',tick_vec,'xticklabels',tick_labels)
        set(gca,'ytick',tick_vec,'yticklabels',tick_labels)
        
        xlabel('\mum')
        ylabel('\mum')
        ylabel(h,'Bcd concentration (au)')
        
        text(2,4,[num2str(round(time_vec(f)/60,2)) ' minutes'],'Fontsize',14,'Color','w')
        
        saveas(snip_fig2,[cf_mean_dir filesep  prefix '_mean_contour_' sprintf('%03d',frame_vec(f)) '.tif'])
        
        pause(0.25)
%         close all

        snip_fig3 = figure('Visible','off');
        cm1 = brewermap(1e3,'Greens');
        cm2 = flipud(brewermap(32,'RdYlBu'));
        colormap(cm2)              
          
        hold on
        % plot protein heatmap
        pc = pcolor(particle_slice(:,:,f));
        pc.EdgeAlpha = 0.1;
        % plot spot position
        scatter(rel_spot_x(f),rel_spot_y(f),10 + (sz_val/fluo_max*70)^1.3,'MarkerFaceColor',mc,'MarkerEdgeColor','k','MarkerFaceAlpha',1,'LineWidth',1.5)
        % plot circle indicating integration radisu
        th = 0:pi/50:2*pi;
        xunit = int_radius * cos(th) + rel_spot_x(f);
        yunit = int_radius * sin(th) + rel_spot_y(f);
        plot(xunit, yunit,'Color','k');
        
        % make colorbar
        caxis([pt_min,pt_max])
        h = colorbar;
        
        % formatng stuff
        set(gca,'Fontsize',14);
        xlim([1 snip_size])
        ylim([1 snip_size])
        box on
        
        % labels
        tick_vec = round([1/PixelSize 3/PixelSize 5/PixelSize]);
        tick_labels = [-2 0 2];
        set(gca,'xtick',tick_vec,'xticklabels',tick_labels)
        set(gca,'ytick',tick_vec,'yticklabels',tick_labels)
        
        xlabel('\mum')
        ylabel('\mum')
        ylabel(h,'Bcd concentration (au)')
        
        text(2,4,[num2str(round(time_vec(f)/60,2)) ' minutes'],'Fontsize',14,'Color','w')
        
        saveas(snip_fig3,[hm_dir filesep prefix '_heatmap_' sprintf('%03d',frame_vec(f)) '.tif'])
        
        close all
    end
end


