clear
close all

projectName = 'Bcd-GFP_hbMS2-mCh_Airy_fast_int';

liveProject = LiveEnrichmentProject(projectName);
resultsRoot = [liveProject.dataPath filesep];
FigurePath = [liveProject.figurePath filesep 'nucleus_movies'];
RefPath = [liveProject.dataPath 'refFrames' filesep];

% load data
load([resultsRoot 'spot_struct.mat'])  

% call header function
[liveProject, ~, dataName, hasAPInfo, has3DSpotInfo, hasProteinInfo, hasNucleusProbFiles] = headerFunction(projectName);


% determine protein channel 
proteinChannel = liveProject.includedExperiments{1}.inputChannels;

% select candidate nuclei for movies
% we need to do some homework here: want to identify true starts and stops
% in particles with reasonable segmentation that last a long time
minSpotDP = 75;
n_lead_frames = 3;
z_offset = 1;
% n_vec = [spot_struct.N];
% candidate_indices = find(n_vec>=minSpotDP);
spot_struct_trunc = struct;
qc_kernel = ones(1,5);
i_pass = 1;
transfer_fields = {'fluo','time','frames','yPosNucleus','xPosNucleus','APPosNucleus','FrameQCFlags'};
for i = 1:length(spot_struct)
    FrameQCFlags = spot_struct(i).FrameQCFlags;
    FrameQCFlags_1 = FrameQCFlags == 1;
    FrameQCFlags_2 = conv(qc_kernel,FrameQCFlags_1,'full');
    FrameQCFlags_2 = FrameQCFlags_2(3:end-2) >= 4;
    
    if nansum(FrameQCFlags_2) >= minSpotDP
        % calculate start and stop indices
        first_i = max([1,find(FrameQCFlags_2,1)-n_lead_frames]);
        last_i = min([length(FrameQCFlags_2),find(FrameQCFlags_2,1,'last')+n_lead_frames]);
        % add fields
        for f = 1:length(transfer_fields)
            vec = spot_struct(i).(transfer_fields{f});
            spot_struct_trunc(i_pass).(transfer_fields{f}) = vec(first_i:last_i);
        end
        % basic ID info
        spot_struct_trunc(i_pass).particleID = spot_struct(i).particleID;
        spot_struct_trunc(i_pass).setID = spot_struct(i).setID;
        spot_struct_trunc(i_pass).ncID = spot_struct(i).ncID;
        % pull kalman inference for particle position
        xPosInfDS = spot_struct(i).xPosInf(2:2:end);
        spot_struct_trunc(i_pass).xPosParticle = xPosInfDS(first_i:last_i)';
        yPosInfDS = spot_struct(i).yPosInf(2:2:end);
        spot_struct_trunc(i_pass).yPosParticle = yPosInfDS(first_i:last_i)';
        zPosInfDS = spot_struct(i).zPosInf(2:2:end);
        spot_struct_trunc(i_pass).zPosParticle = zPosInfDS(first_i:last_i)';
        spot_struct_trunc(i_pass).origID = i;
        % increment
        i_pass = i_pass + 1;
    end
end

%%
% set basic sampling parameters
nucleus_snip_radius_um = 3; % in microns
PixelSize = liveProject.includedExperiments{1}.pixelSize_um;
nucleus_snip_radius_px = round(nucleus_snip_radius_um/PixelSize);
snip_size = 2*nucleus_snip_radius_px+1;
int_radius = 0.25 / PixelSize
nMovies = 100; % number of nucleus movies


% randomize sampling order (eventually we can do this in parallel)
rng(312)
sample_order = randsample(1:length(spot_struct_trunc),nMovies,false);
spot_struct_input = spot_struct_trunc(sample_order);

% initialize data structure to store movie slices
movieData = struct;

% generate ref vectors
[RefStruct, SetFrameArray, SamplingResults] = generateReferenceVectors(spot_struct_input,RefPath,false,true,[]);
RefStruct.spot_fluo_ref = [spot_struct_input.fluo];
RefStruct.nc_ap_ref = [spot_struct_input.APPosNucleus];
RefStruct.time_ref = [spot_struct_input.time];

% iterate through all relevant set/frame combinations
NIter = size(RefStruct.set_frame_array,1);
tic
parpool(12);
parfor i_stack = 1:NIter
    
    % extract basic info    
    currFrame = RefStruct.set_frame_array(i_stack,2);
    currSetID = RefStruct.set_frame_array(i_stack,1);
    currExperiment = liveProject.includedExperiments{currSetID};
    Prefix = currExperiment.Prefix;
    % load protein stack
    proteinPath = [currExperiment.preFolder  Prefix '_' sprintf('%03d',currFrame) '_ch0' num2str(proteinChannel) '.tif'];
    proteinStack = imreadStack(proteinPath);
    
    % get indices of nucleus/frames in this stack
    frame_set_indices = find(RefStruct.setID_ref==currSetID & RefStruct.frame_ref==currFrame);
    
    % initialize fields
    movieData(i_stack).sliceArray = zeros(snip_size,snip_size,length(frame_set_indices));
    movieData(i_stack).x_pos_vec = RefStruct.spot_x_ref(frame_set_indices);
    movieData(i_stack).y_pos_vec = RefStruct.spot_y_ref(frame_set_indices);
    movieData(i_stack).ap_pos_vec = RefStruct.nc_ap_ref(frame_set_indices);
    movieData(i_stack).time_vec = RefStruct.time_ref(frame_set_indices);
    movieData(i_stack).frame_vec = RefStruct.frame_ref(frame_set_indices);
    movieData(i_stack).fluo_vec = RefStruct.spot_fluo_ref(frame_set_indices);
    movieData(i_stack).particleID_vec = RefStruct.particleID_ref(frame_set_indices);
    movieData(i_stack).x_nucleus = RefStruct.nc_x_ref(frame_set_indices);
    movieData(i_stack).y_nucleus = RefStruct.nc_y_ref(frame_set_indices);
    
    % iterate through individual nuclei
    for i_nucleus = 1:length(frame_set_indices)
        % nucleus position
        x_nucleus = RefStruct.nc_x_ref(frame_set_indices(i_nucleus));
        y_nucleus = RefStruct.nc_y_ref(frame_set_indices(i_nucleus));
        
        % spot position
%         x_spot = RefStruct.spot_x_ref(frame_set_indices(i_nucleus));
%         y_spot = RefStruct.spot_y_ref(frame_set_indices(i_nucleus));
        z_spot = max(min([currExperiment.zDim+2-1,round(RefStruct.spot_z_ref(frame_set_indices(i_nucleus)))]),2);
        
        % generate coordinate vectors for indexing
        x_vec_full = x_nucleus-nucleus_snip_radius_px:x_nucleus+nucleus_snip_radius_px;
        x_trim = x_vec_full(x_vec_full>=1 & x_vec_full<=currExperiment.xDim);
        y_vec_full = y_nucleus-nucleus_snip_radius_px:y_nucleus+nucleus_snip_radius_px;
        y_trim = y_vec_full(y_vec_full>=1 & y_vec_full<=currExperiment.yDim);
        
        % extract snip
        nucleus_snip = zeros(length(y_vec_full),length(x_vec_full));
        nucleus_snip(ismember(y_vec_full,y_trim),ismember(x_vec_full,x_trim)) = proteinStack(y_trim,x_trim,z_spot);
        
        % store info
        movieData(i_stack).sliceArray(:,:,i_nucleus) = nucleus_snip;
    end
        
end 
toc

%% save data
disp('Saving...')
save([resultsRoot 'nucleusMovieData.mat'],'movieData','-v7.3') 
disp('Done.')
%% make figures

% first, we need to reconfigure stuff to be in particle-centric
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
% parpool(12)
for p = 13%:length(particle_index)
  
    % make subdirectories
    prefix = ['particle_' num2str(particle_index(p))];
    cf_dir = [FigurePath filesep prefix '_contour'];
    mkdir(cf_dir);
    cf_mean_dir = [FigurePath filesep prefix '_mean_contour'];
    mkdir(cf_mean_dir);
    im_dir = [FigurePath filesep prefix '_heatmap'];
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
    pt_min = prctile(particle_slice(:),50);
    pt_max = prctile(particle_slice(:),99.9);
    
    pt_min_mean = prctile(particle_slice_mean(:),5);
    pt_max_mean = prctile(particle_slice_mean(:),99.9);
%     mc_vec = [];
    for f = 1:20%1:length(frame_vec)
      
         
        
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
            sz_val = spot_fluo(f);
            c_val = ceil(1e3*spot_fluo(f)/fluo_max);
            c_val = max([1 min(1e3,c_val)]);
            mc = cm1(c_val,:);
        end  
        
        hold on
        % plot protein heatmap
        contourf(particle_slice(:,:,f));
        % plot spot position
        scatter(rel_spot_x(f),rel_spot_y(f),20 + sz_val/fluo_max*150,'MarkerFaceColor',mc,'MarkerEdgeColor','k','MarkerFaceAlpha',1,'LineWidth',1.5)
        % plot circle indicating integration radisu
%         viscircles([rel_spot_x(f),rel_spot_y(f)],int_radius,'Color',[0.5 0.5 0.5 0.5],'LineStyle','--')
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
        
%         saveas(snip_fig1,[cf_dir filesep prefix '_contour_' sprintf('%03d',frame_vec(f)) '.tif'])
        
        %%%%%%%%%%%%%%%%%%
        % time-averaged contour
        snip_fig2 = figure;%('Visible','off');
       
        colormap(cm2)
        % calculate color            
        lead_vec = max([1 f-n_avg_frames]):f; 
        
        sz_vec = spot_fluo(lead_vec);
        sz_vec(isnan(sz_vec))=0;
        sz_vec(sz_vec<0)=0;
%         c_val_vec = ceil(1e3*spot_fluo(lead_vec)/fluo_max);
%         c_val_vec(isnan(c_val_vec)) = 1;
%         c_val_vec(c_val_vec > 1e3) = 1e3;
%         c_val_vec(c_val_vec < 1) = 1;
      
        
        hold on        
        contourf(particle_slice_mean(:,:,f))
        plot(rel_spot_x(lead_vec),rel_spot_y(lead_vec),'Color','k')
        scatter(rel_spot_x(lead_vec),rel_spot_y(lead_vec),20 + sz_vec/fluo_max*150,'MarkerFaceColor',gray,'MarkerEdgeAlpha',0.25,'MarkerFaceAlpha',0.25,'LineWidth',1,'MarkerEdgeColor','k')
        scatter(rel_spot_x(f),rel_spot_y(f),20 + sz_val/fluo_max*150,'MarkerFaceColor',mc,'MarkerEdgeColor','k','MarkerFaceAlpha',1,'LineWidth',1.5)
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
        
%         saveas(snip_fig2,[cf_mean_dir filesep  prefix '_mean_contour_' sprintf('%03d',frame_vec(f)) '.tif'])
        
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
        scatter(rel_spot_x(f),rel_spot_y(f),20 + sz_val/fluo_max*150,'MarkerFaceColor',mc,'MarkerEdgeColor','k','MarkerFaceAlpha',1,'LineWidth',1.5)
        % plot circle indicating integration radisu
%         viscircles([rel_spot_x(f),rel_spot_y(f)],int_radius,'Color',[0.5 0.5 0.5 0.5],'LineStyle','--')
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
        
%         saveas(snip_fig3,[cf_dir filesep prefix '_heatmap_' sprintf('%03d',frame_vec(f)) '.tif'])
        pause(0.1)
    end
end


