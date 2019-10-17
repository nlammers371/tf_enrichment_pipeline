% Script to make movie showing how averaging over time and space reveals a
% local enrichment signal
clear
close all
addpath('utilities')

distLim = 0.8; % min distance of spot from edge permitted (um)
colormap_heat = viridis(128); %specifies the colormap used to make heatmaps

rawEnrichHeatmap_ub = 1.75;
rawEnrichHeatmap_lb = 1.25;
% absDiffEnrichHeatmap_ub = 0.5;
% absDiffEnrichHeatmap_lb = 0;
% relEnrichHeatmap_ub = 1.2;
% relEnrichHeatmap_lb = 1;

project = 'Dl-Ven_snaBAC-mCh';

dropboxFolder = 'E:\Nick\LivemRNA\Dropbox (Personal)\';
dataPath = [dropboxFolder '\ProcessedEnrichmentData\' project '\'];
figPath = ['E:\Meghan\Dropbox' '\LocalEnrichmentFigures\PipelineOutput\' project '\'];

writePath = [figPath '\avg_snips_movie_frames\'];
mkdir(writePath)

% paperFigPath = [figPath 'paperFigs/'];
% basicFigPath = [figPath 'basicFigs/'];

% Load analysis data
load([dataPath 'nucleus_struct_protein.mat'], 'nucleus_struct_protein');

%% Extract and Process Snip Stacks

% Extract snip stacks
PixelSize = nucleus_struct_protein(1).PixelSize;
spot_protein_snips = cat(3,nucleus_struct_protein.spot_protein_snips);
null_protein_snips = cat(3,nucleus_struct_protein.edge_null_protein_snips);

% Make r reference array where each element is its distance, in um, from 
% the center pixel in the array
snip_size = size(spot_protein_snips,1);
[y_ref, x_ref] = meshgrid(1:snip_size,1:snip_size);
r_ref = sqrt((x_ref - ceil(snip_size/2)).^2 + (y_ref - ceil(snip_size/2)).^2)*PixelSize;

spot_protein_vec = [nucleus_struct_protein.spot_protein_vec];
null_protein_vec = [nucleus_struct_protein.edge_null_protein_vec];
dist_vec = [nucleus_struct_protein.spot_edge_dist_vec]*PixelSize;

% Apply distance filter
spot_protein_snips_dist = spot_protein_snips(:,:,dist_vec>=distLim);
null_protein_snips_dist = null_protein_snips(:,:,dist_vec>=distLim);

% Randomize snip orientation (via reflections)
inv_mat = [fliplr(1:snip_size); 1:snip_size]' ;
spot_protein_snips_reflected = NaN(size(spot_protein_snips_dist));
null_protein_snips_reflected = NaN(size(spot_protein_snips_dist));
for i = 1:size(spot_protein_snips_dist,3)
    h = ceil(rand()*2);
    v = ceil(rand()*2);
    spot_protein_snips_reflected(:,:,i) = spot_protein_snips_dist(inv_mat(:,v),inv_mat(:,h),i);
    null_protein_snips_reflected(:,:,i) = null_protein_snips_dist(inv_mat(:,v),inv_mat(:,h),i);
end

% Shuffle the order of the snips in the stack
index_vec = 1:size(spot_protein_snips_reflected,3);
rand_order_vec = randsample(index_vec,numel(index_vec));
spot_protein_snips_mixed = spot_protein_snips_reflected(:,:,rand_order_vec);
null_protein_snips_mixed = null_protein_snips_reflected(:,:,rand_order_vec);   %Putting null snips in same order as spot snips

spot_protein_full_mean = nanmean(spot_protein_snips_mixed,3);
null_protein_full_mean = nanmean(null_protein_snips_mixed,3);

%% Make Averaged Heatmap Movie

% Set movie parameters
frame_steps = unique(round(logspace(0,2,50)*8));    %movie with 5 to 800 frames averaged
visibleOn = false; % Don't want to display figures as they're made

clabel_raw = 'Dorsal total protein (au)';
title_raw = 'Dorsal at {\itsnail} - total protein (au)';
% clabel_rel = 'Dorsal total protein (au)';
% title_rel = 'Dorsal at {\itsnail} - relative enrichment';
% clabel_absDiff = 'Dorsal absolute enrichment (au)';
% title_absDiff = 'Dorsal at {\itsnail} - absolute enrichment (au)';

% Create & open the video writer with 2 fps
meanFrameWriter = VideoWriter([writePath 'meanSnips.avi'],'Uncompressed AVI');
meanFrameWriter.FrameRate = 2;
open(meanFrameWriter);

for n = 1:numel(frame_steps)  %(end_frame/frame_incr + 1)
    %
    % Make averaged frame version of the movie
    %
    curr_max_frame = frame_steps(n);%min([frame_incr*n,numel(index_vec)]);
    mean_frame_raw = nanmean(spot_protein_snips_mixed(:,:,1:curr_max_frame),3);
    null_mean_frame = nanmean(null_protein_snips_mixed(:,:,1:curr_max_frame),3);
    mean_frame_rel = mean_frame_raw ./ null_mean_frame; %Relative difference (fold) enrichment
    mean_frame_absDiff = mean_frame_raw - null_mean_frame;  %Absolute difference enrichment
    
    % Make averaged figures for raw enrichment (total protein)
    temp_fig_raw = makeHeatmapPlots(mean_frame_raw, visibleOn, ...
                  '', clabel_raw, colormap_heat,PixelSize,...
                  rawEnrichHeatmap_lb,rawEnrichHeatmap_ub);
    set(gca,'xcolor','black','ycolor','black')
    set(gcf,'color','white');
    text(0.8,1.9,[num2str(frame_steps(n)) ' samples'],'Color','black', ...
        'BackgroundColor','white', 'FontSize',15,'FontName','Lucida Sans')
    % Write the current frame to the movie
    writeVideo(meanFrameWriter, getframe(gcf));
    close all
    
end
close(meanFrameWriter);

%% Make single-frame heatmap movie
frame_incr = 25;
end_frame = 1000;
n_frames = ceil(size(spot_protein_snips_mixed,3)/frame_incr);

% Create & open the video writer with 2 fps
singleFrameWriter = VideoWriter([writePath 'singleFrameSnips.avi'],'Uncompressed AVI');
singleFrameWriter.FrameRate = 2;
open(singleFrameWriter);
for n = 1:(end_frame/frame_incr + 1)
    %
    % Make single-frame version of the movie
    %
    curr_single_frame_raw = min([numel(index_vec),n*frame_incr]);
    spot_single_frame_raw = spot_protein_snips_mixed(:,:,curr_single_frame_raw);
    null_single_frame = null_protein_snips_mixed(:,:,curr_single_frame_raw);
    single_frame_rel = spot_single_frame_raw ./ null_single_frame; %Relative difference (fold) enrichment
    single_frame_absDiff = spot_single_frame_raw - null_single_frame;  %Absolute difference enrichment
    
    % Make single-frame figures for raw enrichment (total protein)
    fig_raw_single = makeHeatmapPlots(spot_single_frame_raw, visibleOn, ...
                  '',clabel_raw, colormap_heat,PixelSize,0,4);
    set(gca,'xcolor','black','ycolor','black')
    set(gcf,'color','white');
    % Write the current frame to the movie
    writeVideo(singleFrameWriter, getframe(gcf));
    close all
end
close(singleFrameWriter);