% Script to examine movements of spots within nuclei
clear 
close all
% Define ID and path variables
project = 'Bcd_GFP_hb_mCherry_Zoom3x_LowPower';
DataPath = ['../../dat/' project '/'];

% Load data structures
load([DataPath  'nucleus_struct.mat']);
load([DataPath  'ref_frame_struct.mat']);

% Get path to ref images
RefPath = [DataPath '/mf_images/'];
nc_set_vec = [nucleus_struct.setID];
set_index = unique(nc_set_vec);

% Get pixel sizes
px_sizes = [];
for s = 1:numel(set_index)
    px = [nucleus_struct(nc_set_vec==set_index(s)).PixelSize];
    px_sizes = [px_sizes px(1)];
end

snip_size = round(5 / mean(px_sizes));
% Vectors to index ref frame
ref_frame_sets = [ref_frame_struct.set]; 
ref_frame_frames = [ref_frame_struct.frame];

src_cell = {1,numel(set_index)};
pt_id_vec = [nucleus_struct.ParticleID];
for s = 1:numel(set_index)
    ind = find([nucleus_struct.setID]==set_index(s)&~isnan(pt_id_vec),1);   
    src = nucleus_struct(ind).source_path;    
    src_cell{s} = src;     
end
meta_iter = 1;
particle_trace_struct = struct;
for i = 100:200%numel(nucleus_struct)
    % basic id info
    ParticleID = nucleus_struct(i).ParticleID;  
    if isnan(ParticleID)
        continue
    end
    nc_ind = find(pt_id_vec==ParticleID);
    setID = floor(ParticleID);
    src = src_cell{set_index==floor(ParticleID)};    
    % get frame and position info
    frame_vec = nucleus_struct(i).frames;     
    nc_x_vec = nucleus_struct(i).xPos;     
    nc_y_vec = nucleus_struct(i).yPos;
    pt_x_vec = nucleus_struct(i).xPosParticle;     
    pt_y_vec = nucleus_struct(i).yPosParticle;
    pt_z_vec = nucleus_struct(i).brightestZs;
    
    pt_indices = find(~isnan(pt_x_vec));
    % initialize arrays to store particle info 
    nc_boundary_cell = cell(numel(pt_indices),2);
    iter = 0;
    % Iterate through frames
    for f = pt_indices
        iter = iter + 1;
        frame = frame_vec(f);        
        % load
        ref_name = [RefPath '/ref_frame_' sprintf('%03d',frame) '_set' sprintf('%03d',setID) '.mat'];
        load(ref_name,'rf');    
        
        x_particle = pt_x_vec(f);
        y_particle = pt_y_vec(f);
        x_nucleus = nc_x_vec(f);
        y_nucleus = nc_y_vec(f);
                
        ID = rf(y_particle,x_particle); 
        if isnan(ID) 
            continue
        end
        ref_mask = rf==ID;
        B = bwboundaries(ref_mask);
        nc_boundary_cell{iter,1} = B{1}(:,1) - y_nucleus;
        nc_boundary_cell{iter,2} = B{1}(:,2) - x_nucleus;        
    end
    particle_trace_struct(meta_iter).xPosParticle = pt_x_vec(~isnan(pt_x_vec)) - nc_x_vec(~isnan(pt_x_vec));
    particle_trace_struct(meta_iter).yPosParticle = pt_y_vec(~isnan(pt_x_vec)) - nc_y_vec(~isnan(pt_x_vec));    
    particle_trace_struct(meta_iter).frames = frame_vec(~isnan(pt_x_vec));    
    particle_trace_struct(meta_iter).nc_boundary_cell = nc_boundary_cell;
    meta_iter = meta_iter + 1;
end
% make plots
cm = jet(128);
for i = 1:numel(particle_trace_struct)
    xp_vec = particle_trace_struct(i).xPosParticle;
    yp_vec = particle_trace_struct(i).yPosParticle;
    frame_vec = particle_trace_struct(i).frames;
    nc_boundary_cell = particle_trace_struct(i).nc_boundary_cell;
    
    inc_vec = floor(frame_vec / max(frame_vec) * 128);
    pt_fig = figure;
    hold on
    plot(xp_vec,yp_vec,'Color','black')
    for j = 1:size(nc_boundary_cell,1)
        nc_x = nc_boundary_cell{j,2};
        nc_y = nc_boundary_cell{j,1};
        if isempty(nc_x)
            continue
        end
        plot(nc_x, nc_y, 'Color', [cm(inc_vec(j),:) .2],'LineWidth',2)
        scatter(xp_vec(j), yp_vec(j), 10, 'MarkerFaceColor', cm(inc_vec(j),:),...
                'MarkerEdgeAlpha', 0)
    end
    
end
    