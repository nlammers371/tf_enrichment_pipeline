% Script to investigate spatio-temporal dynamics of protein distributions
% within nuclei. Ultimate goal is to find a way to infer position of active
% loci using protein channel alone
clear
close all

% define core ID variables
project = 'Dl-Ven_snaBAC-mCh';
% project = 'Dl-Ven_hbP2P-mCh';
dropboxFolder =  'E:\Nick\LivemRNA\Dropbox\';
dataPath = [dropboxFolder 'ProcessedEnrichmentData\' project '/'];
rawPath = 'E:\LocalEnrichment\Data\PreProcessedData\';
refFramePath = [dataPath 'refFrames/'];
ref_frame_list = dir([refFramePath 'nc_ref*']);
proteinChannel = 1;

% specify number of snips
n_snips = 10000;

% Load trace data
load([dataPath '/nucleus_struct.mat'],'nucleus_struct')
load([dataPath '/set_key.mat'],'set_key')
figPath = [dropboxFolder 'LocalEnrichmentFigures/spot_location_prediction/'];
mkdir(figPath)
addpath('./utilities')
% get MCP channel
options = [1 2];
mcp_channel = options(options~=proteinChannel);


%%% Generate reference vectors
nc_x_ref_vec = [nucleus_struct.xPos];
nc_y_ref_vec = [nucleus_struct.yPos];
spot_x_ref_vec = [nucleus_struct.xPosParticle];
spot_y_ref_vec = [nucleus_struct.yPosParticle];
spot_z_ref_vec = [nucleus_struct.zPosParticle];
fluo_ref_vec = [nucleus_struct.fluo];
frame_ref_vec = [nucleus_struct.frames];
time_ref_vec = [nucleus_struct.time];
mf_protein_vec = NaN(size(fluo_ref_vec));
set_ref_vec = NaN(size(fluo_ref_vec));
particle_ref_vec = NaN(size(fluo_ref_vec));
iter = 1;
for i = 1:numel(nucleus_struct)
    fluo = nucleus_struct(i).fluo;
    ptID = nucleus_struct(i).ParticleID;
    set_ref_vec(iter:iter+numel(fluo)) = floor(ptID);
    particle_ref_vec(iter:iter+numel(fluo)-1) = ptID;
    iter = iter + numel(fluo);
end

set_index = unique(set_ref_vec);
% determine size of neighborhood to use
px_size = nucleus_struct(1).PixelSize;
nb_size = round(3 ./ px_size);
xDim = nucleus_struct(1).xDim;
yDim = nucleus_struct(1).yDim;
zDim = nucleus_struct(1).zDim;
% filter for spots that meet criteria
analysis_ft = spot_z_ref_vec > 1 & spot_z_ref_vec < zDim & nc_x_ref_vec > nb_size ...
    & nc_x_ref_vec <= xDim-nb_size & nc_y_ref_vec > nb_size & nc_y_ref_vec ...
    <= yDim-nb_size & fluo_ref_vec > prctile(fluo_ref_vec,25);
% select indices to use for testing and training
rng(123);
index_list = randsample(find(analysis_ft),n_snips,false);
nc_x_ref_vec = nc_x_ref_vec(index_list);
nc_y_ref_vec = nc_y_ref_vec(index_list);
spot_x_ref_vec = spot_x_ref_vec(index_list);
spot_y_ref_vec = spot_y_ref_vec(index_list);
spot_z_ref_vec = spot_z_ref_vec(index_list);
time_ref_vec = time_ref_vec(index_list);
fluo_ref_vec = fluo_ref_vec(index_list);
particle_ref_vec = particle_ref_vec(index_list);
set_ref_vec = set_ref_vec(index_list);
frame_ref_vec = frame_ref_vec(index_list);
% make nucleus mask matrix
nc_mask = false(2*nb_size+1,2*nb_size+1);
nc_rad = 23;
[x_ref,y_ref] = meshgrid(1:2*nb_size+1,1:2*nb_size+1);
r_mat = sqrt((x_ref-nb_size-1).^2+(y_ref-nb_size-1).^2);
nc_mask(r_mat<=nc_rad) = true;
% initialize arrays to store info
nc_protein_array = NaN(2*nb_size+1,2*nb_size+1,1,numel(index_list));
nc_mcp_array = NaN(2*nb_size+1,2*nb_size+1,1,numel(index_list));
spot_x_vec = NaN(1,numel(index_list));
spot_y_vec = NaN(1,numel(index_list));
spot_time_vec = NaN(1,numel(index_list));
spot_fluo_vec = NaN(1,numel(index_list));
spot_frame_vec = NaN(1,numel(index_list));
spot_particle_vec = NaN(1,numel(index_list));
% get unique list of set-frame combinations
set_frame_array = unique([set_ref_vec' frame_ref_vec'],'row','stable');
iter = 1;
for i = 1:size(set_frame_array,1)    
    tic
    setID = set_frame_array(i,1);
    frame = set_frame_array(i,2);  
    
    % get nucleus
    frame_set_filter = set_ref_vec==setID&frame_ref_vec==frame;
    nc_x = nc_x_ref_vec(frame_set_filter);
    nc_y = nc_y_ref_vec(frame_set_filter);    
  
    % particle positions    
    time_vec = time_ref_vec(frame_set_filter);
    fluo_vec = fluo_ref_vec(frame_set_filter);
    spot_x = spot_x_ref_vec(frame_set_filter);
    spot_y = spot_y_ref_vec(frame_set_filter);        
    spot_z = spot_z_ref_vec(frame_set_filter);  
    particle_ids = particle_ref_vec(frame_set_filter);
    src = set_key(set_key.setID==setID,:).prefix{1};
    % load stacks   
    mcp_stack = load_stacks(rawPath, src, frame, mcp_channel);
    protein_stack = load_stacks(rawPath, src, frame, proteinChannel);  
    for j = 1:numel(nc_x)
        xn = nc_x(j);
        yn = nc_y(j);
        zp = spot_z(j);
        xp = spot_x(j);
        yp = spot_y(j);
        % calculate relative spot positions
        dx = xp-xn;
        dy = yp-yn;
        xp_rel = dx+nb_size+1;
        yp_rel = dy+nb_size+1;  
        rad_rel = sqrt(dx^2+dy^2);
        if rad_rel <= nc_rad
            % take slice
            protein_slice = protein_stack(yn-nb_size:yn+nb_size,xn-nb_size:xn+nb_size,zp);
            mcp_slice = mcp_stack(yn-nb_size:yn+nb_size,xn-nb_size:xn+nb_size,zp);
            nc_protein_array(:,:,1,iter) = protein_slice;
            nc_mcp_array(:,:,1,iter) = mcp_slice;
            % store characteristics
            pt_sm = imgaussfilt(protein_slice,7);
            mf_protein_vec(iter) = pt_sm(nb_size+1,nb_size+1);
            spot_x_vec(iter) = xp_rel;
            spot_y_vec(iter) = yp_rel;
            spot_time_vec(iter) = time_vec(j);
            spot_fluo_vec(iter) = fluo_vec(j);
            spot_frame_vec(iter) = frame;
            spot_particle_vec(iter) = particle_ids(j);
        end
        % increment
        iter = iter + 1;
        if mod(iter,5) ==0
            disp(iter)
        end
    end    
end    
% store
training_struct = struct;
training_struct.nc_mask = nc_mask;
training_struct.nc_protein_array = nc_protein_array;
training_struct.nc_mcp_array = nc_mcp_array;
training_struct.mf_protein_vec = mf_protein_vec;
training_struct.spot_x_vec = spot_x_vec;
training_struct.spot_y_vec = spot_y_vec;
training_struct.spot_time_vec = spot_time_vec;
training_struct.spot_particle_vec = spot_particle_vec;
training_struct.spot_frame_vec = spot_frame_vec;
training_struct.spot_fluo_vec = spot_fluo_vec;
% save
save([dataPath 'spot_loc_train_set.mat'],'training_struct')