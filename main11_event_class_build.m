% Script to generate training set for event prediction
clear
close all
addpath('utilities')
% define core ID variables
project = 'Dl-Ven_snaBAC-mCh';
% project = 'Dl-Ven_hbP2P-mCh';
dropboxFolder =  'E:\Nick\LivemRNA\Dropbox\';
dataPath = [dropboxFolder 'ProcessedEnrichmentData\' project '/'];
rawPath = 'E:\LocalEnrichment\Data\PreProcessedData\';
proteinChannel = 1;
% load input-output data set
K = 3;
w = 7;
dT = 20;
load([dataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'],'hmm_input_output')
load([dataPath 'nucleus_struct_protein.mat'],'nucleus_struct_protein')
load([dataPath '/set_key.mat'],'set_key')
%%% get image dimensions
px_size = nucleus_struct_protein(1).PixelSize;
z_size = nucleus_struct_protein(1).zStep;
xDim = nucleus_struct_protein(1).xDim;
yDim = nucleus_struct_protein(1).yDim;
zDim = nucleus_struct_protein(1).zDim;
xy_rad = round(.75 / px_size);
z_rad = round(.5 / z_size);
% specify indexing vectors
particle_vec_hmm = [hmm_input_output.ParticleID];
time_vec_raw = [nucleus_struct_protein.time];
spot_x_ref = [nucleus_struct_protein.xPosParticle];
spot_y_ref = [nucleus_struct_protein.yPosParticle];
ctrl_x_ref = [nucleus_struct_protein.edge_null_x_vec];
ctrl_y_ref = [nucleus_struct_protein.edge_null_y_vec];
spot_z_ref = [nucleus_struct_protein.zPosParticle];
particle_vec_raw = NaN(size(time_vec_raw));
set_vec_raw = NaN(size(time_vec_raw));
frame_vec_raw = [nucleus_struct_protein.frames];
% initialize class vectors
event_class_vec = NaN(size(time_vec_raw));
event_size_vec = NaN(size(time_vec_raw));
prev_event_dur_vec = NaN(size(time_vec_raw));
burst_class_vec = NaN(size(time_vec_raw));
mf_protein_vec =  NaN(size(time_vec_raw));
iter = 1;
for i = 1:numel(nucleus_struct_protein)
    fluo = nucleus_struct_protein(i).fluo;
    pID = nucleus_struct_protein(i).ParticleID;
    setID = nucleus_struct_protein(i).setID;
    particle_vec_raw(iter:iter+numel(fluo)-1) = pID;
    set_vec_raw(iter:iter+numel(fluo)-1) = setID;
    iter = iter + numel(fluo);
end

% iterate through input/output data set and document events
for i = 1:numel(hmm_input_output)
    % extract HMM data
    change_points = hmm_input_output(i).change_points;
    z_diff_vec = hmm_input_output(i).z_diff_vec;
    z_vec = hmm_input_output(i).z_vec;
    prev_dur_vec = hmm_input_output(i).z_dur_lead_vec;
    next_dur_vec = hmm_input_output(i).z_dur_lag_vec;
    next_sz_vec = hmm_input_output(i).sz_lag_vec;
    % extend burst start labels
    zd_padded = [NaN(1,4) z_diff_vec NaN(1,4)];
    z_shift_array = vertcat(zd_padded,circshift(zd_padded,1));
    z_diff_new = nanmax(z_shift_array(:,5:end-4));
    time = hmm_input_output(i).time;
    mf_protein = hmm_input_output(i).mf_protein;    
    ptID = hmm_input_output(i).ParticleID;
    % extract raw protein time data
    pt_ids = find(particle_vec_raw==ptID);
    for id = pt_ids
        pt_time = time_vec_raw(id);
        [minTime, mi] = min(abs(pt_time-time));
        if minTime <= dT
            prev_event_dur_vec(id) = prev_dur_vec(mi);
            event_class_vec(id) = z_diff_new(mi);
            burst_class_vec(id) = z_vec(mi);
            mf_protein_vec(id) = mf_protein(mi);
            event_size_vec(id) = next_dur_vec(mi)*next_sz_vec(mi);
        end
    end
end

% remove entries for which either class or data are Missing
nan_ft = ~isnan(event_class_vec)...
    &spot_z_ref>z_rad&spot_z_ref<=(zDim-z_rad)...
    &spot_x_ref>xy_rad&spot_x_ref<=(xDim-xy_rad)...
    &spot_y_ref>xy_rad&spot_y_ref<=(xDim-xy_rad)...
    &ctrl_x_ref>xy_rad&ctrl_x_ref<=(xDim-xy_rad)...
    &ctrl_y_ref>xy_rad&ctrl_y_ref<=(xDim-xy_rad);

% filter sets
spot_x_ref = spot_x_ref(nan_ft);
spot_y_ref = spot_y_ref(nan_ft);
spot_z_ref = spot_z_ref(nan_ft);
ctrl_x_ref = ctrl_x_ref(nan_ft);
ctrl_y_ref = ctrl_y_ref(nan_ft);
event_class_vec = event_class_vec(nan_ft);
burst_class_vec = burst_class_vec(nan_ft);
event_size_vec = event_size_vec(nan_ft);
mf_protein_vec = mf_protein_vec(nan_ft);
prev_event_dur_vec = prev_event_dur_vec(nan_ft);
time_vec_raw = time_vec_raw(nan_ft);
particle_vec_raw = particle_vec_raw(nan_ft);
set_vec_raw = set_vec_raw(nan_ft);
frame_vec_raw = frame_vec_raw(nan_ft);
%%% Sample protein
% get unique list of set-frame combinations
mcp_channel = double(proteinChannel==1) + 1;
set_frame_array = unique([set_vec_raw' frame_vec_raw'],'row');
% initialize arrays
xy_dim = 2*xy_rad+1;
z_dim = 2*z_rad+1;
locus_protein_stack = NaN(xy_dim,xy_dim,z_dim,sum(nan_ft));
locus_mcp_stack = NaN(xy_dim,xy_dim,z_dim,sum(nan_ft));
control_protein_stack = NaN(xy_dim,xy_dim,z_dim,sum(nan_ft));
iter = 1;
for i = 1:size(set_frame_array,1)    
    tic
    setID = set_frame_array(i,1);
    frame = set_frame_array(i,2);          
    % get nucleus
    frame_set_filter = set_vec_raw==setID&frame_vec_raw==frame;    
  
    % particle positions        
    spot_x_vec = spot_x_ref(frame_set_filter);
    spot_y_vec = spot_y_ref(frame_set_filter); 
    ctrl_x_vec = ctrl_x_ref(frame_set_filter);
    ctrl_y_vec = ctrl_y_ref(frame_set_filter); 
    spot_z_vec = spot_z_ref(frame_set_filter);  
    particle_id_vec = particle_vec_raw(frame_set_filter);
    src = set_key(set_key.setID==setID,:).prefix{1};
    % load stacks
    tic
    mcp_stack = double(load_stacks(rawPath, src, frame, mcp_channel));
    protein_stack = double(load_stacks(rawPath, src, frame, proteinChannel));
    toc
    % iterate through particles
    for j = 1:numel(spot_x_vec)        
        xp = spot_x_vec(j);
        yp = spot_y_vec(j);
        zp = spot_z_vec(j);
        locus_protein_stack(:,:,:,iter) = ...
            protein_stack(yp-xy_rad:yp+xy_rad,xp-xy_rad:xp+xy_rad,zp-z_rad:zp+z_rad);
        locus_mcp_stack(:,:,:,iter) = ...
            mcp_stack(yp-xy_rad:yp+xy_rad,xp-xy_rad:xp+xy_rad,zp-z_rad:zp+z_rad);
        xc = ctrl_x_vec(j);
        yc = ctrl_y_vec(j);
        control_protein_stack(:,:,:,iter) = ...
            protein_stack(yc-xy_rad:yc+xy_rad,xc-xy_rad:xc+xy_rad,zp-z_rad:zp+z_rad);
        % increment
        iter = iter + 1;        
    end     
    disp(i)
end

% store in structure
training_struct = struct;
training_struct.event_class_vec = event_class_vec;
training_struct.burst_class_vec = burst_class_vec;
training_struct.mf_protein_vec = mf_protein_vec;
training_struct.event_size_vec = event_size_vec;
training_struct.prev_event_dur_vec = prev_event_dur_vec;
training_struct.particle_vec_protein = particle_vec_raw;
training_struct.time_vec = time_vec_raw;
training_struct.locus_protein_stack = locus_protein_stack;
training_struct.control_protein_stack = control_protein_stack;
training_struct.locus_mcp_stack = locus_mcp_stack;
% save structure
save([dataPath 'event_training_data.mat'],'training_struct','-v7.3') 

