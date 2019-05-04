% Script to calculate average spot profile from data to use for protein
% sampling
clear
close all

project = 'Dl-Ven x snaBAC';
rawPath = 'E:\LocalEnrichment\Data\PreProcessedData\';
dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
dataPath = [dropboxFolder 'ProcessedEnrichmentData\' project '/'];

% sampling parameters
n_spots = 1000;
mcp_channel = 2;
snip_size = 25;
stack_size = 5;
% load nucleus structure
load([dataPath 'nucleus_struct.mat'])
load([dataPath '/set_key.mat'],'set_key')
% remove all frames that do not contain a segmented particle or that
% contain a particle that fails QC standards
% nucleus_struct = nucleus_struct([nucleus_struct.qc_flag]==1);
clean_fields = {'xPos','yPos','xPosParticle','yPosParticle','zPosParticle','fluo','time','frames'};
for i = 1:numel(nucleus_struct)
    fluo = nucleus_struct(i).fluo;
    nan_ft = ~isnan(fluo);
    for j = 1:numel(clean_fields)
        vec = nucleus_struct(i).(clean_fields{j});
        nucleus_struct(i).(clean_fields{j}) = vec(nan_ft);
    end
end

set_ref = [];
frame_ref = [];
for i = 1:numel(nucleus_struct)
    frame_ref = [frame_ref nucleus_struct(i).frames];
    set_ref = [set_ref repelem(nucleus_struct(i).setID,numel(nucleus_struct(i).frames))];
end
spot_x_ref = [nucleus_struct.xPosParticle];
spot_y_ref = [nucleus_struct.yPosParticle];
spot_z_ref = [nucleus_struct.zPosParticle];

set_frame_array = unique([set_ref' frame_ref'],'row');
set_index = unique(set_ref);
%%% make source key
src_cell = {1,numel(set_index)};
pt_id_vec = [nucleus_struct.ParticleID];
for s = 1:numel(set_index)
    ind = find([nucleus_struct.setID]==set_index(s),1);   
    src = nucleus_struct(ind).source_path;    
    src_cell{s} = src;     
end
% get dim info
xDim = nucleus_struct(1).xDim;
yDim = nucleus_struct(1).yDim;
zDim = nucleus_struct(1).zDim;

psf_cell = cell(1,n_spots);
n_sampled = 0;
rng(123); % for replicability
sampling_order = randsample(1:size(set_frame_array,1),size(set_frame_array,1),false);
% sample spots
for i = sampling_order    
    tic
    exit_flag = 0;
    setID = set_frame_array(i,1);
    frame = set_frame_array(i,2); 

    frame_set_filter = set_ref==setID&frame_ref==frame;
  
    % particle positions        
    spot_x_vec = spot_x_ref(frame_set_filter);
    spot_y_vec = spot_y_ref(frame_set_filter);        
    spot_z_vec = spot_z_ref(frame_set_filter);      
    src = set_key(set_key.setID==setID,:).prefix{1};
    
    % load stacks
    mcp_stack = load_stacks(rawPath, src, frame, mcp_channel);
    if isempty(mcp_stack)
        continue
    end
    for j = 1:numel(spot_x_vec)
        x_spot = spot_x_vec(j);
        y_spot = spot_y_vec(j);
        z_spot = spot_z_vec(j);
        spot_stack = NaN(2*snip_size + 1,2*snip_size + 1,2*stack_size+1);
        % pull sample
        x_range = max(1,x_spot-snip_size):min(xDim,x_spot+snip_size);
        y_range = max(1,y_spot-snip_size):min(yDim,y_spot+snip_size);
        z_range = max(1,z_spot-stack_size):min(zDim,z_spot+stack_size);
        x_range_full = x_spot-snip_size:x_spot+snip_size;
        y_range_full = y_spot-snip_size:y_spot+snip_size; 
        z_range_full = z_spot-stack_size:z_spot+stack_size; 
        
        spot_stack(ismember(y_range_full,y_range),ismember(x_range_full,x_range),ismember(z_range_full,z_range)) = ...
            mcp_stack(y_range,x_range,z_range);
        
        psf_cell{j} = spot_stack;
        n_sampled = n_sampled + 1;
        
        if n_sampled > n_spots
            exit_flag = 1;
            break
        end
    end
    if exit_flag == 1
        break;
    end
end
%%
% find average psf 
mean_spot = nanmean(cat(4,psf_cell{:}),4);
mean_spot(isnan(mean_spot)) = 0;
mean_psf = mean_spot / nansum(mean_spot(:));


% fit gaussian to mean spot
[mesh_y,mesh_x, mesh_z] = meshgrid(1:size(mean_psf,2), 1:size(mean_psf,1), 1:size(mean_psf, 3));
mesh_x = mesh_x - snip_size - 1;
mesh_y = mesh_y - snip_size - 1;
mesh_z = mesh_z - stack_size - 1;

% Single 3D generalized gaussian function
single3DGaussian = @(params) params(3)*...params(1)^2 * params(2) / (2*pi)^1.5 * ...
    exp(-.5*( params(1)*(mesh_x).^2 + params(1)*(mesh_y).^2 + params(2)*(mesh_z).^2 )) + params(4);

objective_fun = @(params) single3DGaussian(params) - (mean_spot);

centroid_guess = [size(mean_psf, 2)/2, size(mean_psf, 1)/2, size(mean_psf,3)/2];

initial_parameters = [1/2,1/2,1.2,1];
ub = [1 1 Inf Inf];%f snip_size snip_size snip_size];
lb = [0 0 0 0];% -snip_size -snip_size -snip_size];
%%% params and fits: %%%
%fitting options
lsqOptions=optimset('maxfunevals',10000,'TolFun',1e-4,'maxiter',10000);

fit = lsqnonlin(objective_fun, initial_parameters,lb,ub,lsqOptions);


% save dim params
psf_dims = struct;
psf_dims.xy_sigma = fit(1)^-.5;
psf_dims.z_sigma = fit(2)^-.5;

save([dataPath 'psf_dims.mat'],'psf_dims')
