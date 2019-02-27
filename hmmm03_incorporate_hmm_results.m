% hmmm03_incorporate_hmm_results(project,wInf,KInf)
%
% DESCRIPTION
% Script to conduct HMM inference
%
% ARGUMENTS
% project: master ID variable 
%
% wInf: memory used for inference
%
% KInf: number of states used for inference
%
% OPTIONS
% dropboxFolder: Path to data folder where you wish to save
%                pipeline-generated data sets and figures. If this
%                var is not specified, output will be saved one level
%                above git repo in the folder structure

% OUTPUT: 

function hmm_input_output = hmmm03_incorporate_hmm_results(project,wInf,KInf,varargin)

close all
warning('off','all') %Shut off Warnings

tWindow = 50*60; % determines width of sliding window
%%%%% These options will remain fixed for now
clipped = 1; % if 0 use "full" trace with leading and trailing 0's
fluo_field = 1; % specify which fluo field to (1 or 3)
off_traces_flag = 0; % if 1 filter for only traces that are observed to turn off
                     % if 2 filter and back-align
clipped_ends = 0; % if one, remove final w time steps from traces 
dpBootstrap = 1;
dataPath = ['../dat/' project '/'];
alphaFrac = 1302 / 6000;

% basic protein extraction params
ROIRadius = .25; % radus (um) of region used to query and compare TF concentrations
distLim = .6;
ControlType = 'edge';
%%%%%%%%%%%%%%
for i = 1:numel(varargin)    
    if strcmpi(varargin{i},'dropboxFolder')
        dataPath = [varargin{i+1} '\ProcessedEnrichmentData\' project '/'];
    end
%     if ischar(varargin{i})
%         if ismember(varargin{i},{'nBoots','savio','K','minDp','tWindow','sampleSize','maxWorkers','alphaFrac','dpBootstrap'})
%             eval([varargin{i} '=varargin{i+1};']);
%         end
%     end
end

load([dataPath '/nucleus_struct_protein.mat'],'nucleus_struct_protein') % load data

Tres = nucleus_struct_protein(i).TresInterp; % Time Resolution
alpha = alphaFrac*wInf;
pixelSize = nucleus_struct_protein(1).PixelSize;
d_type = '';
if dpBootstrap
    d_type = '_dp';
end

% Set write path (inference results are now written to external directory)
hmm_suffix =  ['/hmm_inference/w' num2str(wInf) '_t' num2str(Tres)...
    '_alpha' num2str(round(alpha*10)) '_f' num2str(fluo_field) '_cl' num2str(clipped) ...
    '_no_ends' num2str(clipped_ends) '_off' num2str(off_traces_flag) '/K' num2str(KInf) ...
    '_tw' num2str(tWindow/60) d_type '_1/']; 

file_list = dir([dataPath hmm_suffix '*.mat']);

inference_results = [];
for i = 1:numel(file_list)
    load([dataPath hmm_suffix file_list(i).name])
    if numel(fieldnames(output)) > 3
        inference_results = [inference_results output];
    end
end

%% extract protein info
% Generate distance vector for filtering snip stacks
dist_vec = [nucleus_struct_protein.(['spot_' ControlType '_dist_vec'])];
time_vec = [nucleus_struct_protein.time];
% frame_vec = [nucleus_struct_protein.frames];
% ap_vec = [nucleus_struct_protein.ap_vector];
% fluo_vec = [nucleus_struct_protein.fluo];
mf_protein_vec = [nucleus_struct_protein.protein];
% make distance_filter
dist_filter = dist_vec*pixelSize >= distLim;
set_vec = [];
particle_vec = [];
snip_filter = [];
for i = 1:numel(nucleus_struct_protein)
    snip_frame_vec = nucleus_struct_protein(i).snip_frame_vec;  
    nc_frames = nucleus_struct_protein(i).frames;
    snip_filter = [snip_filter ismember(nc_frames,snip_frame_vec)];
    set_vec = [set_vec repelem(nucleus_struct_protein(i).setID,sum(ismember(nc_frames,snip_frame_vec)))];
    particle_vec = [particle_vec repelem(nucleus_struct_protein(i).ParticleID,sum(ismember(nc_frames,snip_frame_vec)))];
end  
snip_filter = snip_filter == 1;
% filter vectors
% dist_vec = dist_vec(snip_filter&dist_filter);
% fluo_vec = fluo_vec(snip_filter&dist_filter);
% frame_vec = frame_vec(snip_filter&dist_filter);
time_vec = time_vec(snip_filter&dist_filter);
mf_protein_vec = mf_protein_vec(snip_filter&dist_filter);

% Snip stacks
spot_protein_snips = cat(3,nucleus_struct_protein.spot_protein_snips);
null_protein_snips = cat(3,nucleus_struct_protein.([ControlType '_null_protein_snips']));
% spot_mcp_snips = cat(3,nucleus_struct_protein.spot_mcp_snips);
% null_mcp_snips = cat(3,nucleus_struct_protein.([ControlType '_null_mcp_snips']));

% Make r reference array
snip_size = size(spot_protein_snips,1);
[y_ref, x_ref] = meshgrid(1:snip_size,1:snip_size);
r_ref = sqrt((x_ref-ceil(snip_size/2)).^2 + (y_ref-ceil(snip_size/2)).^2)*pixelSize;


r_nan = ones(size(r_ref));
r_nan(r_ref>ROIRadius) = NaN;
r_ft = repmat(r_nan,1,1,size(spot_protein_snips(:,:,dist_filter(snip_filter)),3));
spot_protein_vec = reshape(nanmean(nanmean(r_ft.*spot_protein_snips(:,:,dist_filter(snip_filter)),1),2),1,[]);
null_protein_vec = reshape(nanmean(nanmean(r_ft.*null_protein_snips(:,:,dist_filter(snip_filter)),1),2),1,[]);
% rescale mf reference vec
mf_protein_vec = mf_protein_vec * nanmean(null_protein_vec) / nanmean(mf_protein_vec) ;
particle_vec = particle_vec(dist_filter(snip_filter));

%% now extract corresponding hmm traces
particle_index = unique(particle_vec);
hmm_particle_vec = [inference_results.particle_ids];
index_vec = [];
sub_index_vec = [];
for i = 1:numel(inference_results)
    index_vec = [index_vec repelem(i,numel(inference_results(i).fluo_data))];
    sub_index_vec = [sub_index_vec 1:numel(inference_results(i).fluo_data)];
end

hmm_input_output = [];
for i = 1:numel(particle_index)
    % take average soft decoded result for particle
    ParticleID = particle_index(i);
    particle_ft = hmm_particle_vec==ParticleID;
    indices = index_vec(particle_ft);
    sub_indices = sub_index_vec(particle_ft);
    if ~isempty(sub_indices)
        temp = struct;
        temp.time = inference_results(indices(1)).time_data{sub_indices(1)};
        temp.fluo = inference_results(indices(1)).fluo_data{sub_indices(1)};
        z_mat = zeros(numel(temp.time),KInf);
        r_mat = zeros(numel(temp.time),KInf);
        zz_mat = zeros(KInf,KInf,numel(temp.time)-1);
    else
        continue;
    end
    for j = 1:numel(indices)
        [r,ri] = sort(inference_results(indices(j)).r);
        z = exp(inference_results(indices(j)).soft_struct.p_z_log_soft{sub_indices(j)});
        zz = exp(inference_results(indices(j)).soft_struct.p_zz_log_soft{sub_indices(j)});
        z_mat = z_mat + z(ri,:)';
        zz_mat = zz_mat + zz(ri,ri,:);
        r_mat = r_mat + z(ri,:)'.*r';
    end
    temp.z_mat = z_mat / numel(indices);
    temp.r_mat = r_mat / numel(indices);
    temp.zz_mat = zz_mat / numel(indices);
    % interpolate protein info to match frames of hmm
    sp_pt = spot_protein_vec(particle_vec==ParticleID);
    nn_pt = null_protein_vec(particle_vec==ParticleID);
    mf_pt = mf_protein_vec(particle_vec==ParticleID);
    tt_pt = time_vec(particle_vec==ParticleID);
    
    temp.spot_protein = interp1(tt_pt(~isnan(sp_pt)),sp_pt(~isnan(sp_pt)),temp.time);
    temp.null_protein = interp1(tt_pt(~isnan(sp_pt)),nn_pt(~isnan(sp_pt)),temp.time);
    temp.mf_protein = interp1(tt_pt(~isnan(sp_pt)),mf_pt(~isnan(sp_pt)),temp.time);

    % record general info for later use
    temp.ParticleID = ParticleID;
    temp.ROIRadius = ROIRadius;
    
    hmm_input_output = [hmm_input_output temp];
end
% save results
save([dataPath 'hmm_input_output_w' num2str(wInf) '_K' num2str(KInf) '.mat'],'hmm_input_output')