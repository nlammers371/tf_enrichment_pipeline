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
%
% controlProject: specifies a project to use as an external control
%
% OUTPUT: hmm_input_output, structure containing vectors of protein and MS2
% intensities, along with corresponding HMM-decoded activity trajectories

function hmm_input_output = hmmm03_incorporate_hmm_results(project,wInf,KInf,varargin)

close all

tWindow = 50*60; % determines width of sliding window
%%%%% These options will remain fixed for now
clipped = 1; % if 0 use "full" trace with leading and trailing 0's
fluo_field = 1; % specify which fluo field to (1 or 3)

dpBootstrap = 0;
dataPath = ['../dat/' project '/'];
alphaFrac = 1302 / 6000;

% basic protein extraction params
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
    '/K' num2str(KInf) '_tw' num2str(tWindow/60) d_type '_1/']; 

file_list = dir([dataPath hmm_suffix '*.mat']);

% compile list of inference results
inference_results = [];
for i = 1:numel(file_list)
    load([dataPath hmm_suffix file_list(i).name])
    if numel(fieldnames(output)) > 3
        inference_results = [inference_results output];
    end
end

%%% extract protein info

% Generate distance vector for filtering snip stacks
particle_vec = [nucleus_struct_protein.ParticleID];
time_vec = [nucleus_struct_protein.time];
fluo_vec = [nucleus_struct_protein.fluo];
spot_protein_vec = [nucleus_struct_protein.spot_protein_vec];
spot_mcp_vec = [nucleus_struct_protein.spot_protein_vec];
null_protein_vec = [nucleus_struct_protein.(['null_' ControlType '_protein_vec'])];
mf_protein_vec = [nucleus_struct_protein.null_mf_protein_vec];
mf_count_vec = [nucleus_struct_protein.mf_null_px_counts];
% rescale mf reference vec

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
    ff_pt = fluo(particle_vec==ParticleID);
    mcp_pt = spot_mcp_vec(particle_vec==ParticleID);
    sp_pt = spot_protein_vec(particle_vec==ParticleID);
    nn_pt = null_protein_vec(particle_vec==ParticleID);
    mf_pt = mf_protein_vec(particle_vec==ParticleID);
    mf_ct = mf_count_vec(particle_vec==ParticleID);
    tt_pt = time_vec(particle_vec==ParticleID);
    
    % checks using mcp channel to ensure that we are correctly matching
    % particles and time frames
    temp.mcp_check = interp1(tt_pt(~isnan(mcp_pt)),mcp_pt(~isnan(mcp_pt)),temp.time);
    temp.fluo_check = interp1(tt_pt(~isnan(ff_pt)),ff_pt(~isnan(ff_pt)),temp.time);
    temp.mf_counts = interp1(tt_pt(~isnan(ff_pt)),mc_ct(~isnan(mf_pt)),temp.time);
    % protein information
    temp.spot_protein = interp1(tt_pt(~isnan(sp_pt)),sp_pt(~isnan(sp_pt)),temp.time);
    temp.null_protein = interp1(tt_pt(~isnan(nn_pt)),nn_pt(~isnan(nn_pt)),temp.time);
    temp.mf_protein = interp1(tt_pt(~isnan(mf_pt)),mf_pt(~isnan(mf_pt)),temp.time);

    % record general info for later use
    temp.ParticleID = ParticleID;
    temp.ROIRadius = ROIRadius;
    
    hmm_input_output = [hmm_input_output temp];
end
% save results
save([dataPath 'hmm_input_output_w' num2str(wInf) '_K' num2str(KInf) '.mat'],'hmm_input_output')