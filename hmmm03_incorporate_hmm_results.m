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
dataRoot = ['../dat/'];
alphaFrac = 1302 / 6000;

% basic protein extraction params
ControlType = 'edge';
controlProject = '';
%%%%%%%%%%%%%%
for i = 1:numel(varargin)    
    if strcmpi(varargin{i},'dropboxFolder')
        dataRoot = [varargin{i+1} '\ProcessedEnrichmentData\'];
    end
    if ischar(varargin{i})
        if ismember(varargin{i},{'dpBootstrap','controlProject'})
            eval([varargin{i} '=varargin{i+1};']);
        end
    end
end
d_type = '';
if dpBootstrap
    d_type = '_dp';
end
project_cell = {project};
id_cell = {'target'};
if ~isempty(controlProject)
    project_cell{2} = controlProject;
    id_cell{2} = 'control';
end
master_struct = struct;

for num = 1:numel(project_cell)
    prj = project_cell{num};
    id = id_cell{num};
    master_struct(num).project = prj;
    master_struct(num).ID = id;
    
    dataPath = [dataRoot prj '/'];
    load([dataPath '/nucleus_struct_protein.mat'],'nucleus_struct_protein') % load data
    % check for necessary fields
    analysis_fields = {'TresInterp','fluo_interp','time_interp'};
    if ~isfield(nucleus_struct_protein,analysis_fields{1})
        warning('Interpolation fields missing. Adding now.')
        load([dataPath '/nucleus_struct.mat'],'nucleus_struct') 
        for i = 1:numel(nucleus_struct)
            for a = 1:numel(analysis_fields)
                nucleus_struct_protein(i).(analysis_fields{a}) = nucleus_struct(i).(analysis_fields{a});
            end
        end
    end
    Tres = nucleus_struct_protein(i).TresInterp; % Time Resolution
    alpha = alphaFrac*wInf;

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
    particle_vec = [];
    for i = 1:numel(nucleus_struct_protein)
        particle_vec = [particle_vec repelem(nucleus_struct(i).ParticleID,numel(nucleus_struct_protein(i).time))];
    end
    time_vec = [nucleus_struct_protein.time];
    fluo_vec = [nucleus_struct_protein.fluo];
    spot_protein_vec = [nucleus_struct_protein.spot_protein_vec];
    spot_mcp_vec = [nucleus_struct_protein.spot_mcp_vec];
    null_protein_vec = [nucleus_struct_protein.([ControlType '_null_protein_vec'])];
    mf_protein_vec = [nucleus_struct_protein.mf_null_protein_vec];
    mf_count_vec = [nucleus_struct_protein.mf_null_px_counts];
    orig_protein_vec = [nucleus_struct_protein.protein];
    % rescale mf reference vec

    %%% now extract corresponding hmm traces
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
        temp = struct;
        % extract relevant vectors
        ff_pt = fluo_vec(particle_vec==ParticleID);
        mcp_pt = spot_mcp_vec(particle_vec==ParticleID);
        sp_pt = spot_protein_vec(particle_vec==ParticleID);
        nn_pt = null_protein_vec(particle_vec==ParticleID);
        mf_pt = mf_protein_vec(particle_vec==ParticleID);
        or_pt = orig_protein_vec(particle_vec==ParticleID);
        mf_ct = mf_count_vec(particle_vec==ParticleID);
        tt_pt = time_vec(particle_vec==ParticleID);
        if ~isempty(sub_indices) && sum(~isnan(mf_pt)) > 1 && sum(~isnan(nn_pt)) > 1        
            temp.time = inference_results(indices(1)).time_data{sub_indices(1)};
            temp.fluo = inference_results(indices(1)).fluo_data{sub_indices(1)};
            z_mat = zeros(numel(temp.time),KInf);
            r_mat = zeros(numel(temp.time),KInf);
            r_inf = zeros(1,KInf);
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
            r_inf = r_inf + r';
        end
        temp.z_mat = z_mat / numel(indices);
        temp.r_mat = r_mat / numel(indices);
        temp.zz_mat = zz_mat / numel(indices);
        temp.r = r_inf / numel(indices); 
        
        % checks using mcp channel to ensure that we are correctly matching
        % particles and time frames
        temp.mcp_check = interp1(tt_pt(~isnan(mcp_pt)),mcp_pt(~isnan(mcp_pt)),temp.time);
        temp.fluo_check = interp1(tt_pt(~isnan(ff_pt)),ff_pt(~isnan(ff_pt)),temp.time);
        temp.mf_counts = interp1(tt_pt(~isnan(mf_pt)),mf_ct(~isnan(mf_pt)),temp.time);
        % protein information
        temp.spot_protein = interp1(tt_pt(~isnan(sp_pt)),sp_pt(~isnan(sp_pt)),temp.time);
        temp.null_protein = interp1(tt_pt(~isnan(nn_pt)),nn_pt(~isnan(nn_pt)),temp.time);
        temp.mf_protein = interp1(tt_pt(~isnan(mf_pt)),mf_pt(~isnan(mf_pt)),temp.time);
        temp.orig_protein = interp1(tt_pt(~isnan(or_pt)),or_pt(~isnan(or_pt)),temp.time);
        % reset values to NaN that are too far removed from true ref point
        input_times = tt_pt(~isnan(sp_pt)&~isnan(nn_pt));
        dt_vec = NaN(size(temp.time));
        for t = 1:numel(dt_vec)
            dt_vec(t) = min(abs(temp.time(t)-input_times));
        end
        dt_filter = dt_vec > 60;
        temp.mcp_check(dt_filter) = NaN;
        temp.fluo_check(dt_filter) = NaN;
        temp.fluo(dt_filter) = NaN;
        temp.mf_counts(dt_filter) = NaN;
        temp.spot_protein_all = temp.spot_protein;
        temp.spot_protein(dt_filter) = NaN;
        temp.null_protein_all = temp.null_protein;
        temp.null_protein(dt_filter) = NaN;
        temp.mf_protein_all = temp.mf_protein;
        temp.mf_protein(dt_filter) = NaN;
        % record general info for later use
        temp.ParticleID = ParticleID; 

        hmm_input_output = [hmm_input_output temp];
    end
    master_struct(num).hmm_input_output = hmm_input_output;
end
% save results
save([dataPath 'hmm_input_output_w' num2str(wInf) '_K' num2str(KInf) '.mat'],'master_struct')