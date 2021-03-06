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

function hmm_input_output = hmmm03_incorporate_hmm_results(project,w,KInf,varargin)

close all
% min_time = 8*60;
%%%%% These options will remain fixed for now
dpBootstrap = 0;
% dataRoot = ['../dat/'];
alphaFrac = 1302 / 6000;
dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
dataPath = [dropboxFolder 'ProcessedEnrichmentData\' project '/'];


%%%%%%%%%%%%%%
for i = 1:numel(varargin)    
    if strcmpi(varargin{i},'dropboxFolder')
        dataRoot = [varargin{i+1} 'ProcessedEnrichmentData\'];
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
% extract 1c variables 
% dataPath = [dataRoot project '/'];
load([dataPath '/nucleus_struct_protein.mat'],'nucleus_struct_protein') % load data

minDP = nucleus_struct_protein(1).minDP;
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
Tres = nucleus_struct_protein(1).TresInterp; % Time Resolution
alpha = alphaFrac*w;
% generate alpha kernel 
alpha_kernel = NaN(1,w);
for i = 1:w
    if i < alpha
        alpha_kernel(i) = ((i-1) / alpha  + .5 * 1/alpha )*Tres;
    elseif i > alpha && (i-1) < alpha
        alpha_kernel(i) = Tres*(1 - .5*(alpha-i+1)*(1-(i-1)/alpha));
    else
        alpha_kernel(i) = Tres;
    end
end

% Set write path (inference results are now written to external directory)
hmm_suffix =  ['hmm_inference/w' num2str(w) '_K' num2str(KInf) '/']; 
file_list = dir([dataPath hmm_suffix 'hmm_results*.mat']);
if numel(file_list) > 1
    warning('multiple inference files detected. Ignoring all but first')
end

inference_results = load([dataPath hmm_suffix file_list(1).name]);
inference_results = inference_results.output;
% check for soft fit structure
soft_fit_flag = 1;
if exist([dataPath hmm_suffix 'soft_fit_struct.mat']) > 0
    fit_props = dir([dataPath hmm_suffix 'soft_fit_struct.mat']);
    fit_date = datenum(fit_props(1).date);
    hmm_date = datenum(file_list(1).date);
    if fit_date > hmm_date
        soft_fit_flag = 0;
        load([dataPath hmm_suffix 'soft_fit_struct.mat']);
    end
end
qc_indices = find([nucleus_struct_protein.qc_flag]==1);
if soft_fit_flag
    disp('conducting single trace fits...')
    A_log = log(inference_results.A_mat);
    v = inference_results.r*Tres;
    sigma = sqrt(inference_results.noise);
    pi0_log = log(inference_results.pi0); 
    eps = 1e-4;

    
    fluo_values = cell(numel(qc_indices),1);
    rm_indices = [];
    for i = 1:numel(qc_indices)
        fluo = nucleus_struct_protein(qc_indices(i)).fluo_interp;
        if numel(fluo) < minDP
            error('problem with qc flag')
        end
        fluo_values{i} = fluo;
    end    
    fluo_values = fluo_values(~ismember(qc_indices,rm_indices));
    qc_indices = qc_indices(~ismember(qc_indices,rm_indices));
    tic 
    local_em_outputs = local_em_MS2_reduced_memory (fluo_values, ...
                  v', sigma, pi0_log, A_log, KInf, w, alpha, 1, eps);
    toc
    soft_fit_struct = local_em_outputs.soft_struct;
    save([dataPath hmm_suffix 'soft_fit_struct.mat'],'soft_fit_struct')
end

%%% extract protein info
particle_vec = [];
for i = 1:numel(nucleus_struct_protein)
    particle_vec = [particle_vec repelem(nucleus_struct_protein(i).ParticleID,numel(nucleus_struct_protein(i).time))];
end
time_vec = [nucleus_struct_protein.time];
fluo_vec = [nucleus_struct_protein.fluo];
spot_protein_vec = [nucleus_struct_protein.spot_protein_vec];
spot_mcp_vec = [nucleus_struct_protein.spot_mcp_vec];
null_protein_vec = [nucleus_struct_protein.edge_null_protein_vec];
mf_protein_vec = [nucleus_struct_protein.mf_null_protein_vec];
serial_protein_vec = [nucleus_struct_protein.serial_null_protein_vec];

%%% generate indexing vectors
particle_index = [nucleus_struct_protein.ParticleID];

%%% now extract corresponding hmm traces
hmm_input_output = [];
iter = 1;
for i = qc_indices%1:numel(particle_index)        
    % take average soft decoded result for particle
    ParticleID = particle_index(i);
    temp = struct;
    % extract relevant vectors from protein struct
    ff_pt = fluo_vec(particle_vec==ParticleID);
    mcp_pt = spot_mcp_vec(particle_vec==ParticleID);
    sp_pt = spot_protein_vec(particle_vec==ParticleID);
    nn_pt = null_protein_vec(particle_vec==ParticleID);
    mf_pt = mf_protein_vec(particle_vec==ParticleID);
    sr_pt = serial_protein_vec(particle_vec==ParticleID);       
    tt_pt = time_vec(particle_vec==ParticleID);
    % check to see if we ran inference for this one

    temp.time = nucleus_struct(i).time_interp;%inference_results(indices(1)).time_data{sub_indices(1)}; % doesn't matter which duplicate we reference
    temp.fluo = nucleus_struct(i).fluo_interp;%inference_results(indices(1)).fluo_data{sub_indices(1)};    
    % extract useful values
    [r,ri] = sort(inference_results.r);
    z = exp(soft_fit_struct.p_z_log_soft{iter});    
    temp.z_mat = z(ri,:)';    
    temp.r_mat = z(ri,:)'.*r';
    temp.r_inf = r';
    temp.r_vec = sum(temp.r_mat,2);
    [~,z_vec] = max(temp.z_mat,[],2);
    temp.z_vec = z_vec; 
    % make predicted fluo vec
    fluo_hmm = conv(temp.r_vec,alpha_kernel);
    temp.fluo_hmm = fluo_hmm(1:numel(temp.r_vec));        
    if  sum(~isnan(sr_pt)&~isnan(mf_pt)) > 2
        % checks using mcp channel to ensure that we are correctly matching
        % particles and time frames
        temp.mcp_check = interp1(tt_pt(~isnan(mcp_pt)),mcp_pt(~isnan(mcp_pt)),temp.time);
        temp.fluo_check = interp1(tt_pt(~isnan(ff_pt)),ff_pt(~isnan(ff_pt)),temp.time);
        % protein information
        temp.spot_protein = interp1(tt_pt(~isnan(sp_pt)),sp_pt(~isnan(sp_pt)),temp.time);        
        temp.mf_protein = interp1(tt_pt(~isnan(mf_pt)),mf_pt(~isnan(mf_pt)),temp.time);  
        temp.null_protein = interp1(tt_pt(~isnan(nn_pt)),nn_pt(~isnan(nn_pt)),temp.time);
        temp.serial_protein = interp1(tt_pt(~isnan(sr_pt)),sr_pt(~isnan(sr_pt)),temp.time);
        % reset values to NaN that are too far removed from true ref point
        input_times = tt_pt(~isnan(sp_pt)&~isnan(sr_pt));
        gap_times = tt_pt(~isnan(sp_pt)&isnan(sr_pt)); % look for assymmetries btw spot and control channels
        dt_vec_gap = zeros(size(temp.time));
        dt_vec_imbalance = zeros(size(temp.time));
        for t = 1:numel(dt_vec_gap)
            dt_vec_gap(t) = min(abs(temp.time(t)-input_times)); 
            if ~isempty(gap_times)
                dt_vec_imbalance(t) = min(abs(temp.time(t)-gap_times));
            end
        end
        temp.dt_filter_gap = dt_vec_gap > 60;           
        temp.dt_filter_imb = dt_vec_imbalance < Tres;
        % record general info for later use
        temp.ParticleID = ParticleID; 
        temp.Tres = Tres;
        hmm_input_output  = [hmm_input_output temp];
    end
    iter = iter + 1;
end
% find nearest neighbor particles
% generate array of average protein levels for each nucleus
time_vec = unique([hmm_input_output.time]);
time_vec = time_vec(~isnan(time_vec));
mf_array = NaN(numel(time_vec),numel(hmm_input_output));
start_time_vec = NaN(size(hmm_input_output));
stop_time_vec = NaN(size(hmm_input_output));
for i = 1:numel(hmm_input_output)
    t_vec = hmm_input_output(i).time;
    mf_vec = hmm_input_output(i).mf_protein;
    mf_array(ismember(time_vec,t_vec),i) = mf_vec;
    start_time_vec(i) = min(t_vec);
    stop_time_vec(i) = max(t_vec);
end
% now find closest match for each nucleus
for i = 1:numel(hmm_input_output)
    mf_i = mf_array(:,i);       
    dt_mf_vec = nanmean(abs(mf_array-mf_i));    
    dt_mf_vec(start_time_vec>start_time_vec(i)|stop_time_vec<stop_time_vec(i)) = NaN;
    dt_mf_vec(i) = NaN;
    [~, best_ind] = nanmin(dt_mf_vec);
    % record vales 
    time_swap = hmm_input_output(best_ind).time;  
    time_i = hmm_input_output(i).time; 
    % fill
    s_pt = NaN(size(time_i));
    s_pt(ismember(time_i,time_swap)) = hmm_input_output(best_ind).spot_protein(ismember(time_swap,time_i));
    mf_pt = NaN(size(time_i));
    mf_pt(ismember(time_i,time_swap)) = hmm_input_output(best_ind).mf_protein(ismember(time_swap,time_i));
    % record
    hmm_input_output(i).swap_ind = best_ind;
    hmm_input_output(i).swap_spot_protein = s_pt;
    hmm_input_output(i).swap_mf_protein = mf_pt;
end

% save results
save([dataPath 'hmm_input_output_w' num2str(w) '_K' num2str(KInf) '.mat'],'hmm_input_output')