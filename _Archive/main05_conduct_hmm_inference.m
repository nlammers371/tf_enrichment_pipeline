% main05_conduct_hmm_inference(project,varargin)
%
% DESCRIPTION
% Script to conduct HMM inference
%
% ARGUMENTS
% project: master ID variable 
%
% modelPath: file path to folder containing hmmm scripts
%
% w: Integer corresponding number of time steps for Pol II to transcribe
% gene
%
%
% OPTIONS
% dropboxFolder: Path to data folder where you wish to save
%                pipeline-generated data sets and figures. If this
%                var is not specified, output will be saved one level
%                above git repo in the folder structure
% savio: if 1, indicates we are running inference on savio cluster
% K: number of states
% minDp: min data points needed to be included in inferece
%
% OUTPUT: nucleus_struct_protein: compiled data set with protein samples

function output = main05_conduct_hmm_inference(project, DropboxFolder, varargin)

close all
warning('off','all') %Shut off Warnings
addpath('utilities');

% default path to model scripts
modelPath = './utilities';
savioFlag = 0; % set to 1 if using savio computing cluster
K = 3; % Number of states in system
w = 7; % "memory"--length of construct in time steps
ProteinBinFlag = false;
dpBootstrap = 1;
% minDp = 10;
if ProteinBinFlag
    nBoots = 2; % will run multiple instances on savio
else
    nBoots = 5;
end
n_protein_bins = 20;
sampleSize = 8000;
maxWorkers = 12;
alphaFrac = 1302 / 6000;
% default paths
[~ , DataPath, ~] =   header_function(DropboxFolder, project); 

%%%%% These options generally remain fixed 
clipped = 1; % if 0 use "full" trace with leading and trailing 0's
n_localEM = 25; % set num local runs
n_steps_max = 500; % set max steps per inference
eps = 1e-4; % set convergence criteria
min_dp_per_inf = 1000; % inference will be aborted if fewer present 
%%%%%%%%%%%%%%

for i = 1:numel(varargin)    
    if ischar(varargin{i}) && i ~= numel(varargin)        
        eval([varargin{i} '=varargin{i+1};']);        
    end
end

addpath(modelPath); % Route to utilities folder
% check that we have proper fields
if ~dpBootstrap
    warning('Bootstrap option not selected. Setting nBoots to 1')
    nBoots = 1;
end

%-------------------------------System Vars-------------------------------%
if contains(project,'hbP2P')
    alphaFrac = 1275 / 4670;
elseif contains(project,'snaBAC')
    alphaFrac = 1302 / 6444;
end
alpha = alphaFrac*w;

%----------------------------Bootstrap Vars-------------------------------%                                       
% max num workers
if savioFlag
    maxWorkers = 24;
end

%----------------------------Set Write Paths------------------------------%
d_type = '';
if dpBootstrap
    d_type = '_dp';
end

% Set write path (inference results are now written to external directory)
% out_suffix =  ['/hmm_inference_protein/w' num2str(w) '_K' num2str(K) '/']; 
if ProteinBinFlag
    load([DataPath '/nucleus_struct_protein.mat'],'nucleus_struct_protein') % load data
    analysis_struct = nucleus_struct_protein;
    out_suffix =  ['/hmm_inference_protein/w' num2str(w) '_K' num2str(K) '/']; 
else
    load([DataPath '/nucleus_struct.mat'],'nucleus_struct') % load data
    out_suffix =  ['/hmm_inference_mf/w' num2str(w) '_K' num2str(K) '/']; 
    analysis_struct = nucleus_struct;
end

% set write path
if savioFlag
    out_prefix = '/global/scratch/nlammers/'; %hmmm_data/inference_out/';
else    
    out_prefix = DataPath;
end
outDir = [out_prefix out_suffix];
mkdir(outDir);

Tres = analysis_struct(1).TresInterp; % Time Resolution
% filter for quality traces of sufficient length
trace_struct_filtered = [];
for i = 1:length(analysis_struct)
    temp = struct;
    time = analysis_struct(i).time_interp;
    fluo = analysis_struct(i).fluo_interp;
    if sum(~isnan(fluo)) >= 2 % NL: this is a bit redundant. Leaving for now
        temp.fluo = fluo;
        temp.time = time;
        if ProteinBinFlag
            temp.mf_protein = nanmean(analysis_struct(i).mf_null_protein_vec);
        end
        temp.qc_flag = analysis_struct(i).qc_flag;
        temp.ParticleID = analysis_struct(i).ParticleID;        
        trace_struct_filtered = [trace_struct_filtered temp];
    end
end
trace_struct_filtered = trace_struct_filtered([trace_struct_filtered.qc_flag]==1);

if ProteinBinFlag 
    % generate list of average protein levels
    mf_index = [trace_struct_filtered.mf_protein];
    % generate protein groupings    
    q_vec = linspace(0,1,n_protein_bins+1);        
    mf_prctile_vec = quantile(mf_index,q_vec);    
    % assign traces to groups    
    id_vec = discretize(mf_index,mf_prctile_vec);
    for i = 1:numel(trace_struct_filtered)
        trace_struct_filtered(i).mf_protein_bin = id_vec(i);
    end
end

% define iteration wrapper
iter_list = 1;
iter_ref_index = ones(size(trace_struct_filtered));
if ProteinBinFlag
    iter_list = 1:numel(mf_prctile_vec)-1;
    iter_ref_index = [trace_struct_filtered.mf_protein_bin];
end

%%% Conduct Inference

% iterate through designated groups
for t = 1:length(iter_list)
    iter_filter = iter_ref_index == t;
    for b = 1:nBoots
        iter_start = now;
        local_struct = struct;    
        output = struct;
        
        % Use current time as unique inference identifier 
        inference_id = num2str(round(10e5*now));

        % Generate filenames            
        fName_sub = ['hmm_results_t' inference_id '.mat'];                
        out_file = [outDir '/' fName_sub];  
        
        % Extract fluo_data        
        inference_set = trace_struct_filtered(iter_filter);        
        skip_flag = 0;
        set_size = length([inference_set.fluo]);                 
        if isempty(inference_set)
            skip_flag = 1;
        elseif set_size < min_dp_per_inf                    
            skip_flag = 1;                    
        end
        if skip_flag
            warning('Too few data points. Skipping')                
        else 
            sample_index = 1:length(inference_set);
            if dpBootstrap                        
                ndp = 0;    
                sample_ids = [];                    
                %Reset bootstrap size to be on order of set size for small bins
                if set_size < sampleSize
                    sampleSize = ceil(set_size/100)*100;
                end
                while ndp < sampleSize
                    tr_id = randsample(sample_index,1);
                    sample_ids = [sample_ids tr_id];
                    ndp = ndp + length(inference_set(tr_id).time);
                end
                fluo_data = cell([length(sample_ids), 1]);    
                time_data = cell([length(sample_ids), 1]);    
                sample_particles = [inference_set(sample_ids).ParticleID];
                for tr = 1:length(sample_ids)
                    fluo_data{tr} = inference_set(sample_ids(tr)).fluo;                    
                    time_data{tr} = inference_set(sample_ids(tr)).time;                    
                end            
            else % Take all relevant traces if not bootstrapping
                fluo_data = cell([length(inference_set), 1]);    
                time_data = cell([length(inference_set), 1]);    
                for tr = 1:length(inference_set)
                    fluo_data{tr} = inference_set(tr).fluo;
                    time_data{tr} = inference_set(tr).time;                    
                end
                sample_particles = [inference_set.ParticleID];
            end
            % Random initialization of model parameters
            param_init = initialize_random (K, w, fluo_data);
            % Approximate inference assuming iid data for param initialization                
            local_iid_out = local_em_iid_reduced_memory(fluo_data, param_init.v, ...
                                param_init.noise, K, w, alpha, 500, 1e-4);
            noise_iid = 1/sqrt(exp(local_iid_out.lambda_log));
            v_iid = exp(local_iid_out.v_logs);            
            p = gcp('nocreate');
            if isempty(p)
                parpool(maxWorkers); %6 is the number of cores the Garcia lab server can reasonably handle per user.
            elseif p.NumWorkers > maxWorkers
                delete(gcp('nocreate')); % if pool with too many workers, delete and restart
                parpool(maxWorkers);
            end
            parfor i_local = 1:n_localEM % Parallel Local EM                
                % Random initialization of model parameters
                param_init = initialize_random_with_priors(K, noise_iid, v_iid);
                % Get Intial Values
                pi0_log_init = log(param_init.pi0);
                A_log_init = log(param_init.A);
                v_init = param_init.v;                        
                noise_init = param_init.noise;
                %--------------------LocalEM Call-------------------------%
                local_out = local_em_MS2_reduced_memory(fluo_data, ...
                    v_init, noise_init, pi0_log_init', A_log_init, K, w, ...
                    alpha, n_steps_max, eps);                    
                %---------------------------------------------------------%                
                % Save Results 
                local_struct(i_local).inference_id = inference_id;
                local_struct(i_local).subset_id = i_local;
                local_struct(i_local).logL = local_out.logL;                
                local_struct(i_local).A = exp(local_out.A_log);
                local_struct(i_local).v = exp(local_out.v_logs).*local_out.v_signs;
                local_struct(i_local).r = exp(local_out.v_logs).*local_out.v_signs / Tres;                                
                local_struct(i_local).noise = 1/exp(local_out.lambda_log);
                local_struct(i_local).pi0 = exp(local_out.pi0_log);
                local_struct(i_local).total_steps = local_out.n_iter;               
                local_struct(i_local).soft_struct = local_out.soft_struct;               
            end
            [~, max_index] = max([local_struct.logL]); % Get index of best result                    
            % Save parameters from most likely local run
            output.pi0 =local_struct(max_index).pi0;                        
            output.r = local_struct(max_index).r(:);          
            output.noise = local_struct(max_index).noise;
            output.A = local_struct(max_index).A(:);
            output.A_mat = local_struct(max_index).A;            
            % get soft-decoded structure
            output.soft_struct = local_struct(max_index).soft_struct;
            % Info about run time
            output.total_steps = local_struct(max_index).total_steps;                                  
            output.total_time = 100000*(now - iter_start);            
            % other inference characteristics            
            output.protein_bin_flag = ProteinBinFlag;
            if ProteinBinFlag
                output.protein_bin = t;
                output.protein_bin_list = iter_list;
                output.protein_bin_edges = mf_prctile_vec;
            end
            output.dp_bootstrap_flag = dpBootstrap;   
            output.iter_id = b;            
            output.clipped = clipped;            
            output.particle_ids = sample_particles;
            if dpBootstrap                                    
                output.N = ndp;
            end
            output.w = w;
            output.alpha = alpha;
            output.deltaT = Tres; 
            output.sampleSize = sampleSize; 
            % save inference data used
            output.fluo_data = fluo_data;
            output.time_data = time_data;
        end
        output.skip_flag = skip_flag;
        save([out_file '.mat'], 'output');           
    end  
end
 
