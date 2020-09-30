% Script to call primary cpHMM wrapper function

clear
close all
warning('off','all') %Shut off Warnings

% get inference options
inferenceOptions = determineInferenceOptions;

% project identifier
projectName = '2xDl-Ven_hbP2P-mCh';
% default path to model scripts
modelPath = './utilities';

 
%%%%%%%%%%%%%%



if inferenceOptions.ProteinBinFlag && inferenceOptions.savioFlag
    nBoots = 1; % will run multiple instances on savio
else  
    nBoots = 5;
end

% DataPath = ['S:\Nick\Dropbox\ProcessedEnrichmentData\' project '\'];%
DataPath = ['../../dat/tf_enrichment/' projectName '/'];

addpath(modelPath); % Route to utilities folder
% check that we have proper fields
if ~inferenceOptions.dpBootstrap
    warning('Bootstrap option not selected. Setting nBoots to 1')
    nBoots = 1;
end

%-------------------------------System Vars-------------------------------%
if contains(projectName,'hbP2P')
    alphaFrac = 1275 / 4670;
elseif contains(projectName,'snaBAC')
    alphaFrac = 1302 / 6444;
end
alpha = alphaFrac*inferenceOptions.nSteps;

%----------------------------Set Write Paths------------------------------%
d_type = '';
if inferenceOptions.dpBootstrap
    d_type = '_dp';
end
% Set write path (inference results are now written to external directory)
if inferenceOptions.fluo3DFlag
    fluo_suffix = 'f3D';
else
    fluo_suffix = 'f2D';
end
% if protein_bin_flag    
%     load([DataPath '/nucleus_struct_protein.mat'],'nucleus_struct_protein') % load data
%     analysis_struct = nucleus_struct_protein;
%     clear nucleus_struct_protein
if inferenceOptions.ProteinBinFlag
    out_suffix =  ['/hmm_inference_protein/w' num2str(inferenceOptions.nSteps) '_inferenceOptions.nStates' num2str(inferenceOptions.nStates) '_' fluo_suffix '/']; 
else
  out_suffix =  ['/hmm_inference_mf/w' num2str(inferenceOptions.nSteps) '_inferenceOptions.nStates' num2str(inferenceOptions.nStates) '_' fluo_suffix '/'];   
end
load([DataPath '/nucleus_struct.mat'],'nucleus_struct') % load data
analysis_struct = nucleus_struct;
clear nucleus_struct
% end

% set write path
out_prefix = ['/global/scratch/nlammers/' projectName '/']; %hmmm_data/inference_out/';
% outDir = [DataPath out_suffix];%
outDir = [out_prefix out_suffix];
mkdir(outDir);

Tres = analysis_struct(1).TresInterp; % Time Resolution
% filter for quality traces of sufficient length
trace_struct_filtered = [];
for i = 1:length(analysis_struct)
    temp = struct;
    if inferenceOptions.fluo3DFlag
        fluo = analysis_struct(i).fluo3D_interp;
    else
        fluo = analysis_struct(i).fluo_interp;
    end
    time = analysis_struct(i).time_interp;    
    if sum(~isnan(fluo)) >= 2 % NL: this is a bit redundant. Leaving for now
        temp.fluo = fluo;
        temp.time = time;
        if inferenceOptions.ProteinBinFlag
            temp.mf_protein = nanmean(analysis_struct(i).raw_nc_protein);
        end
        temp.qc_flag = analysis_struct(i).qc_flag;
        temp.ParticleID = analysis_struct(i).ParticleID;  
        temp.N = numel(fluo);
        trace_struct_filtered = [trace_struct_filtered temp];
    end
end
trace_struct_filtered = trace_struct_filtered([trace_struct_filtered.qc_flag]==1);

if inferenceOptions.ProteinBinFlag 
    % estimate number of bins 
    if inferenceOptions.automaticBinning
        nTotal = sum([trace_struct_filtered.N]);
        n_protein_bins = ceil(nTotal/inferenceOptions.SampleSize);
    end
    % generate list of average protein levels
    mf_list = [trace_struct_filtered.mf_protein];
    % generate protein groupings    
    q_vec = linspace(0,1,n_protein_bins+1);        
    mf_prctile_vec = quantile(mf_list,q_vec);    
    % assign traces to groups    
    id_vec = discretize(mf_list,mf_prctile_vec);
    for i = 1:numel(trace_struct_filtered)
        trace_struct_filtered(i).mf_protein_bin = id_vec(i);
    end
end

% define iteration wrapper
iter_list = 1;
iter_ref_index = ones(size(trace_struct_filtered));
if inferenceOptions.ProteinBinFlag
    iter_list = 1:numel(mf_prctile_vec)-1;
    iter_ref_index = [trace_struct_filtered.mf_protein_bin];
end

%%% Conduct Inference
rng('shuffle'); % ensure we don't repeat 
% iterate through designated groups
for t = 1:length(iter_list)
    iter_filter = iter_ref_index == t;%ismember(iter_ref_index,[t-1, t, t+1]);
    for b = 1:nBoots
        iter_start = now;
        local_struct = struct;    
        output = struct;
        
        % Use current time as unique inference identifier 
        inference_id = num2str(round(10e5*now));

        % Generate filenames            
        fName_sub = ['hmm_results_t' inference_id];                
        out_file = [outDir '/' fName_sub];  
        
        % Extract fluo_data        
        inference_set = trace_struct_filtered(iter_filter);        
        skip_flag = 0;
        set_size = length([inference_set.fluo]);                 
        if isempty(inference_set)
            skip_flag = 1;
        elseif set_size < inferenceOptions.minDPperInf                    
            skip_flag = 1;                    
        end
        if skip_flag
            warning('Too few data points. Skipping')                
        else 
            sample_index = 1:length(inference_set);
            if inferenceOptions.dpBootstrap                        
                ndp = 0;    
                sample_ids = [];                    
                %Reset bootstrap size to be on order of set size for small bins
                if set_size < inferenceOptions.SampleSize
                    inferenceOptions.SampleSize = ceil(set_size/100)*100;
                end
                while ndp < inferenceOptions.SampleSize
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
            param_init = initialize_random (inferenceOptions.nStates, inferenceOptions.nSteps, fluo_data);
            % Approximate inference assuming iid data for param initialization                
            local_iid_out = local_em_iid_reduced_memory(fluo_data, param_init.v, ...
                                param_init.noise, inferenceOptions.nStates, inferenceOptions.nSteps, alpha, 500, 1e-4);
            noise_iid = 1/sqrt(exp(local_iid_out.lambda_log));
            v_iid = exp(local_iid_out.v_logs);            
            p = gcp('nocreate');
            if isempty(p)
                parpool(inferenceOptions.maxWorkers); %6 is the number of cores the Garcia lab server can reasonably handle per user.
            elseif p.NumWorkers > inferenceOptions.maxWorkers
                delete(gcp('nocreate')); % if pool with too many workers, delete and restart
                parpool(inferenceOptions.maxWorkers);
            end
            parfor i_local = 1:inferenceOptions.n_localEM % Parallel Local EM                
                % Random initialization of model parameters
                param_init = initialize_random_with_priors(inferenceOptions.nStates, noise_iid, v_iid);
                % Get Intial Values
                pi0_log_init = log(param_init.pi0);
                A_log_init = log(param_init.A);
                v_init = param_init.v;                        
                noise_init = param_init.noise;
                %--------------------LocalEM Call-------------------------%
                local_out = local_em_MS2_reduced_memory(fluo_data, ...
                    v_init, noise_init, pi0_log_init', A_log_init, inferenceOptions.nStates, inferenceOptions.nSteps, ...
                    alpha, inferenceOptions.nStinferenceOptions.epsMax, inferenceOptions.eps);                    
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
                local_struct(i_local).total_stinferenceOptions.eps = local_out.n_iter;               
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
            output.total_stinferenceOptions.eps = local_struct(max_index).total_stinferenceOptions.eps;                                  
            output.total_time = 100000*(now - iter_start);            
            % other inference characteristics            
            output.protein_bin_flag = inferenceOptions.ProteinBinFlag;
            if inferenceOptions.ProteinBinFlag
                output.protein_bin = t;
                output.protein_bin_list = iter_list;
                output.protein_bin_edges = mf_prctile_vec;
            end
            output.dp_bootstrap_flag = inferenceOptions.dpBootstrap;   
            output.iter_id = b;                        
            output.particle_ids = sample_particles;
            if inferenceOptions.dpBootstrap                                    
                output.N = ndp;
            end
            output.w = inferenceOptions.nSteps;
            output.alpha = alpha;
            output.deltaT = Tres;
            output.inferenceOptions.fluo3DFlag = inferenceOptions.fluo3DFlag;
            output.sampleSize = inferenceOptions.SampleSize; 
            % save inference data used
            output.fluo_data = fluo_data;
            output.time_data = time_data;
        end
        output.skip_flag = skip_flag;
%         disp('saving...')
        save([out_file '.mat'], 'output');           
    end  
end
 
