clear
close all

addpath(genpath('./lib'))
addpath(genpath('../../utilities'))

% Load data
dataRoot = 'S:\Nick\Dropbox (Personal)\ProcessedEnrichmentData\';
if ~exist(dataRoot)
    dataRoot = 'C:\Users\nlamm\Dropbox (Personal)\ProcessedEnrichmentData\';
end
readPath = [dataRoot 'combinedOptoSets_v2' filesep];
load([readPath 'inference_data.mat'],'inference_data')

% set write path
writePath = [readPath 'cpHMM_results' filesep];
mkdir(writePath)

sampleSize = 2500;
nBoots = 10;
inferenceOptions.sampleSize = sampleSize;
inferenceOptions.nBoots = nBoots;
inferenceOptions.nStates = 2;
inferenceOptions.Tres = 20;
inferenceOptions.nSteps = 7;
inferenceOptions.alpha = 1.4;
inferenceOptions.eps = 1e-4;
inferenceOptions.maxWorkers = 25;
inferenceOptions.n_localEM = 25; % set num local runs
inferenceOptions.nStepsMax = 500; % set max stinferenceOptions.eps per inference

% group data by ap position
ap_bins = [-0.025 0.025];
ap_axis = ap_bins(1:end-1) + diff(ap_bins)/2;

ap_group_vec = discretize(inference_data.mean_ap_vec,ap_bins);
ap_index = unique(ap_group_vec);
ap_index = ap_index(~isnan(ap_index));
zero_flag_vec = inference_data.zero_flag_vec;
project_id_vec = inference_data.project_id_vec;
project_index_vec = unique(project_id_vec);

% get all possible combinations
elements = {ap_index [0 1] project_index_vec};
combCell = cell(1, numel(elements));
[combCell{:}] = ndgrid(elements{:});
combCell = cellfun(@(x) x(:), combCell,'uniformoutput',false); %there may be a better way to do this
project_ap_zero_array = [combCell{:}]; 

% calculate kniprs bins for each group
inference_struct = struct;
group_i = 1;
for i = 1:size(project_ap_zero_array)
    % get basic ID info
    ap_id = project_ap_zero_array(i,1);
    p_id = project_ap_zero_array(i,3);
    z_id = project_ap_zero_array(i,2);
    % filter for applicable traces
    trace_filter = find(ap_group_vec==ap_id & project_id_vec==p_id & (~zero_flag_vec | z_id));
    % get number of data points
    nDP = length([inference_data.fluo_vec_cell{trace_filter}]);
    % calculate how many bins this permits 
    nBins = floor(nDP/sampleSize);
    % get quantiles
    p_vec = linspace(0,1,nBins+1);
    kni_vec = inference_data.mean_knirps_vec(trace_filter);
    quant_vec = quantile(kni_vec,p_vec);
    kni_groups = discretize(kni_vec,quant_vec);
    % now divide up the traces
    inference_struct(i).fluo_data = cell(1,nBins);
    inference_struct(i).kni_data = cell(1,nBins);
    inference_struct(i).time_data = cell(1,nBins);
    inference_struct(i).particle_id_data = cell(1,nBins);
    inference_struct(i).particle_sub_id_data = cell(1,nBins);
    inference_struct(i).group_id_vec = NaN(1,nBins);
    for n = 1:nBins
        inference_struct(i).fluo_data{n} = inference_data.fluo_vec_cell(trace_filter(kni_groups==n));
        inference_struct(i).kni_data{n} = inference_data.knirps_vec_cell(trace_filter(kni_groups==n));
        inference_struct(i).time_data{n} = inference_data.time_vec_cell(trace_filter(kni_groups==n));
        inference_struct(i).particle_id_data{n} = inference_data.particle_id_vec(trace_filter(kni_groups==n));
        inference_struct(i).particle_sub_id_data{n} = inference_data.rep_id_vec(trace_filter(kni_groups==n));
        
        inference_struct(i).group_id_vec(n) = group_i;
        group_i = group_i + 1;
    end
    inference_struct(i).ap_id = ap_id;
    inference_struct(i).p_id = p_id;
    inference_struct(i).z_id = z_id;
    inference_struct(i).knirps_bins = quant_vec;
    inference_struct(i).ap_pos = ap_axis(ap_id);
end    

save([readPath 'inference_info.mat'],'inference_struct')

for i = 1:2%length(inference_struct)
        
    knirps_group_vec = 1:length(inference_struct(i).fluo_data);
    
    for k = 1:length(knirps_group_vec)
        group_i = inference_struct(i).group_id_vec(k);
        % iterate through bootstraps
        for b = 1:nBoots

            % Initialize structures to store results
            local_struct_temp = struct;  % stores results from individual iterations
            output = struct;  % structure that will be saved to file               

            % Extract subset of traces relevant to this subgroup       
            fluo_data_full = inference_struct(i).fluo_data{k};
            sample_index = 1:length(fluo_data_full);
        
            % take bootstrap sample             
            sample_ids = [];                                           
            ndp = 0;
            % randomly draw traces
            if length(sample_index) > 1
                while ndp < sampleSize
                    tr_id = randsample(sample_index,1);
                    sample_ids = [sample_ids tr_id];
                    ndp = ndp + length(fluo_data_full{tr_id});
                end
            else
                sample_ids = sample_index;
            end
            % add them to data cells
            fluo_data = fluo_data_full(sample_ids);    
            time_data = inference_struct(i).time_data{k}(sample_ids);
            kni_data = inference_struct(i).kni_data{k}(sample_ids);
            sample_particles = inference_struct(i).particle_id_data{k}(sample_ids);
            particle_sub_id_data = inference_struct(i).particle_sub_id_data{k}(sample_ids);
            % Random initialization of model parameters
            param_init = initialize_random(inferenceOptions.nStates, inferenceOptions.nSteps, fluo_data);

            % Approximate inference assuming iid data for param initialization                
            local_iid_out = local_em_iid_reduced_memory(fluo_data, param_init.v, ...
                                param_init.noise,inferenceOptions.nStates, inferenceOptions.nSteps, inferenceOptions.alpha, 500, 1e-4);

            noise_iid = 1/sqrt(exp(local_iid_out.lambda_log));
            v_iid = exp(local_iid_out.v_logs);  

            % create parallel pool if one does not already exist
            p = gcp('nocreate');
            if isempty(p)
                parpool(inferenceOptions.maxWorkers); %6 is the number of cores the Garcia lab server can reasonably handle per user.
            elseif p.NumWorkers > inferenceOptions.maxWorkers
                delete(gcp('nocreate')); % if pool with too many workers, delete and restart
                parpool(inferenceOptions.maxWorkers);
            end

            % conduct cpHMM inference
            parfor i_local = 1:inferenceOptions.n_localEM % Parallel Local EM 

                % Random initialization of model parameters
                param_init = initialize_random_with_priors(inferenceOptions.nStates, noise_iid, v_iid);

                % Get Intial Values
                pi0_log_init = log(param_init.pi0);
                A_log_init = log(param_init.A);
                v_init = param_init.v;                        
                noise_init = param_init.noise;

                %--------------------LocalEM Call-------------------------%
                local_out = local_em_MS2_reduced_memory_truncated(fluo_data, ...
                      v_init, noise_init, pi0_log_init', A_log_init, inferenceOptions.nStates, inferenceOptions.nSteps, ...
                  inferenceOptions.alpha, inferenceOptions.nStepsMax, inferenceOptions.eps);  
                
                %---------------------------------------------------------%                
                % Save Results                 
                local_struct_temp(i_local).subset_id = i_local;
                local_struct_temp(i_local).logL = local_out.logL;                
                local_struct_temp(i_local).A = exp(local_out.A_log);
                local_struct_temp(i_local).v = exp(local_out.v_logs).*local_out.v_signs;
                local_struct_temp(i_local).r = exp(local_out.v_logs).*local_out.v_signs / inferenceOptions.Tres;                                
                local_struct_temp(i_local).noise = 1/exp(local_out.lambda_log);
                local_struct_temp(i_local).pi0 = exp(local_out.pi0_log);
                local_struct_temp(i_local).total_time = local_out.n_iter;               
                local_struct_temp(i_local).soft_struct = local_out.soft_struct;               
            end

            % Record output
            [maxL, max_index] = max([local_struct_temp.logL]); % Get index of best result  

            % Save parameters from most likely local run
            output.pi0 =local_struct_temp(max_index).pi0;                        
            output.r = local_struct_temp(max_index).r(:);          
            output.noise = local_struct_temp(max_index).noise;
            output.A = local_struct_temp(max_index).A(:);
            output.A_mat = local_struct_temp(max_index).A;  
            output.max_logL = maxL;
            output.logL_results = [local_struct_temp.logL];

            % get soft-decoded structure
            output.soft_struct = local_struct_temp(max_index).soft_struct;                                                                     

            % other inference characteristics                        
            output.timeBin = 1;
            output.apBin = inference_struct(i).ap_id;
            output.additionalBin = k;
            output.groupID = group_i;                                        
            
            output.truncInference = 1;
            output.iter_id = b;                        
            output.particle_ids = sample_particles;    
            output.particle_sub_ids = particle_sub_id_data;    
            
            output.N = ndp;

            % save inference data used
            output.fluo_data = fluo_data;
            output.kni_mean = nanmean([kni_data{:}]);
            output.time_data = time_data;

            % Determine unique filename and sace

            % Generate filenames            
            fName_sub = ['hmm_results_group' sprintf('%03d',group_i) '_rep'];

            % generate random string
            rand_string = strrep(num2str(randsample(1:9,5,true)),' ','');
            % save
            out_file = [writePath fName_sub rand_string];          
            save([out_file '.mat'], 'output');           
        end        
    end
end