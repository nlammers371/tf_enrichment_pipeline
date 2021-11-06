function cpHMMInferenceGrouped(InputDataPath,OutputDataPath,modelSpecs,varargin)

  % Script to call primary cpHMM wrapper function

  % get inference options
  inferenceOptions = determineInferenceOptions(modelSpecs,varargin{:});
  
  % add path to HMM functions
  addpath(genpath(inferenceOptions.modelPath));

  %% %%%%%%%%%%%%%%%%%%%%%% Load trace data set %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if ~inferenceOptions.ProteinBinFlag
    if isempty(inferenceOptions.dt)
        load([InputDataPath '/spot_struct.mat'],'spot_struct') % load data
    else
        load([InputDataPath '/spot_struct_dt',num2str(inferenceOptions.dt),'.mat'],'spot_struct') % load data
    end
    analysis_traces = spot_struct;
    clear spot_struct
  else
    if ~isempty(inferenceOptions.dt)
        error('Flexible dt is not currently supported for Protein analysis.')
    end
    try
        load([InputDataPath '/spot_struct_protein.mat'],'spot_struct_protein') % load data
        analysis_traces = spot_struct_protein;
    catch
        warning("Couldn't find spot_struct_protein, loading spot_struct instead")
        load([InputDataPath '/spot_struct.mat'],'spot_struct') % load data
        analysis_traces = spot_struct;
    end
    
    clear spot_struct_prote
  end
  
  % check for consistency
  if ~isfield(analysis_traces,'APPosParticle') && inferenceOptions.apBinFlag
    warning('AP-binned option is selected, but data do not contain AP position info')
    return
  end
  %% %%%%%%%%%%%%%%%%%%%%%% Generate output directory %%%%%%%%%%%%%%%%%%%%%%%
  % Set write path (inference results are now written to external directory)
  if inferenceOptions.fluo3DFlag
      fluoSuffix = 'f3D';
  else
      fluoSuffix = 'f2D';
  end
  addSuffix = '';
  if ~isempty(inferenceOptions.AdditionalGroupingVariable)
    addSuffix = ['_' inferenceOptions.AdditionalGroupingVariable];
  end
  
  % generate directory    
  outSuffix =  ['cpHMM_results' filesep 'w' num2str(inferenceOptions.nSteps) '_K' num2str(inferenceOptions.nStates) '_p' ...
    num2str(inferenceOptions.ProteinBinFlag) '_ap' num2str(length(inferenceOptions.apBins)-1) ...
    '_t' num2str(length(inferenceOptions.timeBins)) '_' fluoSuffix addSuffix]; 

  if isfield(inferenceOptions, 'dt')
      if ~isempty(inferenceOptions.dt)
          outSuffix = [outSuffix, '_dt', num2str(inferenceOptions.dt)];
      end
  end
  
  outSuffix = [outSuffix filesep];
  

  % set write path
  outPrefix = [OutputDataPath filesep];
  outDir = [outPrefix outSuffix];
  mkdir(outDir);

  %% %%%%%%%%%%%%%%%%%%%%%% Process/filter trace data %%%%%%%%%%%%%%%%%%%%%%%
  [trace_struct_filtered, indexInfo, inferenceOptions] = filterTraces(inferenceOptions,analysis_traces);

  %% %%%%%%%%%%%%%%%%%%%%%% Conduct cpHMM Inference %%%%%%%%%%%%%%%%%%%%%%%%%
  rng('shuffle'); % ensure we don't repeat bootstrap samples across different replicates
  if ~inferenceOptions.savioFlag
    h = waitbar(0,'Conducting cpHMM inference...');
  end
  % save inference options file
  inferenceOptions.indexInfo = indexInfo;
  save([outDir 'inferenceOptions.mat'],'inferenceOptions');
  
  % iterate through designated groups
  inferenceIDVecShuffled =  randsample(1:length(indexInfo.indexVecUnique),length(indexInfo.indexVecUnique),false);
  for t = inferenceIDVecShuffled
    
      % find subset of eligible traces
      iter_filter = indexInfo.indexList == indexInfo.indexVecUnique(t);
      timeBin = indexInfo.time_group_vec(t);
      
      % iterate through bootstraps
      for b = 1:inferenceOptions.nBoots

          % Initialize structures to store results
          local_struct_temp = struct;  % stores results from individual iterations
          output = struct;  % structure that will be saved to file               

          % Extract subset of traces relevant to this subgroup       
          inference_set = trace_struct_filtered(iter_filter);                
          set_size = length([inference_set.fluo]);  

          skip_flag = 0;
          if set_size < inferenceOptions.minDPperInf                    
              skip_flag = 1;                    
              warning('Too few data points. Skipping')                                    
          else 
              sample_index = 1:length(inference_set);

              %% take bootstrap sample
              ndp = 0;    
              sample_ids = [];                    

              %Reset bootstrap size to be on order of set size for small bins            
              inferenceOptions.SampleSize(t) = min([inferenceOptions.SampleSize ceil(set_size/100)*100]);

              % randomly draw traces
              if length(sample_index) > 1
                  while ndp < inferenceOptions.SampleSize(t)
                      tr_id = randsample(sample_index,1);
                      sample_ids = [sample_ids tr_id];
                      ndp = ndp + length(inference_set(tr_id).time);
                  end
              else
                  sample_ids = sample_index;
              end
              % add them to data cells
              fluo_data = cell([length(sample_ids), 1]);    
              time_data = cell([length(sample_ids), 1]);    
              sample_particles = [inference_set(sample_ids).particleID];
              for tr = 1:length(sample_ids)
                  fluo_data{tr} = inference_set(sample_ids(tr)).fluo;                    
                  time_data{tr} = inference_set(sample_ids(tr)).time;                    
              end            

              %% Random initialization of model parameters
              param_init = initialize_random(inferenceOptions.nStates, inferenceOptions.nSteps, fluo_data);

              % Approximate inference assuming iid data for param initialization                
              local_iid_out = local_em_iid_reduced_memory(fluo_data, param_init.v, ...
                                  param_init.noise, inferenceOptions.nStates, inferenceOptions.nSteps, inferenceOptions.alpha, 500, 1e-4);

              noise_iid = 1/sqrt(exp(local_iid_out.lambda_log));
              v_iid = exp(local_iid_out.v_logs);  

              %% create parallel pool if one does not already exist
              p = gcp('nocreate');
              if isempty(p)
                  parpool(inferenceOptions.maxWorkers); %6 is the number of cores the Garcia lab server can reasonably handle per user.
              elseif p.NumWorkers > inferenceOptions.maxWorkers
                  delete(gcp('nocreate')); % if pool with too many workers, delete and restart
                  parpool(inferenceOptions.maxWorkers);
              end

              %% conduct cpHMM inference
              parfor i_local = 1:inferenceOptions.n_localEM % Parallel Local EM 

                  % Random initialization of model parameters
                  param_init = initialize_random_with_priors(inferenceOptions.nStates, noise_iid, v_iid);

                  % Get Intial Values
                  pi0_log_init = log(param_init.pi0);
                  A_log_init = log(param_init.A);
                  v_init = param_init.v;                        
                  noise_init = param_init.noise;

                  %--------------------LocalEM Call-------------------------%
                  if ~inferenceOptions.truncInference(timeBin)
                    local_out = local_em_MS2_reduced_memory(fluo_data, ...
                        v_init, noise_init, pi0_log_init', A_log_init, inferenceOptions.nStates, inferenceOptions.nSteps, ...
                        inferenceOptions.alpha, inferenceOptions.nStepsMax, inferenceOptions.eps);                    
                  else
                    local_out = local_em_MS2_reduced_memory_truncated(fluo_data, ...
                          v_init, noise_init, pi0_log_init', A_log_init, inferenceOptions.nStates, inferenceOptions.nSteps, ...
                      inferenceOptions.alpha, inferenceOptions.nStepsMax, inferenceOptions.eps);  
                  end
                  %---------------------------------------------------------%                
                  % Save Results                 
                  local_struct_temp(i_local).subset_id = i_local;
                  local_struct_temp(i_local).logL = local_out.logL;                
                  local_struct_temp(i_local).A = exp(local_out.A_log);
                  local_struct_temp(i_local).v = exp(local_out.v_logs).*local_out.v_signs;
                  local_struct_temp(i_local).r = exp(local_out.v_logs).*local_out.v_signs / inferenceOptions.Tres;                                
                  local_struct_temp(i_local).noise = 1/exp(local_out.lambda_log);
                  local_struct_temp(i_local).pi0 = exp(local_out.pi0_log);
                  local_struct_temp(i_local).total_stinferenceOptions.eps = local_out.n_iter;               
                  local_struct_temp(i_local).soft_struct = local_out.soft_struct;               
              end

              %% Record output
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
              output.timeBin = indexInfo.time_group_vec(t);
              output.apBin = indexInfo.ap_group_vec(t);
              output.additionalBin = indexInfo.additional_group_vec(t);
              output.groupID = t;                            
              
              output.additionalBinVar = inferenceOptions.AdditionalGroupingVariable;
              output.intensityBin = indexInfo.intensity_group_vec(t);
              output.intensityBinVar = inferenceOptions.intensityBinVar;
              output.ProteinBinFlag = inferenceOptions.ProteinBinFlag;
              output.FluoBinFlag = inferenceOptions.FluoBinFlag;
              
              output.truncInference = inferenceOptions.truncInference(timeBin);
              output.iter_id = b;                        
              output.particle_ids = sample_particles;            
              output.N = ndp;

              % save inference data used
              output.fluo_data = fluo_data;
              output.time_data = time_data;
          end
          output.skip_flag = skip_flag;

          %% Determine unique filename and sace

          % Generate filenames            
          fName_sub = ['hmm_results_group' sprintf('%03d',t) '_rep'];
%           file_list = dir([outDir fName_sub '*']);
          % Get largest sub-id
%           if isempty(file_list) 
%             repNum = 1;
%           else
%             repNumList = zeros(size(file_list));
%             start_index = length(fName_sub)+1;
%             for f = 1:length(file_list)
%               repNumList(f) = str2double(file_list(f).name(start_index:start_index+2));
%             end
%             repNum = max(repNumList)+1;
%           end
          % generate random string
          rand_string = strrep(num2str(randsample(1:9,5,true)),' ','');
          % save
          out_file = [outDir fName_sub rand_string];          
          save([out_file '.mat'], 'output');           
      end 
      
      if ~inferenceOptions.savioFlag
        waitbar(t/length(inferenceIDVecShuffled),h)
      end
  end
  if ~inferenceOptions.savioFlag
    delete(h);
  end

