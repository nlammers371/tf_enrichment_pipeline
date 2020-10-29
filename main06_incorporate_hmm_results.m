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

function main06_incorporate_hmm_results(projectName,varargin)

close all
addpath(genpath('utilities'))

makeLongFormSet = 0;
random_seed = 432;
myCluster = parcluster('local');
nWorkersMax = ceil(myCluster.NumWorkers/2);
bootstrap_flag = 0;
n_bootstraps = 1;

% check for optional inputs
for i = 1:(numel(varargin)-1)  
    if i ~= numel(varargin)        
        eval([varargin{i} '=varargin{i+1};']);                
    end    
end

% get path to results
liveProject = LiveProject(projectName);
resultsDir = [liveProject.dataPath filesep 'cpHMM_results' filesep];

% get list of all inference subdirectories. By default, we'll generate
% summaries for all non-empty inference sub-directory
infDirList = dir([resultsDir 'w*']);% get list of all inference subdirectories. By default, we'll generate
   
for inf = 1:length(infDirList)
  
    % get path to results
    resultsPath = [infDirList(inf).folder filesep];
    
    % load inference options
    load([resultsPath infDirList(inf).name filesep 'inferenceOptions.mat'])
    if length(fieldnames(inferenceOptions))==1
        inferenceOptions = inferenceOptions.inferenceOptions;
    end
    
    % extract key inference parameters     
    nStates = inferenceOptions.nStates;
    nSteps = inferenceOptions.nSteps;
    alpha = inferenceOptions.alpha;
    Tres = inferenceOptions.Tres;
    timeGrid = 0:Tres:60*60;
    maxDT = 1.2*Tres; % maximum distance from observed data point
    
    % load compiled inference files
    load([resultsPath 'compiledResults_' infDirList(inf).name '.mat']);
    
    % load corresponding trace structure
    if ~inferenceOptions.ProteinBinFlag
        load([liveProject.dataPath filesep 'nucleus_struct.mat'])
        analysis_traces = nucleus_struct;
        clear nucleus_struct;
    else
        load([liveProject.dataPath filesep 'spot_struct_protein.mat'])
        analysis_traces = spot_struct_protein;
        clear spot_struct_protein;
    end
           
    % generate alpha kernel for estimating predicted HMM fluroescence
    alpha_kernel = NaN(1,nSteps);
    for i = 1:nSteps
        if i < alpha
            alpha_kernel(i) = ((i-1) / alpha  + .5 * 1/alpha )*Tres;
        elseif i > alpha && (i-1) < alpha
            alpha_kernel(i) = Tres*(1 - .5*(alpha-i+1)*(1-(i-1)/alpha));
        else
            alpha_kernel(i) = Tres;
        end
    end

    % check for existence of fit structure
    trace_fit_flag = 1;
    inf_files = dir([resultsPath infDirList(inf).name filesep 'hmm_results*.mat']);
    if exist([resultsPath 'singleTraceFits_' infDirList(inf).name '.mat']) > 0
        fit_props = dir([resultsPath 'singleTraceFits_' infDirList(inf).name '.mat']);
        fit_date = datenum(fit_props(1).date);
        hmm_date = datenum(inf_files(1).date);
        if fit_date > hmm_date
            trace_fit_flag = 0;
        end
    end
    
    % get list of particle IDs
    trace_particle_index = [analysis_traces.particleID];
    
    % perform viterbi trace decoding if necessary
    if trace_fit_flag    
        
        groupID_vec = compiledResults.groupID_index;
        groupParticles = compiledResults.particle_ids;
%         nParticles = length([groupParticles{:}]);
        
        % initialize trace structure
        singleTraceFits = struct; 
        
        % initialize random number generator
        rng(random_seed);
        
        % loop through inference groups
        wb = waitbar(0, 'Conducting single trace fits...');
        nextIndex = 1;
        for g = 1:length(groupID_vec)
            infIDs_to_use = NaN;
            
            if bootstrap_flag
                infIDs_to_use = randsample(find(~compiledResults.outlier_flags{g}),n_boostraps,true);
            end
            
            for j = 1:length(infIDs_to_use)
                infID = infIDs_to_use(j);
                if bootstrap_flag                    
                    A_log = log(compiledResults.A_results{g}(:,:,infID));
                    v = compiledResults.r_results{g}(:,infID)*Tres;
                    sigma = compiledResults.sigma_results{g}(infID);
                    pi0_log = log(compiledResults.pi0_results{g}(:,infID)); 
                else
                    A_log = log(compiledResults.A_array_mean(:,:,g));
                    v = compiledResults.r_array_mean(:,g)*Tres;
                    sigma = compiledResults.sigma_vec_mean(g);
                    pi0_log = log(compiledResults.pi0_array_mean(:,g));
                end

                % obtain subset of valid traces                
                fit_trace_indices = find(ismember(trace_particle_index,groupParticles{g}));                             
                fit_trace_particles = trace_particle_index(ismember(trace_particle_index,groupParticles{g}));
                
                fluo_values = cell(length(fit_trace_indices),1);                
                for i = 1:numel(fit_trace_indices)
                    if ~inferenceOptions.fluo3DFlag
                        fluo = analysis_traces(fit_trace_indices(i)).fluoInterp;
                    else
                        fluo = analysis_traces(fit_trace_indices(i)).fluo3DInterp;
                    end
%                     start_i = find(~isnan(fluo),1);
%                     stop_i = find(~isnan(fluo),1,'last');
%                     fluo = fluo(start_i:stop_i);
                    fluo_values{i} = fluo;
                end    

                % initialize pool                   
                p = gcp('nocreate'); % If no pool, do not create new one.
                if isempty(p)
                    parpool(nWorkersMax);
                elseif p.NumWorkers~= nWorkersMax
                    delete(p);
                    parpool(nWorkersMax);
                end       
%                 disp('conducting viterbi trace fits...')
                
                v_fits = struct;
                parfor f = 1:length(fluo_values)
        %             waitbar(f/numel(fluo_values),h);
                    viterbi_out = viterbi (fluo_values{f}, v', sigma, pi0_log,A_log, nStates, nSteps, alpha);
                    fnames = fieldnames(viterbi_out);
                    for fn = 1:numel(fnames)
                        v_fits(f).(fnames{fn}) = viterbi_out.(fnames{fn});
                    end        
                end
                
                % map back to trace fit structure               
                for f = 1:length(fluo_values)                    
                    fnames = fieldnames(v_fits(f));
                    for fn = 1:numel(fnames)
                        singleTraceFits(nextIndex).(fnames{fn}) = v_fits(f).(fnames{fn});
                    end 
                    % add additional info
                    singleTraceFits(nextIndex).time = analysis_traces(fit_trace_indices(f)).timeInterp;
                    singleTraceFits(nextIndex).fluo_exp = fluo_values{f};
                    singleTraceFits(nextIndex).particleID = fit_trace_particles(f);
                    singleTraceFits(nextIndex).groupID = groupID_vec(g);
                    singleTraceFits(nextIndex).bootID = j;
                    singleTraceFits(nextIndex).infID = infID;
                    if isfield(analysis_traces,'APPosParticleInterp')
                        singleTraceFits(nextIndex).APPos = analysis_traces(fit_trace_indices(f)).APPosParticleInterp;
                    end
%                     % NL: this is specific to Augusto's Eve data
%                     if isfield(analysis_traces,'APPosParticleNorm')
%                         singleTraceFits(nextIndex).APPosNorm = analysis_traces(fit_trace_indices(f)).APPosParticleInterp;
%                     end
                    
                    nextIndex = length(singleTraceFits) + 1;
                end                

            end   
            waitbar(g/length(groupID_vec),wb);
        end
        delete(wb);
        save([resultsPath 'singleTraceFits_' infDirList(inf).name '.mat'],'singleTraceFits')
    else
        load([resultsPath 'singleTraceFits_' infDirList(inf).name '.mat'],'singleTraceFits')
    end
  
    % generate longform dataset
    if true%makeLongFormSet        
        disp('generating longform table...')
        % generate longform dataset interpolated to regular time grid
        keepFields = {'xPosNucleus','yPosNucleus','xPosParticle','yPosParticle','ncID','particleID','qcFlag'};
        extendFlags = [0 0 0 0 1 1 1];
        if isfield(analysis_traces,'APPosParticleNorm')
            keepFields(end+1) = {'APPosParticleNorm'};
            extendFlags(end+1) = 0;
        end
        
        
        analysis_traces_interp = struct;
        % First, generate regularized time field
        for a = 1:length(analysis_traces)
            analysis_traces_interp(a).time = timeGrid(timeGrid>=analysis_traces(a).time(1)&timeGrid<=analysis_traces(a).time(end));
        end
        
        % now, interpolate/extend selected fields
        for a = 1:length(analysis_traces)
            timeRaw = analysis_traces(a).time;
            timeClean = analysis_traces_interp(a).time;
            for k = 1:length(keepFields)
                vec = analysis_traces(a).(keepFields{k});
                if ~extendFlags(k)
                    nan_ft = ~isnan(vec);
                    if sum(nan_ft) > 1
                        analysis_traces_interp(a).(keepFields{k}) = interp1(timeRaw(nan_ft),double(vec(nan_ft)),timeClean);
                    elseif sum(nan_ft) == 1
                        vecNew = NaN(size(timeClean));
                        [~, mi] = min(abs(timeClean-timeRaw(nan_ft)));
                        vecNew(mi) = vec(nan_ft);
                        analysis_traces_interp(a).(keepFields{k}) = vecNew;
                    else                      
                        analysis_traces_interp(a).(keepFields{k}) = NaN(size(timeClean));
                    end
                else
                    analysis_traces_interp(a).(keepFields{k}) = repelem(vec,length(timeClean));
                end
            end
        end
        % add previously interpolated fields
        prevInterpFields = {'fluoInterp','APPosParticleInterp'};
        for a = 1:length(analysis_traces)
            time1 = analysis_traces(a).timeInterp;
            time2 = analysis_traces_interp(a).time;            
            toFilter = ismember(time2,time1);
            fromFilter = ismember(time1,time2);
            if ~all(fromFilter) && ~all(isnan(time1))
              error('inconsistent time frames')
            end
            for k = 1:length(prevInterpFields)
                vecNew = NaN(size(time2));
                vecOrig = analysis_traces(a).(prevInterpFields{k});
                vecNew(toFilter) = vecOrig(fromFilter);
                analysis_traces_interp(a).(prevInterpFields{k}) = vecNew;
            end
        end
                
        % lastly, add new fields from viterbi fits
        fit_particles = [singleTraceFits.particleID];
        fitFields = {'z_viterbi','fluo_viterbi'};
        newNames = {'promoter_state','predicted_fluo'};
        for a = 1:length(analysis_traces_interp)
            fitIndex = find(fit_particles==analysis_traces(a).particleID);
                                   
            refTime = analysis_traces_interp(a).time;
            
            if ~isempty(fitIndex)
                fitTime = singleTraceFits(fitIndex).time;
                toFilter = ismember(refTime,fitTime);
                if ~all(ismember(fitTime,refTime))
                  error('issue with time cross-referencing')
                end
                for f = 1:length(fitFields)
                    vecNew = NaN(size(refTime));
                    vecFit = singleTraceFits(fitIndex).(fitFields{f});
                    vecNew(toFilter) = vecFit;
                    analysis_traces_interp(a).(newNames{f}) = vecNew;
                end
                
                % add additional grouping info if appropriate
%                 if ~isempty(inferenceOptions.AdditionalGroupingVariable)
%                    analysis_traces_interp(a).(inferenceOptions.AdditionalGroupingVariable) = repelem(...
%                       singleTraceFits(fitIndex).groupID,length(refTime));
%                 end
            else
                for f = 1:length(fitFields)                    
                    analysis_traces_interp(a).(newNames{f}) = NaN(size(refTime));
                end
                
%                 if ~isempty(inferenceOptions.AdditionalGroupingVariable)
%                    analysis_traces_interp(a).(inferenceOptions.AdditionalGroupingVariable) = NaN(...
%                       size(refTime));
%                 end
            end
            
            
        end               
                
        % make longform table
        resultsTable = struct;
        fnames = fieldnames(analysis_traces_interp);
        for f = 1:length(fnames)
            resultsTable.(fnames{f}) = [analysis_traces_interp.(fnames{f})]';
        end
        resultsTable.particleID(isnan(resultsTable.fluoInterp)) = NaN;
        resultsTable = struct2table(resultsTable);
        
        % try to infer normalized nucleus AP by building model of norm AP
        % assignment
        setVec = floor(resultsTable.ncID);
        setIndex = unique(setVec(~isnan(setVec)));
        resultsTable.apPosNormNucleus = NaN(size(setVec));
        for s = 2
            setFilter = setVec==setIndex(s);
            % build linear model
            apNormTrue = resultsTable.APPosParticleNorm(setFilter)';
            inputs = [resultsTable.xPosParticle(setFilter) resultsTable.yPosParticle(setFilter) ...
                      resultsTable.time(setFilter)];
            apModel = fitlm(inputs, apNormTrue);
            % use model to make normalized nucleus position variable
            nucleusInputs = [resultsTable.xPosNucleus(setFilter) resultsTable.yPosNucleus(setFilter) ...
                      resultsTable.time(setFilter)];
            resultsTable.apPosNormNucleus(setFilter) = predict(apModel,nucleusInputs);                    
        end
        % save
        disp('saving...')
        writetable(resultsTable,[resultsPath 'singleTraceFits_' infDirList(inf).name '_longform.csv'])
        disp('done.')
    end

end
%     trace_particle_index_orig = [analysis_traces.ParticleID];
%     %%% now extract corresponding hmm traces
%     disp('building input/output dataset...')
%     hmm_input_output = [];
%     for inf = 1:numel(singleTraceFits)
% 
%         soft_fits = singleTraceFits(inf).soft_fits;
%         viterbi_fits = singleTraceFits(inf).viterbi_fits;
%         inference_id_vec = singleTraceFits(inf).inference_id_vec;
%         inference_particles = singleTraceFits(inf).particle_index;
%         fit_trace_indices = find(ismember(trace_particle_index,inference_particles));
%         for i = fit_trace_indices
%             % initialize temporary structure to store results
%             ParticleID = trace_particle_index(i);   
%             traceIndex = find(inference_particles==ParticleID);
%             if isempty(traceIndex)
%               error('uh oh')
%             end
%             temp = struct;
%             % extract relevant vectors from protein struct    
%             % these quantities have not been interpolated
%             if fluo_dim == 3   
%                 ff_pt = analysis_traces(i).fluo3D;
%                 master_fluo = analysis_traces(i).fluo3D_interp;            
%             elseif fluo_dim == 2
%                 ff_pt = analysis_traces(i).fluo;
%                 master_fluo = analysis_traces(i).fluo_interp;            
%             end        
%             if protein_dim == 3
%                 sp_pt = analysis_traces(i).spot_protein_vec_3d;
%                 sr_pt = analysis_traces(i).serial_null_protein_vec_3d; 
%             elseif protein_dim == 2
%                 sp_pt = analysis_traces(i).spot_protein_vec;
%                 sr_pt = analysis_traces(i).serial_null_protein_vec; 
%             end
%             mcp_pt = analysis_traces(i).spot_mcp_vec;         
%             nn_pt = analysis_traces(i).edge_null_protein_vec;
%             mf_pt_mf = analysis_traces(i).mf_null_protein_vec;          
%             tt_pt = analysis_traces(i).time;
% 
%             x_pt = analysis_traces(i).xPosParticle;  
%             y_pt = analysis_traces(i).yPosParticle;  
%             ap_pt = analysis_traces(i).APPosParticle;  
% 
% 
%             if sum(~isnan(mf_pt_mf)) > minDP && sum(~isnan(sr_pt)) > minDP && sum(~isnan(sp_pt)) > minDP              
%                 % extract interpolated fluorescence and time vectors
%                 master_time = analysis_traces(i).time_interp;
% 
%                 % check for mismatch between analysis_traces and
%                 % analysis_traces...this is due to a dumb mistake on my
%                 % part
%                 time_vec_orig = analysis_traces(trace_particle_index_orig==ParticleID).time_interp;
%                 if ~isequal(master_time,time_vec_orig)
%                   master_time = time_vec_orig;
%                   analysis_traces(i).time_interp = master_time;
% 
%                   master_fluo = analysis_traces(trace_particle_index_orig==ParticleID).fluo_interp;
%                   analysis_traces(i).fluo_interp = master_fluo;
%                 end
%                 % extract position vectors (used for selecting nearest neighbor)
%                 x_nc = double(analysis_traces(i).xPos);
%                 y_nc = double(analysis_traces(i).yPos);
%                 temp.xPosMean = nanmean(x_nc(~isnan(ff_pt)));
%                 temp.yPosMean = nanmean(y_nc(~isnan(ff_pt)));
% 
%                 % record time, space, and fluo vars
%                 start_i = find(~isnan(master_fluo),1);
%                 stop_i = find(~isnan(master_fluo),1,'last');
%                 temp.time = master_time(start_i:stop_i);
%                 temp.fluo = master_fluo(start_i:stop_i);   
%                 temp.fluo_raw = ff_pt;      
% 
%                 % extract useful HMM inference parameters             
%                 inference_id = inference_id_vec(traceIndex); % inference id
%                 [r,ri] = sort(inference_results(inference_id).r); % enforce rank ordering of states
%                 z = exp(soft_fits{traceIndex});    
%                 temp.z_mat = z(ri,:)';    
%                 temp.r_mat = z(ri,:)'.*r';
%                 temp.r_inf = r';
%                 temp.r_vec = sum(temp.r_mat,2)';
%                 [~,z_vec] = max(temp.z_mat,[],2);
%                 temp.z_vec = z_vec; 
%                 if length(z_vec)~=length(temp.fluo)
%                   error('goddammit')
%                 end
%                 % extract viterbi fits
%                 temp.z_viterbi = viterbi_fits(traceIndex).z_viterbi;
%                 temp.f_viterbi = viterbi_fits(traceIndex).fluo_viterbi;
% 
%                 % make predicted fluo vec (for consistency checks)
%                 fluo_hmm = conv(temp.r_vec,alpha_kernel);
%                 temp.fluo_hmm = fluo_hmm(1:numel(temp.r_vec));        
% 
%                 % checks using mcp channel to ensure that we are correctly matching
%                 % particles and time frames
%                 temp.mcp_check = interp1(tt_pt(~isnan(mcp_pt)),mcp_pt(~isnan(mcp_pt)),temp.time);
%                 temp.fluo_check = interp1(tt_pt(~isnan(ff_pt)),ff_pt(~isnan(ff_pt)),temp.time);
% 
%                 % record raw data vectors
%                 temp.spot_protein_raw = sp_pt;        
%                 temp.mf_protein_raw = mf_pt_mf;
%                 temp.null_protein_raw = nn_pt;
%                 temp.serial_protein_raw = sr_pt;  
%                 temp.time_raw = tt_pt;  
%                 temp.xPos_raw = x_nc;
%                 temp.yPos_raw = y_nc;
% 
%                 % interpolate protein information such that it coincides with HMM
%                 % inference results    
%                 temp.spot_protein = interp1(tt_pt(~isnan(sp_pt)),sp_pt(~isnan(sp_pt)),temp.time);                    
%                 temp.serial_protein = interp1(tt_pt(~isnan(sr_pt)),sr_pt(~isnan(sr_pt)),temp.time);            
%                 temp.null_protein = interp1(tt_pt(~isnan(nn_pt)),nn_pt(~isnan(nn_pt)),temp.time);
%                 temp.mf_protein = interp1(tt_pt(~isnan(mf_pt_mf)),mf_pt_mf(~isnan(mf_pt_mf)),temp.time);
% 
%                 % interpolate position info
%                 temp.xPos = interp1(tt_pt(~isnan(x_nc)),x_nc(~isnan(x_nc)),temp.time);        
%                 temp.yPos = interp1(tt_pt(~isnan(y_nc)),y_nc(~isnan(y_nc)),temp.time);
% 
%                 temp.xPosParticle = interp1(tt_pt(~isnan(x_pt)),x_pt(~isnan(x_pt)),temp.time);        
%                 temp.yPosParticle = interp1(tt_pt(~isnan(y_pt)),y_pt(~isnan(y_pt)),temp.time);
%                 temp.apPosParticle = interp1(tt_pt(~isnan(ap_pt)),ap_pt(~isnan(ap_pt)),temp.time);
%                 % generate flag var indicating interpolated obs that are too far from 
%                 % true points
%                 input_times = tt_pt(~isnan(sp_pt));
%                 dt_vec_gap = NaN(size(temp.time));
%                 for t = 1:numel(dt_vec_gap)
%                     dt_vec_gap(t) = min(abs(temp.time(t)-input_times));   
%                 end
%                 temp.dt_filter_gap = dt_vec_gap > maxDT;            
%                 % record general info for later use
%                 temp.ParticleID = ParticleID; 
%                 temp.Tres = Tres;
%                 temp.maxDT = maxDT;
%                 temp.InferenceID = inf;    
%                 hmm_input_output  = [hmm_input_output temp];
%             end
%             % increment
%     %         iter = iter + 1;
%         end
%     end
% 
%     % find nearest neighbor particles
%     % use name nearest neighbor for each bootstrap instance
%     n_unique = numel(hmm_input_output) / n_boots;%numel(inference_results);
%     start_time_vec = NaN(1,n_unique);
%     stop_time_vec = NaN(1,n_unique);
%     set_vec = NaN(1,n_unique);
%     for i = 1:n_unique
%         dt_flag = hmm_input_output(i).dt_filter_gap;
%         t_vec = hmm_input_output(i).time(~dt_flag);
%         start_time_vec(i) = min(t_vec);
%         stop_time_vec(i) = max(t_vec);
%         set_vec(i) = floor(hmm_input_output(i).ParticleID);
%     end
% 
%     % xy nearest neighbor calculations
%     dist_mat_x = pdist2([hmm_input_output(1:n_unique).xPosMean]',[hmm_input_output(1:n_unique).xPosMean]');
%     dist_mat_y = pdist2([hmm_input_output(1:n_unique).yPosMean]',[hmm_input_output(1:n_unique).yPosMean]');
%     dist_mat_r = sqrt(dist_mat_x.^2 + dist_mat_y.^2);
% 
%     % now find closest match for each nucleus
%     for i = 1:n_unique
%         % require that mat trace is (1) from same set, (2) starts and ends
%         % within 3 min of locus trace
%         setID = set_vec(i);  
%         option_filter = ((start_time_vec-start_time_vec(i)) <= 3*60) & ...
%             ((stop_time_vec-stop_time_vec(i)) >= -3*60) & set_vec==setID;        
% 
%         %%% Spatial Nearest Neighbor   
%         time_i = hmm_input_output(i).time;
%         dist_vec = dist_mat_r(i,:);               
%         dist_vec(~option_filter) = NaN;
%         dist_vec(i) = NaN; % remove self
%         [best_r, best_ind_dist] = nanmin(dist_vec);
%         % record vales 
%         time_swap_dist = hmm_input_output(best_ind_dist).time;       
%         % fill
%         swap_ft = ismember(time_swap_dist,time_i);
%         target_ft = ismember(time_i,time_swap_dist);
%         s_pt_dist = NaN(size(time_i));
%         s_pt_dist(target_ft) = hmm_input_output(best_ind_dist).spot_protein(swap_ft);    
%         mf_pt_dist = NaN(size(time_i));
%         mf_pt_dist(target_ft) = hmm_input_output(best_ind_dist).mf_protein(swap_ft);    
%         r_fluo_dist = NaN(size(time_i));
%         r_fluo_dist(target_ft) = hmm_input_output(best_ind_dist).r_vec(swap_ft);
%         dt_filter_dist = true(size(time_i));
%         dt_filter_dist(target_ft) = hmm_input_output(best_ind_dist).dt_filter_gap(swap_ft);
% 
%         % assign to ALL copies
%         for ind = i:n_unique:length(hmm_input_output)
%     %         ind = (inf-1)*n_unique + i;
%             % record dist
%             hmm_input_output(ind).nn_best_r = best_r;
%             hmm_input_output(ind).dist_swap_ind = best_ind_dist;
%             hmm_input_output(ind).dist_swap_spot_protein = s_pt_dist;
%             hmm_input_output(ind).dist_swap_mf_protein = mf_pt_dist;
%             hmm_input_output(ind).dist_swap_hmm = r_fluo_dist;
%             hmm_input_output(ind).dist_swap_dt_filter_gap = dt_filter_dist;
%         end
%     end
% 
%     % save results
%     save([DataPath 'hmm_input_output_w' num2str(w) '_K' num2str(nSteps) '_f' num2str(fluo_dim)  'D_p' num2str(protein_dim) 'D.mat'],'hmm_input_output')
% end