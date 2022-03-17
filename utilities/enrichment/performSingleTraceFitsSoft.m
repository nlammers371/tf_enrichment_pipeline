function singleTraceFitsSS = performSingleTraceFitsSoft(compiledResults, inferenceOptions, bootstrap_flag, ...
                                      analysis_traces, trace_particle_index, nWorkersMax, resultsPath, infName)

                          
    % get basic inference characteristics
    nStates = inferenceOptions.nStates;
    nSteps = inferenceOptions.nSteps;
    alpha = inferenceOptions.alpha;
    Tres = inferenceOptions.Tres;
    eps = 1e-4;
    % get grouping vectors
    groupID_vec = compiledResults.groupID_index;
    groupParticles = compiledResults.particle_ids;

    % initialize trace structure
    singleTraceFitsSS = struct; 

%     % initialize random number generator
%     rng(random_seed);

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

            soft_fits = struct;
            r_vec = compiledResults.r_array_mean(:,g);
            parfor f = 1:length(fluo_values)
    %             waitbar(f/numel(fluo_values),h);
%                 viterbi_out = viterbi (fluo_values{f}, v', sigma, pi0_log,A_log, nStates, nSteps, alpha);
                local_em_outputs = local_em_MS2_reduced_memory (fluo_values(f), ...
                                v', sigma, pi0_log, A_log, nStates, nSteps, alpha, 1, eps);
                
                soft_fits(f).state_prob_array = exp(local_em_outputs.soft_struct.p_z_log_soft{1});
                soft_fits(f).transition_prob_array = exp(local_em_outputs.soft_struct.p_zz_log_soft{1});
                soft_fits(f).r = r_vec;
            end

            % map back to trace fit structure               
            for f = 1:length(fluo_values)                    
                fnames = fieldnames(soft_fits(f));
                for fn = 1:numel(fnames)
                    singleTraceFitsSS(nextIndex).(fnames{fn}) = soft_fits(f).(fnames{fn});
                end 
                % add additional info
                singleTraceFitsSS(nextIndex).time = analysis_traces(fit_trace_indices(f)).timeInterp;
                singleTraceFitsSS(nextIndex).fluo_exp = fluo_values{f};
                singleTraceFitsSS(nextIndex).particleID = fit_trace_particles(f);
                singleTraceFitsSS(nextIndex).groupID = groupID_vec(g);
                singleTraceFitsSS(nextIndex).bootID = j;
                singleTraceFitsSS(nextIndex).infID = infID;
                if isfield(analysis_traces,'APPosParticleInterp')
                    singleTraceFitsSS(nextIndex).APPos = analysis_traces(fit_trace_indices(f)).APPosParticleInterp;
                end

                nextIndex = length(singleTraceFitsSS) + 1;
            end                

        end   
        waitbar(g/length(groupID_vec),wb);
    end
    delete(wb);
    
    save([resultsPath 'singleTraceFitsSS_' infName '.mat'],'singleTraceFitsSS')