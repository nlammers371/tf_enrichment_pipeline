function analysis_traces_interp = addTraceAnalysisFields(analysis_traces,singleTraceFitsViterbi,singleTraceFitsSS,inferenceOptions)


    % generate alpha kernel for estimating predicted HMM fluroescence
    nSteps = inferenceOptions.nSteps;
    alpha =inferenceOptions.alpha;
    Tres = inferenceOptions.Tres;
    
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
    
    % restrict to traces with spot info
    analysis_traces = analysis_traces(~isnan([analysis_traces.particleID]));
    
    % generate longform dataset interpolated to regular time grid
    keepFields = {'xPosNucleus','yPosNucleus','xPosParticle','yPosParticle','ncID','particleID','TraceQCFlag'};
    extendFlags = [0 0 0 0 1 1 1];
    if isfield(analysis_traces,'APPosParticleNorm')
        keepFields(end+1) = {'APPosParticleNorm'};
        extendFlags(end+1) = 0;
    end

    analysis_traces_interp = struct;
    
    % First, generate regularized time field
    for a = 1:length(analysis_traces)
        analysis_traces_interp(a).time = analysis_traces(a).timeInterp;
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
    if isfield(analysis_traces,'APPosParticleInterp')
        prevInterpFields = {'fluoInterp','APPosParticleInterp'};
        newInterpNames = {'fluo','APPosParticle'};
    else
        prevInterpFields = {'fluoInterp'};
        newInterpNames = {'fluo'};
    end
    
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
            analysis_traces_interp(a).(newInterpNames{k}) = vecNew;
        end
    end

    % lastly, add new fields from viterbi fits
    fit_particles = [singleTraceFitsViterbi.particleID];
    nStates = length(singleTraceFitsSS(1).r);
    vFitFields = {'z_viterbi','v_viterbi','fluo_viterbi'};
    newNames = {'promoter_state','initiation_rate','predicted_fluo'};
    for a = 1:length(analysis_traces_interp)
        % add viterbi and soft fit info
        fitIndexViterbi = find(fit_particles==analysis_traces(a).particleID);

        refTime = analysis_traces_interp(a).time;

        if ~isempty(fitIndexViterbi)
            analysis_traces(a).hasHMMInfo = 1;
            fitTime = singleTraceFitsViterbi(fitIndexViterbi).time;
            toFilter = ismember(refTime,fitTime);
            toFilterTr = ismember(refTime(1:end-1),fitTime);
            if ~all(ismember(fitTime,refTime))
              error('issue with time cross-referencing')
            end
            % viterbi fields
            for f = 1:length(vFitFields)
                vecNew = NaN(size(refTime));
                vecFit = singleTraceFitsViterbi(fitIndexViterbi).(vFitFields{f});
                vecNew(toFilter) = vecFit;
                analysis_traces_interp(a).(newNames{f}) = vecNew;
            end
            
            % soft fit fields
            analysis_traces_interp(a).state_prob_array = NaN(length(refTime),nStates);
            analysis_traces_interp(a).state_prob_array(toFilter,:) = singleTraceFitsSS(fitIndexViterbi).state_prob_array';
            
            analysis_traces_interp(a).transition_prob_array = NaN(nStates,nStates,length(refTime)-1);
            analysis_traces_interp(a).transition_prob_array(:,:,toFilterTr) = singleTraceFitsSS(fitIndexViterbi).transition_prob_array;
            
            analysis_traces_interp(a).r_vec_soft = analysis_traces_interp(a).state_prob_array*singleTraceFitsSS(fitIndexViterbi).r;
            
            % get predicted fluo
            fluo_pd_MS2 = conv(alpha_kernel, analysis_traces_interp(a).r_vec_soft,'Full');
            analysis_traces_interp(a).fluo_pd_MS2 = fluo_pd_MS2(1:end-nSteps+1);
            
            full_kernel = Tres*ones(size(alpha_kernel));
            fluo_pd_full = conv(full_kernel, analysis_traces_interp(a).r_vec_soft,'Full');
            analysis_traces_interp(a).fluo_pd_full = fluo_pd_full(1:end-nSteps+1);
            
            
        else
            analysis_traces(a).hasHMMInfo = 0;
            for f = 1:length(vFitFields)                    
                analysis_traces_interp(a).(newNames{f}) = NaN(size(refTime));
            end
        end

    end                   