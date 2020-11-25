function analysis_traces_interp = addTraceAnalysisFields(analysis_traces,singleTraceFits)

    % restrict to traces with spot info
    analysis_traces = analysis_traces(~isnan([analysis_traces.particleID]));
    
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
            analysis_traces(a).hasHMMInfo = 1;
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

        else
            analysis_traces(a).hasHMMInfo = 0;
            for f = 1:length(fitFields)                    
                analysis_traces_interp(a).(newNames{f}) = NaN(size(refTime));
            end
        end

    end                   