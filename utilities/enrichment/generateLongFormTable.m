function resultsTable = generateLongFormTable(analysis_traces,timeGrid,singleTraceFits,resultsPath,infName)

    disp('generating longform table...')
    % generate longform dataset interpolated to regular time grid
    keepFields = {'xPosNucleus','yPosNucleus','xPosParticle','yPosParticle','ncID','particleID','qcFlag'};
    extendFlags = [0 0 0 0 1 1 1];
    if isfield(analysis_traces,'Stripe')
        keepFields(end+1) = {'Stripe'};
        extendFlags(end+1) = 1;
    end
    
    if isfield(analysis_traces,'ectopicFlag')
        keepFields(end+1) = {'ectopicFlag'};
        extendFlags(end+1) = 1;
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

    if isfield(analysis_traces,'APPosParticleNorm')
        for a = 1:length(analysis_traces)
            analysis_traces_interp(a).APPosNorm = analysis_traces(a).APPosParticleNorm;
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

        else
            for f = 1:length(fitFields)                    
                analysis_traces_interp(a).(newNames{f}) = NaN(size(refTime));
            end
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
%     setVec = floor(resultsTable.ncID);
%     setIndex = unique(setVec(~isnan(setVec)));
%     resultsTable.apPosNormNucleus = NaN(size(setVec));
%     for s = 2
%         setFilter = setVec==setIndex(s);
%         % build linear model
%         apNormTrue = resultsTable.APPosParticleNorm(setFilter)';
%         inputs = [resultsTable.xPosParticle(setFilter) resultsTable.yPosParticle(setFilter) ...
%                   resultsTable.time(setFilter)];
%         apModel = fitlm(inputs, apNormTrue);
%         % use model to make normalized nucleus position variable
%         nucleusInputs = [resultsTable.xPosNucleus(setFilter) resultsTable.yPosNucleus(setFilter) ...
%                   resultsTable.time(setFilter)];
%         resultsTable.apPosNormNucleus(setFilter) = predict(apModel,nucleusInputs);                    
%     end
    % save
    disp('saving...')
    writetable(resultsTable,[resultsPath 'singleTraceFits_' infName '_longform.csv'])
    disp('done.')