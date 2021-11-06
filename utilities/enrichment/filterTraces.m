function [trace_struct_filtered, indexInfo, inferenceOptions] = filterTraces(inferenceOptions,analysis_traces)
  % this function filters traces and generates grouping variables as
  % appropriate
  % Supported grouping variables include: time, APPositionParticle, and nuclear_protein_vec
    
  % apply QC filter calculated in main01
  analysis_traces = analysis_traces(~isnan([analysis_traces.TraceQCFlag]));
  
  if inferenceOptions.ProteinBinFlag && inferenceOptions.FluoBinFlag  
    error('Protein and spot fluorescence binning are mutually exclusive inference options')
  end
  
  % generate stripped down data structure
  trace_struct_filtered = struct;
  
  for i = 1:length(analysis_traces)
      time = analysis_traces(i).timeInterp;
      if ~isfield(inferenceOptions,'Tres') && length(time)>=2
          inferenceOptions.Tres = time(2)-time(1);
          break
      end
  end
    
  i_pass = 1;
  for i = 1:length(analysis_traces)    
      if inferenceOptions.fluo3DFlag
          fluo = analysis_traces(i).fluo3DInterp;
      else
          fluo = analysis_traces(i).fluoInterp;
      end
      time = analysis_traces(i).timeInterp;          
      fluo_raw = analysis_traces(i).fluo;
      time_raw = analysis_traces(i).time;
      
      if inferenceOptions.apBinFlag
        ap_interp = analysis_traces(i).APPosParticleInterp;
      else
        ap_interp = ones(size(time));
      end
      
      qcFlag = analysis_traces(i).TraceQCFlag || ~inferenceOptions.useQCFlag;      

      for t = 1:length(inferenceOptions.timeBins)     
          % calculate first time at which trace can start
          start_time = time_raw(1);         
          if inferenceOptions.truncInference(t)
            start_time = start_time + (1+inferenceOptions.nSteps)*inferenceOptions.Tres;
          end
          
          timeBounds = inferenceOptions.timeBins{t};

          time_filter_raw = time_raw >= timeBounds(1) & time_raw < timeBounds(2) & time_raw >= start_time;

          nRaw_flag = (sum(~isnan(fluo_raw(time_filter_raw))) >= inferenceOptions.minDP) || inferenceOptions.ignoreNDP;

          if nRaw_flag && qcFlag            

            time_filter_interp = time >= timeBounds(1) & time < timeBounds(2) & time >= start_time;

            % generate new entry    
            trace_struct_filtered(i_pass).fluo = fluo(time_filter_interp);
            trace_struct_filtered(i_pass).time = time(time_filter_interp);
            if inferenceOptions.FluoBinFlag          
                trace_struct_filtered(i_pass).mean_intensity = nanmean(trace_struct_filtered(i_pass).fluo);
            elseif inferenceOptions.ProteinBinFlag
                if isfield(analysis_traces,'nuclear_protein_vecInterp')
                    trace_struct_filtered(i_pass).mean_intensity = nanmean(analysis_traces(i).nuclear_protein_vecInterp(time_filter_interp));                            
                elseif isfield(analysis_traces,'dark_yap')
                    trace_struct_filtered(i_pass).mean_intensity = analysis_traces(i).dark_yap;                            
                end
            end

            trace_struct_filtered(i_pass).particleID = analysis_traces(i).particleID;  
            trace_struct_filtered(i_pass).N = sum(time_filter_interp);                
            trace_struct_filtered(i_pass).timeBin = t;
            trace_struct_filtered(i_pass).interp_filter = time_filter_interp;
            trace_struct_filtered(i_pass).raw_filter = time_filter_raw;
            
            if ~isempty(inferenceOptions.AdditionalGroupingVariable)
              addVarName = inferenceOptions.AdditionalGroupingVariable;
              trace_struct_filtered(i_pass).(addVarName) = analysis_traces(i).(addVarName);
            end
            
            % determin ap bin
            trace_struct_filtered(i_pass).apBin = discretize(mean(ap_interp(time_filter_interp)),inferenceOptions.apBins);
            
            % increment
            i_pass = i_pass + 1;
          end
      end
  end
    
  % make new indexing vectors
  ap_group_vec = [trace_struct_filtered.apBin];
  time_group_vec = [trace_struct_filtered.timeBin];
  
  % initialize dummy grouping variable  
  additional_group_vec = ones(size(ap_group_vec));
  
  % check for additional binning variable (only 1 addtional variable
  % supported for now)  
  if ~isempty(inferenceOptions.AdditionalGroupingVariable)
    % find corresponding grouper variable
    addVarName = inferenceOptions.AdditionalGroupingVariable;
    additional_group_vec = [trace_struct_filtered.(addVarName)];
    inferenceOptions.additionalGroupIDs = unique(additional_group_vec(~isnan(additional_group_vec)));
  end
  
  nan_filter1 = isnan(ap_group_vec) | isnan(time_group_vec) | isnan(additional_group_vec);
  
  % perform protein binning if appropriate
  if inferenceOptions.ProteinBinFlag || inferenceOptions.FluoBinFlag  
    
    for a = 1:length(inferenceOptions.apBins)-1 % Note: if either option is not activated, there will be only one bin                
        for t = 1:length(inferenceOptions.timeBins)         
            for g = 1:length(inferenceOptions.additionalGroupIDs)
                group_ids = find(ap_group_vec==a & time_group_vec==t & additional_group_vec == inferenceOptions.additionalGroupIDs(g));          
                % estimate number of bins 
                if inferenceOptions.automaticBinning
                    nTotal = 0.98*length([trace_struct_filtered(group_ids).fluo]);
                    n_intensity_bins = max([1 round(nTotal/inferenceOptions.SampleSize)]);
                end

                % generate list of average protein levels
                intensity_list_short = [trace_struct_filtered(group_ids).mean_intensity];
                intensity_list_long = [];
                for gInd = 1:length(group_ids)
                  intensity_list_long = [intensity_list_long repelem(trace_struct_filtered(group_ids(gInd)).mean_intensity,...
                                                                          length(trace_struct_filtered(group_ids(gInd)).fluo))];
                end
                
                % generate protein groupings    
                q_vec = linspace(.01,.99,n_intensity_bins+1);        
                intensity_prctile_vec = quantile(intensity_list_long,q_vec);    
                
                % assign traces to groups    
                id_vec = discretize(intensity_list_short,intensity_prctile_vec);                                
                i_val_vec = intensity_prctile_vec(1:end-1) + diff(intensity_prctile_vec)/2; % intensity bin centroids
                for i = 1:length(id_vec)
                    trace_struct_filtered(group_ids(i)).intensity_bin = id_vec(i);
                    if ~isnan(id_vec(i))
                      trace_struct_filtered(group_ids(i)).intensity_bin_val = i_val_vec(id_vec(i));
                    else
                      trace_struct_filtered(group_ids(i)).intensity_bin_val = NaN;
                    end
                end
            end
        end
        % fill in traces that do not fit in any "Additional group" with
        % NaNs
        for i = find(nan_filter1)
          trace_struct_filtered(i).intensity_bin = NaN;
          trace_struct_filtered(i).intensity_bin_val = NaN;
        end
    end
  else
      for i = 1:length(trace_struct_filtered)
          trace_struct_filtered(i).intensity_bin = 1;
          trace_struct_filtered(i).intensity_bin_val = NaN;
      end
  end
  
  
  % remove traces where one more more ID field in NAN
  intensity_group_vec = [trace_struct_filtered.intensity_bin];
  intensity_value_vec = [trace_struct_filtered.intensity_bin_val];
  
  nan_filter = isnan(intensity_group_vec) | nan_filter1;
  trace_struct_filtered = trace_struct_filtered(~nan_filter);
  ap_group_vec = ap_group_vec(~nan_filter);
  time_group_vec = time_group_vec(~nan_filter);
  intensity_group_vec = intensity_group_vec(~nan_filter);
  intensity_value_vec = intensity_value_vec(~nan_filter);
  additional_group_vec = additional_group_vec(~nan_filter);
  
  if inferenceOptions.singleTraceInference
      trace_id_vec = find(~nan_filter);
  else
      trace_id_vec = ones(size(additional_group_vec));
  end
  % generate indexing structure
  indexInfo = struct;  
  [indexInfo.indexVarArray, mapTo, indexInfo.indexList] = unique([ap_group_vec' time_group_vec' intensity_group_vec' additional_group_vec' trace_id_vec'],'rows');
  indexInfo.indexVecUnique = 1:size(indexInfo.indexVarArray,1);
  indexInfo.ap_group_vec = ap_group_vec(mapTo);
  indexInfo.time_group_vec = time_group_vec(mapTo);  
  indexInfo.intensity_group_vec = intensity_group_vec(mapTo);
  indexInfo.intensity_value_vec = intensity_value_vec(mapTo);
  indexInfo.additional_group_vec = additional_group_vec(mapTo);
  indexInfo.varArrayCols = {'AP','Time',inferenceOptions.intensityBinVar,inferenceOptions.AdditionalGroupingVariable};
  
  % update sample size
  inferenceOptions.SampleSize = repelem(inferenceOptions.SampleSize,size(indexInfo.indexVarArray,1));
  