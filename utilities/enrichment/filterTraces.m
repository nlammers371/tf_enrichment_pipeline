function [trace_struct_filtered, indexInfo] = filterTraces(inferenceOptions,analysis_traces)
  % this function filters traces and generates grouping variables as
  % appropriate
  % Supported grouping variables include: time, APPositionParticle, and nuclear_protein_vec
    
  % apply QC filter calculated in main01
  analysis_traces = analysis_traces([analysis_traces.qc_flag]==1);
  
  % generate stripped down data structure
  trace_struct_filtered = struct;
  
  i_pass = 1;
  for i = 1:length(analysis_traces)    
      if inferenceOptions.fluo3DFlag
          fluo = analysis_traces(i).fluo3D_interp;
      else
          fluo = analysis_traces(i).fluo_interp;
      end
      time = analysis_traces(i).time_interp;    
      if ~isfield(inferenceOptions,'Tres')
        inferenceOptions.Tres = time(2)-time(1);
      end
      time_raw = analysis_traces(i).time;
      ap_raw = analysis_traces(i).APPosParticle;
      
      for a = 1:length(inferenceOptions.APBins)-1 % Note: if either option is no activated, there will be only one bin
          apBounds = inferenceOptions.APBins(a:a+1);

          for t = 1:length(inferenceOptions.timeBins)-1
            
              start_time = time(1);
              if inferenceOptions.truncInference(t)
                start_time = time(inferenceOptions.nSteps+1);
              end

              timeBounds = inferenceOptions.timeBins(t:t+1);

              ap_time_filter_raw = time_raw >= timeBounds(1) & time_raw < timeBounds(2) & time_raw >= start_time &...
                               ap_raw >= apBounds(1) & ap_raw < apBounds(2);

              nRaw = sum(ap_time_filter_raw);

              if nRaw >= inferenceOptions.minDP
                
                time_vec_temp = time_raw(ap_time_filter_raw);
                time_filter_interp = time>=time_vec_temp(1) & time<-time_vec_temp(end);
                
                % generate new entry    
                trace_struct_filtered(i_pass).fluo = fluo(time_filter_interp);
                trace_struct_filtered(i_pass).time = time(time_filter_interp);
                if inferenceOptions.ProteinBinFlag          
                    trace_struct_filtered(i_pass).mf_protein = nanmean(analysis_traces(i).nuclear_protein_vec(ap_time_filter_raw));
                end       
                trace_struct_filtered(i_pass).ParticleID = analysis_traces(i).ParticleID;  
                trace_struct_filtered(i_pass).N = length(fluo);    
                trace_struct_filtered(i_pass).apBin = a;
                trace_struct_filtered(i_pass).timeBin = t;
                trace_struct_filtered(i_pass).interp_filter = time_filter_interp;
                trace_struct_filtered(i_pass).raw_filter = ap_time_filter_raw;

                % increment
                i_pass = i_pass + 1;
              end
          end
      end
  end
  
  % make new indexing vectors
  ap_group_vec = [trace_struct_filtered.apBin];
  time_group_vec = [trace_struct_filtered.timeBin];
  
  
  % perform protein binning if appropriate
  if inferenceOptions.ProteinBinFlag 
    
    for a = 1:length(inferenceOptions.APBins)-1 % Note: if either option is no activated, there will be only one bin                
        for t = 1:length(inferenceOptions.timeBins)-1          

            ap_time_ids = find(ap_group_vec==a & time_group_vec==t);          
            % estimate number of bins 
            if inferenceOptions.automaticBinning
                nTotal = sum([trace_struct_filtered(ap_time_ids).N]);
                n_protein_bins = ceil(nTotal/inferenceOptions.SampleSize);
            end

            % generate list of average protein levels
            mf_list = [trace_struct_filtered(ap_time_ids).mf_protein];
            % generate protein groupings    
            q_vec = linspace(0,1,n_protein_bins+1);        
            nuclear_protein_prctile_vec = quantile(mf_list,q_vec);    
            % assign traces to groups    
            id_vec = discretize(mf_list,nuclear_protein_prctile_vec);

            for i = 1:length(id_vec)
                trace_struct_filtered(ap_time_ids(i)).nuclear_protein_bin = id_vec(i);
                trace_struct_filtered.nuclear_protein_quantiles = nuclear_protein_prctile_vec;
            end
        end
    end
  end
  
  % generate indexing structure
  indexInfo = struct;
  nuclear_protein_vec = [trace_struct_filtered.nuclear_protein_bin];
  [indexInfo.indexVarArray, mapTo, indexInfo.indexList] = unique([ap_group_vec' time_group_vec' nuclear_protein_vec'],'rows');
  indexInfo.indexVec = 1:size(indexInfo.indexVarArray,1);
  indexInfo.ap_group_vec = ap_group_vec(mapTo);
  indexInfo.time_group_vec = time_group_vec(mapTo);  
  indexInfo.protein_groupd_vec = nuclear_protein_vec(mapTo);
  