function spot_struct_protein = mapToSpotStructure(spot_struct_protein,SamplingResults,RefStruct,SetFrameArray)

  disp('mapping parallelized data back to main data structure...')
  SampledFields = fieldnames(SamplingResults);
  
  for i = 1:size(SetFrameArray,1)
    
      % read basic info from set_frame array
      SetID = SetFrameArray(i,1);
      Frame = SetFrameArray(i,2);  
      
      FrameSetIndices = find(RefStruct.setID_ref==SetID&RefStruct.frame_ref==Frame);    
      
      for j = 1:length(FrameSetIndices)
        
        ParticleIndex = RefStruct.particle_index_ref(FrameSetIndices(j));
        ParticleSubIndex = RefStruct.particle_subindex_ref(FrameSetIndices(j));
        
        for sf = 1:length(SampledFields)
          VecTo = spot_struct_protein(ParticleIndex).(SampledFields{sf});
          VecFrom = SamplingResults(i).(SampledFields{sf});
          VecTo(ParticleSubIndex) = VecFrom(j);
          spot_struct_protein(ParticleIndex).(SampledFields{sf}) = VecTo;
        end
      end
  end
  
  % generate interpolated protein fields (keep this limited for now)
  interp_fields = {'nuclear_protein_vec','spot_protein_vec'};
  for i = 1:length(spot_struct_protein)
    timeInterp = spot_struct_protein(i).timeInterp;
    timeRaw = spot_struct_protein(i).time;
    for j = 1:length(interp_fields)
      vec = spot_struct_protein(i).(interp_fields{j});
      nan_ft = ~isnan(vec);
      if sum(nan_ft)>1
          spot_struct_protein(i).([interp_fields{j} 'Interp']) = interp1(timeRaw(nan_ft),vec(nan_ft),timeInterp);
%       elseif sum(nan_ft)==1
%           spot_struct_protein(i).([interp_fields{j} 'Interp']) = vec(nan_ft);
      else
          spot_struct_protein(i).([interp_fields{j} 'Interp']) = NaN(size(timeInterp));
      end
    end
  end
  disp('Done.')