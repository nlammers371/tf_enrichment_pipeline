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
  
  disp('Done.')