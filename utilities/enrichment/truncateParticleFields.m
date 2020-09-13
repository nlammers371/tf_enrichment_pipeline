function nucleus_struct = truncateParticleFields(nucleus_struct,has3DSpotInfo,hasAPInfo)

    nucleus_struct = nucleus_struct(~isnan([nucleus_struct.particleID]));
    nucleus_struct = rmfield(nucleus_struct,'spotFrames');
    
    for i = 1:length(nucleus_struct)
      fluo = nucleus_struct(i).fluo;
      nanFilter = ~isnan(fluo);
      
      nucleus_struct(i).fluo = fluo(nanFilter);
      nucleus_struct(i).time = nucleus_struct(i).time(nanFilter);
      nucleus_struct(i).frames = nucleus_struct(i).frames(nanFilter);
      
      nucleus_struct(i).xPosParticle = nucleus_struct(i).xPosParticle(nanFilter);
      nucleus_struct(i).yPosParticle = nucleus_struct(i).yPosParticle(nanFilter);
      nucleus_struct(i).zPosParticle = nucleus_struct(i).zPosParticle(nanFilter);      
      
      if has3DSpotInfo
        nucleus_struct(i).xPosParticle3D = nucleus_struct(i).xPosParticle3D(nanFilter);
        nucleus_struct(i).yPosParticle3D = nucleus_struct(i).yPosParticle3D(nanFilter);
        nucleus_struct(i).yPosParticle3D = nucleus_struct(i).zPosParticle3D(nanFilter); 
        nucleus_struct(i).fluo3D = nucleus_struct(i).fluo3D(nanFilter);
      end
      if hasAPInfo 
        nucleus_struct(i).APPosParticle = nucleus_struct(i).APPosParticle(nanFilter);
      end

      nucleus_struct(i).zPosParticleDetrended = nucleus_struct(i).zPosParticleDetrended(nanFilter);
      
      nucleus_struct(i).fluoOffset = nucleus_struct(i).fluoOffset(nanFilter);       
    end