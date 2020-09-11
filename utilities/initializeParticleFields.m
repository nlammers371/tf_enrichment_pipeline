function compiledSchnitzCells = initializeParticleFields(compiledSchnitzCells,sWithNC,...
  n_entries)%,threeDFlag,APFlag,DetrendedZFlag)

    compiledSchnitzCells(sWithNC).particleID = NaN;
    compiledSchnitzCells(sWithNC).xPosParticle = NaN(1,n_entries);
    compiledSchnitzCells(sWithNC).yPosParticle = NaN(1,n_entries);
    compiledSchnitzCells(sWithNC).zPosParticle = NaN(1,n_entries);      
    compiledSchnitzCells(sWithNC).spotFrames = NaN(1,n_entries); 
%     if threeDFlag
      compiledSchnitzCells(sWithNC).xPosParticle3D = NaN(1,n_entries);
      compiledSchnitzCells(sWithNC).yPosParticle3D = NaN(1,n_entries);
      compiledSchnitzCells(sWithNC).zPosParticle3D = NaN(1,n_entries); 
      compiledSchnitzCells(sWithNC).fluo3D = NaN(1,n_entries);
%     end
%     if APFlag 
    compiledSchnitzCells(sWithNC).APPosParticle = NaN(1,n_entries);
%     end
%     if DetrendedZFlag
    compiledSchnitzCells(sWithNC).zPosParticleDetrended = NaN(1,n_entries);
%     end
    compiledSchnitzCells(sWithNC).fluo = NaN(1,n_entries); 
    compiledSchnitzCells(sWithNC).fluoOffset = NaN(1,n_entries); 
    compiledSchnitzCells(sWithNC).spotFrames = NaN(1,n_entries); 