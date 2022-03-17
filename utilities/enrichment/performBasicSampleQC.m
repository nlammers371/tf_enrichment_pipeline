function samplingChecks = performBasicSampleQC(spot_struct_protein)
    
    % set thresholds for throwing various warnings
    missThresh = 0.01;
    outsideThresh = 0.01;
    nnThresh = 0.05;
    
    samplingChecks = struct;
    
    % initialize QC flags
    samplingChecks.edgeMissQC = 0;
    samplingChecks.edgeNNQC = 0;
    samplingChecks.serialMissQC = 0;
    samplingChecks.SpotOutsideQC = 0;
    
    % calculate edge control-related stats
    samplingChecks.totalEdgeControlSamples = sum(0<[spot_struct_protein.edge_qc_flag_vec]);
    samplingChecks.fracEdgeControlSamples = mean(0<[spot_struct_protein.edge_qc_flag_vec]);
    samplingChecks.fracEdgeMisses = mean(0==[spot_struct_protein.edge_qc_flag_vec]);
    samplingChecks.fracNearestNeighbor = sum(2==[spot_struct_protein.edge_qc_flag_vec])/samplingChecks.totalEdgeControlSamples;
    
    % serial control-related stats
    samplingChecks.totalSerialControlSamples = sum(0<[spot_struct_protein.serial_qc_flag_vec]);
    samplingChecks.fracSerialControlSamples = mean(0<[spot_struct_protein.serial_qc_flag_vec]);
    samplingChecks.fracSerialMisses = mean(0==[spot_struct_protein.edge_qc_flag_vec]);
    
    % check how many spots fell outside of a nucleus mask 
    samplingChecks.fracSpotsOutsideNuclei = mean(-1==([spot_struct_protein.serial_qc_flag_vec]));
    
    % throw flags if appropriate
    if samplingChecks.fracEdgeMisses > missThresh
      prctMiss = round(100*samplingChecks.fracEdgeMisses,3);
      warning(['Unusually large percentage (' num2str(prctMiss) '%) of edge control spots were not assigned. Check Nucleus segmentation quality, as well as the size of "minSampleSep_um"'])
      samplingChecks.edgeMissQC = 1;    
    end
    
    if samplingChecks.fracSerialMisses > missThresh
      prctMiss = round(100*samplingChecks.fracSerialMisses,3);
      warning(['Unusually large percentage (' num2str(prctMiss) '%) of serial control spots were not assigned. Check Nucleus segmentation quality, as well as the size of "minSampleSep_um"'])
      warning('Also check that "drifTol" parameter calculation in "calculateVirtualSpotDrif" function is working properly')      
      samplingChecks.serialMissQC = 1;      
    end
    
    if samplingChecks.fracNearestNeighbor > nnThresh
      samplingChecks.edgeNNQC = 1;
      prctNN = round(100*samplingChecks.fracNearestNeighbor ,3);
      warning(['Unusually large number (' num2str(prctNN) '%) of edge control spots could not be drawn from same nucleus as spot and were instead drawn from neighboring nuclei. Check Nucleus segmentation quality, as well as the size of "minSampleSep_um"'])
    end
    
    if samplingChecks.fracSpotsOutsideNuclei > outsideThresh
       prctOut = round(100*samplingChecks.fracSpotsOutsideNuclei,3);
       warning(['Unusually large number (' num2str(prctOut) '%) of spots were found to be located outside nuclear masks. Check Nucleus segmentation quality, and make sure the "min_nucleus_radius_um" and "max_nucleus_radius_um" are set appropriately for current dataype'])
       samplingChecks.SpotOutsideQC = 1;
    end