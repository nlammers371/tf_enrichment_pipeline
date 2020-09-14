function refVecStruct = generateReferenceVectors(spot_struct_protein,refPath,use3DSpotInfo,ignoreQC)

% initialize structure to store reference vectors 
refVecStruct = struct;
refVecStruct.refPath = refPath;

% All of these vectors have an element for each spot detection in the set
refVecStruct.frame_ref = [spot_struct_protein.frames];
refVecStruct.nc_x_ref = [spot_struct_protein.xPosNucleus];
refVecStruct.nc_y_ref = [spot_struct_protein.yPosNucleus];

if use3DSpotInfo
    refVecStruct.spot_x_ref = [spot_struct_protein.xPosParticle3D]; 
    refVecStruct.spot_y_ref = [spot_struct_protein.yPosParticle3D];
    refVecStruct.spot_z_ref = [spot_struct_protein.zPosParticle3D];
else
  refVecStruct.spot_x_ref = [spot_struct_protein.xPosParticle];
  refVecStruct.spot_y_ref = [spot_struct_protein.yPosParticle];
  refVecStruct.spot_z_ref = [spot_struct_protein.zPosParticle];
end
  
refVecStruct.setID_ref = [];
refVecStruct.master_nucleusID_ref = []; % list of nucleus indices
refVecStruct.particle_index_ref = []; % list of particle indices within nucleus_struct
refVecStruct.particle_subindex_ref = []; % list of frame indices within each particle
refVecStruct.particleID_ref = [];
refVecStruct.particle_qcFlag_ref = [];
for i = 1:numel(spot_struct_protein)
    ParticleID = spot_struct_protein(i).particleID;    
    refVecStruct.particleID_ref = [refVecStruct.particleID_ref repelem(ParticleID, numel(spot_struct_protein(i).frames))];
    refVecStruct.particle_qcFlag_ref = [refVecStruct.particle_qcFlag_ref repelem(spot_struct_protein(i).qcFlag, numel(spot_struct_protein(i).frames))];
    refVecStruct.setID_ref = [refVecStruct.setID_ref repelem(spot_struct_protein(i).setID,numel(spot_struct_protein(i).frames))];
    refVecStruct.master_nucleusID_ref = [refVecStruct.master_nucleusID_ref repelem(round(spot_struct_protein(i).ncID*1e5),numel(spot_struct_protein(i).frames))];
    refVecStruct.particle_index_ref = [refVecStruct.particle_index_ref repelem(i,numel(spot_struct_protein(i).frames))]; 
    refVecStruct.particle_subindex_ref = [refVecStruct.particle_subindex_ref 1:numel(spot_struct_protein(i).frames)];
end

% option to override qc 
if ignoreQC 
    refVecStruct.particle_qcFlag_ref = true(size(refVecStruct.particle_qcFlag_ref));
end

%%% Generate reference array for set-frame combos
refVecStruct.set_frame_array = unique([refVecStruct.setID_ref' refVecStruct.frame_ref'],'row');