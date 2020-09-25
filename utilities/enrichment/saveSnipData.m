function saveSnipData(samplingInfo,TempSnipStruct,RefStruct,liveProject)

    % initialize struct to store snip data
    snip_data = struct;        
    % get fields
    snip_fields = fieldnames(TempSnipStruct);
    % add protein snip stacks    
    for j = 1:length(snip_fields)
        snip_data.(snip_fields{j}) = TempSnipStruct.(snip_fields{j});
    end    

    % store key ID variables
    snip_data.frame = samplingInfo.Frame;
    snip_data.setID = samplingInfo.SetID;
    snip_data.particle_id_vec = RefStruct.particleID_ref(samplingInfo.FrameSetFilter);

    % indexing vectors    
    snip_data.nc_sub_index_vec = RefStruct.particle_subindex_ref(samplingInfo.FrameSetFilter);
    snip_data.nc_lin_index_vec = RefStruct.particle_index_ref(samplingInfo.FrameSetFilter); 
    snip_data.nc_master_vec = RefStruct.master_nucleusID_ref(samplingInfo.FrameSetFilter);    

    % specify name and save
    snip_name = ['snip_data_F' sprintf('%03d',samplingInfo.Frame) '_S' sprintf('%02d',samplingInfo.SetID)];     
    save([liveProject.dataPath '/snip_fragments_temp/' snip_name '.mat'],'snip_data')
    