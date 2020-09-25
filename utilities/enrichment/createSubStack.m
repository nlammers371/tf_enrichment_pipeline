function samplingSubInfo = createSubStack(samplingInfo,samplingSubInfo,nucleus_mask_id)

    % create mask          
    nucleus_mask_2D = samplingInfo.nc_label_frame == nucleus_mask_id;        
    samplingSubInfo.nucleus_mask_3D_full = repmat(nucleus_mask_2D,1,1,samplingInfo.zDim); 
    samplingSubInfo.nucleus_mask_3D_full(:,:,~ismember(samplingSubInfo.z_range3_full,samplingSubInfo.z_range3)) = false;

    % find bounding box for mask and use this to extract compact
    % area of interest to sample
    nucleus_mask_stats = regionprops3(samplingSubInfo.nucleus_mask_3D_full,'BoundingBox','SubarrayIdx');
    samplingSubInfo.protein_stack = samplingInfo.protein_stack(nucleus_mask_stats.SubarrayIdx{:});
    samplingSubInfo.mcp_stack = samplingInfo.mcp_stack(nucleus_mask_stats.SubarrayIdx{:});
    samplingSubInfo.nucleus_mask_3D = samplingSubInfo.nucleus_mask_3D_full(nucleus_mask_stats.SubarrayIdx{:});

    % Now adjust spot location 
    samplingSubInfo.x_shift =  - nucleus_mask_stats.BoundingBox(1) + 0.5;
    samplingSubInfo.x_spot = samplingSubInfo.x_spot_full + samplingSubInfo.x_shift;
    samplingSubInfo.y_shift =  - nucleus_mask_stats.BoundingBox(2) + 0.5;
    samplingSubInfo.y_spot = samplingSubInfo.y_spot_full + samplingSubInfo.y_shift;
    samplingSubInfo.z_shift =  - nucleus_mask_stats.BoundingBox(3) + 0.5;
    samplingSubInfo.z_spot = samplingSubInfo.z_spot_full + samplingSubInfo.z_shift;

    % create adjusted reference arrays and variables
    samplingSubInfo.yDim = size(samplingSubInfo.mcp_stack,1);
    samplingSubInfo.xDim = size(samplingSubInfo.mcp_stack,2);
    samplingSubInfo.zDim = size(samplingSubInfo.mcp_stack,3);
    [samplingSubInfo.x_ref,samplingSubInfo.y_ref,samplingSubInfo.z_ref] = ...
      meshgrid(1:samplingSubInfo.xDim,1:samplingSubInfo.yDim,1:samplingSubInfo.zDim);

    samplingSubInfo.xy_sigma = samplingInfo.xy_sigma;
    samplingSubInfo.z_sigma = samplingInfo.z_sigma;
    samplingSubInfo.snippet_size = samplingInfo.snippet_size;
    
    % get 2D nucleus centroid
    stats2D = regionprops(nanmax(samplingSubInfo.nucleus_mask_3D,[],3),'Centroid');
    samplingSubInfo.x_nucleus = stats2D(1).Centroid(1);
    samplingSubInfo.y_nucleus = stats2D(1).Centroid(2);