function samplingInfo = performNucleusQC(samplingInfo)

  nucleus_stats_raw = regionprops(samplingInfo.nc_ref_frame_raw==1,'MajorAxisLength','MinorAxisLength','Area');
  nc_label_frame_raw = bwlabel(samplingInfo.nc_ref_frame_raw==1);
  
  % perform basic QC
  areaFlags = [nucleus_stats_raw.Area] < samplingInfo.min_nucleus_area | [nucleus_stats_raw.Area] > samplingInfo.max_nucleus_area;
  shapeFlags = [nucleus_stats_raw.MajorAxisLength] >= 1.5*[nucleus_stats_raw.MinorAxisLength];
  keepIndices = find(~areaFlags&~shapeFlags);
  
  % record
  samplingInfo.nc_ref_frame = ismember(nc_label_frame_raw,keepIndices);
  samplingInfo.nc_dist_frame = bwdist(~samplingInfo.nc_ref_frame);
  samplingInfo.nc_label_frame = bwlabel(samplingInfo.nc_ref_frame);
  
  % generate lookup table
  nucleus_stats = regionprops(samplingInfo.nc_ref_frame==1,'Centroid');
  nucleus_centers = vertcat(nucleus_stats.Centroid);
  nc_x_vec = nucleus_centers(:,1)';%RefStruct.nc_x_ref(FrameSetFilter);
  nc_y_vec = nucleus_centers(:,2)';%RefStruct.nc_y_ref(FrameSetFilter);  
  x_dist_mat = repmat(nc_x_vec,numel(nc_x_vec),1)-repmat(nc_x_vec',1,numel(nc_x_vec));
  y_dist_mat = repmat(nc_y_vec,numel(nc_y_vec),1)-repmat(nc_y_vec',1,numel(nc_y_vec));
  
  % record
  samplingInfo.nc_x_vec = nc_x_vec;
  samplingInfo.nc_y_vec = nc_y_vec;
  samplingInfo.r_dist_mat = sqrt(double(x_dist_mat).^2 + double(y_dist_mat).^2); 