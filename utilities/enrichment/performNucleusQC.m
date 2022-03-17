function samplingInfo = performNucleusQC(samplingInfo)

  nucleus_stats_raw = regionprops(samplingInfo.nc_ref_frame_raw>0,'MajorAxisLength','MinorAxisLength','Area');
  
%   idList = unique(samplingInfo.nc_ref_frame_raw);
%   idList = idList(idList~=0);
%   areaList = NaN(size(idList));
%   for i = 2:length(idList)
%     areaList(i-1) = sum(samplingInfo.nc_ref_frame_raw(:)==idList(i));
%   end
    
  % perform basic QC
  areaFlags = [nucleus_stats_raw.Area] < samplingInfo.min_nucleus_area | [nucleus_stats_raw.Area] > samplingInfo.max_nucleus_area;%
  shapeFlags = [nucleus_stats_raw.MajorAxisLength] >= 2*[nucleus_stats_raw.MinorAxisLength];
  keepIndices = find(~areaFlags&~shapeFlags);
  
  % record binary mask and labeled mask 
  samplingInfo.nc_ref_frame = ismember(bwlabel(samplingInfo.nc_ref_frame_raw>0),keepIndices);
  samplingInfo.nc_label_frame = bwlabel(samplingInfo.nc_ref_frame);%cast(samplingInfo.nc_ref_frame,'uint8').*samplingInfo.nc_ref_frame_raw;
%   samplingInfo.nc_label_frame(samplingInfo.nc_ref_frame) = samplingInfo.nc_ref_frame_raw(samplingInfo.nc_ref_frame);
  
  % shave off edges for distance calculation to avoid occasional instances
  % where two nucleus masks are directly adjacent
%   se = strel('disk',1);
%   eroded_frame = imerode(samplingInfo.nc_ref_frame,se);  
  samplingInfo.nc_dist_frame = bwdist(~samplingInfo.nc_ref_frame);
       

  % generate lookup table
  nucleus_stats = regionprops(samplingInfo.nc_ref_frame==1,'Centroid');
  nucleus_centers = vertcat(nucleus_stats.Centroid);
  nc_x_vec = [];
  nc_y_vec = [];
  if ~isempty(nucleus_centers)
      nc_x_vec = nucleus_centers(:,1)';%RefStruct.nc_x_ref(FrameSetFilter);
      nc_y_vec = nucleus_centers(:,2)';%RefStruct.nc_y_ref(FrameSetFilter);  
  end
  x_dist_mat = repmat(nc_x_vec,numel(nc_x_vec),1)-repmat(nc_x_vec',1,numel(nc_x_vec));
  y_dist_mat = repmat(nc_y_vec,numel(nc_y_vec),1)-repmat(nc_y_vec',1,numel(nc_y_vec));
  
  % record
  samplingInfo.nc_x_vec = nc_x_vec;
  samplingInfo.nc_y_vec = nc_y_vec;
  samplingInfo.r_dist_mat = sqrt(double(x_dist_mat).^2 + double(y_dist_mat).^2); 