function nuclearSegmentation(liveProject, refVecStruct, segmentIndices, NumWorkers)  

  p = gcp('nocreate');
  if isempty(p) || p.NumWorkers~=NumWorkers
    parpool(NumWorkers);
  end
  % initialize arrays to store segmentation info
  nucleus_frame_cell = cell(1,length(segmentIndices));
  spot_frame_cell = cell(1,length(segmentIndices));
  parfor w = 1:length(segmentIndices)
      i = segmentIndices(w);
      currentSetID = refVecStruct.set_frame_array(i,1);
      currentFrame = refVecStruct.set_frame_array(i,2);  
      
      % get nucleus
      frame_set_filter = refVecStruct.setID_ref==currentSetID&refVecStruct.frame_ref==currentFrame;
      nc_x_vec_temp = round(refVecStruct.nc_x_ref(frame_set_filter));
      nc_y_vec_temp = round(refVecStruct.nc_y_ref(frame_set_filter)); 
      
      % indexing vectors    
      nc_master_vec = refVecStruct.master_nucleusID_ref(frame_set_filter);  
      
      % get list of unique indices 
      [nc_master_vec_u,ia,~] = unique(nc_master_vec,'stable');
      
      % unique nucleus vectors
      nc_x_vec_u = round(nc_x_vec_temp(ia));
      nc_y_vec_u = round(nc_y_vec_temp(ia));    
      
      % particle positions        
      spot_x_vec = round(refVecStruct.spot_x_ref(frame_set_filter));
      spot_y_vec = round(refVecStruct.spot_y_ref(frame_set_filter));            

      % Get experiment info
      Prefix = liveProject.includedExperimentNames{currentSetID}; 
      currExperiment = liveProject.includedExperiments{currentSetID};                  
      proteinChannel = currExperiment.inputChannels;
      if length(proteinChannel) > 1
        error(['Problem with Prefix: ' Prefix '. This pipeline does not currently support multiple input channels'])
      end

      % Generate reference vectors
      xDim = currExperiment.xDim;
      yDim = currExperiment.yDim;
      zDim = currExperiment.zDim;        

      % git pixel size info        
      PixelSize = currExperiment.pixelSize_nm / 1e3; % convert to microns
      nucleus_neighborhood_size = round(10 ./ PixelSize);  % determine size of neighborhood to use during nucleus segmentation      
      smoothing_kernel_size = round(1 ./ PixelSize); % size of gaussian smoothing kernel 

      % get protein channel       
      stackPath = [currExperiment.preFolder  Prefix '_' sprintf('%03d',currentFrame) '_ch0' num2str(proteinChannel) '.tif'];
      protein_stack = imreadStack2(stackPath, yDim, xDim, zDim+2);      
      protein_stack = protein_stack(:,:,2:end-1);
      
      % generate protein gradient frame for segmentation
      protein_smooth = imgaussfilt(mean(protein_stack,3),round(smoothing_kernel_size/2));                
      protein_grad = imgradient(protein_smooth);   
      
      % flatten background         
      protein_bkg = imgaussfilt(protein_smooth, round(nucleus_neighborhood_size/2));
      protein_grad_norm = protein_grad ./ protein_bkg;        

      % try local thresholding
      n_x_local = floor(xDim / nucleus_neighborhood_size);
      n_y_local = floor(yDim / nucleus_neighborhood_size);
      x_dim_local = round(xDim / n_x_local);
      y_dim_local = round(yDim / n_y_local);
      protein_bin_clean = false(size(protein_grad_norm));
      for x = 1:n_x_local
          for y = 1:n_y_local
              x_start = (x-1)*x_dim_local+1;
              x_stop = min([x*x_dim_local,xDim]);
              y_start = (y-1)*y_dim_local+1;
              y_stop = min([y*y_dim_local,yDim]);
              pt_section = protein_grad_norm(y_start:y_stop,x_start:x_stop);
              thresh = multithresh(pt_section);
              section_bin = pt_section > thresh; 
              section_bin_clean = bwareaopen(section_bin,smoothing_kernel_size^2);                
              % record
              protein_bin_clean(y_start:y_stop,x_start:x_stop) = bwmorph(section_bin_clean,'hbreak');
          end
      end        
      % label regions
      nc_frame_labels = logical(protein_bin_clean);         
      % frame info
      nc_lin_indices = sub2ind(size(protein_bin_clean),round(nc_y_vec_u),round(nc_x_vec_u));
      % take convex hull
      stats = regionprops(nc_frame_labels,'ConvexHull');
      nc_ref_frame = zeros(size(nc_frame_labels)); 
      for j = 2:numel(stats)
          hull_points = stats(j).ConvexHull;
          mask = poly2mask(hull_points(:,1),hull_points(:,2),yDim,xDim);   
          nc_bin_ids = mask(nc_lin_indices);
          if sum(nc_bin_ids) == 1 % enforce unique
              nc_ref_frame(mask) = nc_master_vec_u(nc_bin_ids);
          end
      end       
      % generate array indicating distance of each pixel from an active locus 
      nc_indices = sub2ind(size(nc_ref_frame),spot_y_vec,spot_x_vec);
      spot_dist_frame_temp = zeros(size(nc_ref_frame));
      spot_dist_frame_temp(nc_indices(~isnan(nc_indices))) = 1;
      spot_dist_frame_temp = bwdist(spot_dist_frame_temp);
      % label regions within designated integration radius of a spot    
      % store arrays
      nucleus_frame_cell{w} = nc_ref_frame;
      spot_frame_cell{w} = spot_dist_frame_temp;      
  end
  disp('saving segmentation results...')
  % save arrays
  for w = 1:numel(segmentIndices)
      i = segmentIndices(w);
      currentSetID = refVecStruct.set_frame_array(i,1);
      currentFrame = refVecStruct.set_frame_array(i,2);  
      nc_ref_name = [refVecStruct.refPath 'nc_ref_frame_set' sprintf('%02d',currentSetID) '_frame' sprintf('%03d',currentFrame) '.mat'];
      nc_ref_frame = nucleus_frame_cell{w};
      save(nc_ref_name,'nc_ref_frame');
      spot_ref_name = [refVecStruct.refPath 'spot_roi_frame_set' sprintf('%02d',currentSetID) '_frame' sprintf('%03d',currentFrame) '.mat'];
      spot_dist_frame = spot_frame_cell{w};
      save(spot_ref_name,'spot_dist_frame');
  end
  disp('done.')
 