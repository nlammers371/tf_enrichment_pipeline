function [serial_control_x, serial_control_y, serial_edge_dist] = drawSerializedControlSpot(spot_struct_protein,spotIndex,spotSubIndex,...
                                  tempParams,nc_dist_frame,spot_dist_frame,spot_nc_mask,driftTol)

  % spot and nucleus distance info
  spot_edge_dist = nc_dist_frame(y_index,x_index);        
  nc_edge_dist_vec = nc_dist_frame(spot_nc_mask);
  spot_struct_protein(spotIndex).spot_edge_dist_vec(spotSubIndex) = spot_edge_dist;        
  spot_sep_vec = spot_dist_frame(spot_nc_mask); 
                                
  % time series info
  frame_vec_temp = spot_struct_protein(spotIndex).frames; 
  currentFrame = frame_vec_temp(spotSubIndex);
  serial_null_x = spot_struct_protein(spotIndex).serial_null_x_vec;
  serial_null_y = spot_struct_protein(spotIndex).serial_null_y_vec; 

  % if this is the first sample for this spot, just find random
  % control snip. This will "seed" subsequent samples
  if all(isnan(serial_null_x))            
      % Take a random sample filter for regions far enough away from locus
      sample_index_vec = 1:numel(spot_sep_vec);          
      criteria_filter = spot_sep_vec >= tempParams.minSampleSep & nc_edge_dist_vec >= tempParams.minEdgeSep;
      sample_index_vec = sample_index_vec(criteria_filter);

      % if candidate found, then proceed. Else look to neighboring nuclei
      if ~isempty(sample_index_vec)
          if length(sample_index_vec) > 1
            new_index = randsample(sample_index_vec,1);
          else
            new_index = sample_index_vec;
          end
          x_pos_vec = x_ref(spot_nc_mask);
          y_pos_vec = y_ref(spot_nc_mask);
          serial_control_x = x_pos_vec(new_index);
          serial_control_y = y_pos_vec(new_index);
          serial_edge_dist = nc_edge_dist_vec(new_index);                           
      else
          error('Unable to draw random sample. Check "PixelSize" and "min_sample_sep" variables')
      end  
  % otherwise, draw snip based on previous location
  else
      prev_frame = find(~isnan(serial_null_x),1,'last');
      n_frames = currentFrame - frame_vec_temp(prev_frame); % used to adjust jump weights
      old_x = double(serial_null_x(prev_frame));
      old_y = double(serial_null_y(prev_frame));    

      % possible locations
      x_pos_vec = x_ref(spot_nc_mask);
      y_pos_vec = y_ref(spot_nc_mask);
      drControl = double(sqrt((old_x-x_pos_vec).^2+(old_y-y_pos_vec).^2));   

      % calculate weights
      wt_vec = exp(-.5*((drControl/double(sqrt(n_frames)*(driftTol))).^2));

      % anything too close to locus or with an edge distance too different from locus is excluded
      wt_vec(spot_sep_vec<minSampleSep|nc_edge_dist_vec < minEdgeSep) = 0;

      % draw sample
      if any(wt_vec>0)
          new_index = randsample(1:numel(x_pos_vec),1,true,wt_vec);
          serial_control_x = x_pos_vec(new_index);
          serial_control_y = y_pos_vec(new_index);
          serial_edge_dist = nc_edge_dist_vec(new_index);       
      else
          new_index = randsample(1:numel(x_pos_vec),1,true);
          serial_control_x = x_pos_vec(new_index);
          serial_control_y = y_pos_vec(new_index);
          serial_edge_dist = nc_edge_dist_vec(new_index);  
      end
  end   