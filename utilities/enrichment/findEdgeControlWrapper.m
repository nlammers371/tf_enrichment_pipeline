function [spot_edge_dist, edge_x, edge_y, qc_flag] = findEdgeControlWrapper(...
                                samplingInfo,nucleus_mask,x_index,y_index,nucleus_id)

  
  % Edge sampling 
  spot_edge_dist = samplingInfo.nc_dist_frame(y_index,x_index);        
  nc_edge_dist_vec = samplingInfo.nc_dist_frame(nucleus_mask);  
  spot_sep_vec = samplingInfo.spot_dist_frame(nucleus_mask);        

  % Now find control "spot" that is same distance from nucleus edge
  % as true spot
  [edge_x, edge_y, qc_flag,~]= find_control_sample(...
    nc_edge_dist_vec, samplingInfo, spot_sep_vec, spot_edge_dist, nucleus_mask);  

  % if initial attempt failed, try nearest neighbor nucleus
  if qc_flag == 0

      % Find nearest neighbor nucleus
      r_vec = samplingInfo.r_dist_mat(:,nucleus_id);
      r_vec(nucleus_id) = Inf;
      [min_dist, closestIndex] = min(r_vec);
      
      if min_dist <= samplingInfo.max_dist_nearest_neighbor
          % extract mask 
          nn_x_pos = min([samplingInfo.xDim max([1 round(samplingInfo.nc_x_vec(closestIndex))])]);
          nn_y_pos = min([samplingInfo.yDim max([1 round(samplingInfo.nc_y_vec(closestIndex))])]);
          nn_nucleus_mask_id = samplingInfo.nc_label_frame(nn_y_pos,nn_x_pos);      
          nn_nucleus_mask = samplingInfo.nc_label_frame == nn_nucleus_mask_id;      

          % make sure size is reasonable                     
          nn_edge_dist_vec = samplingInfo.nc_dist_frame(nn_nucleus_mask);
          nn_sep_vec = samplingInfo.spot_dist_frame(nn_nucleus_mask);

          [edge_x, edge_y, qc_flag,~]= find_control_sample(...
              nn_edge_dist_vec, samplingInfo, nn_sep_vec, spot_edge_dist, nn_nucleus_mask);  

          if qc_flag == 1
              qc_flag = 2;
          end
      end
      
  end 
  