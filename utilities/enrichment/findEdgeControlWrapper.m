function [spot_struct_protein, null_mask] = findEdgeControlWrapper(spot_struct_protein,refVecStruct,x_index,y_index,spot_nc_mask,...
              spotIndex,spotSubIndex,r_dist_mat,frame_set_indices,j_pass)


  % Edge sampling 
  spot_edge_dist = nc_dist_frame(y_index,x_index);        
  nc_edge_dist_vec = nc_dist_frame(spot_nc_mask);
  spot_struct_protein(spotIndex).spot_edge_dist_vec(spotSubIndex) = spot_edge_dist;        
  spot_sep_vec = spot_dist_frame(spot_nc_mask);        

  % Now find control "spot" that is same distance from nucleus edge
  % as true spot
  [spot_struct_protein(spotIndex).edge_null_x_vec(spotSubIndex), spot_struct_protein(spotIndex).edge_null_y_vec(spotSubIndex),...
   spot_struct_protein(spotIndex).edge_null_nc_vec(spotSubIndex), spot_struct_protein(spotIndex).edge_qc_flag_vec(spotSubIndex),~]...
      = find_control_sample(nc_edge_dist_vec, x_ref, y_ref, spot_sep_vec, spot_edge_dist,...
           j_pass, tempParams.minSampleSep, spot_nc_mask,0);  

  % if initial attempt failed, try nearest neighbor nucleus
  null_mask = spot_nc_mask;
  if spot_struct_protein(spotIndex).edge_qc_flag_vec(spotSubIndex) == 0

      % Find nearest neighbor nucleus
      r_vec = r_dist_mat(:,j_pass);
      r_vec(j_pass) = Inf;
      [~, mi] = min(r_vec);

      % get nn nucleus mask   
      x_spot_nn = refVecStruct.spot_x_ref(frame_set_indices(mi));
      y_spot_nn = refVecStruct.spot_y_ref(frame_set_indices(mi)); 

      nn_nc_mask = nc_ref_frame == refVecStruct.master_nucleusID_ref(frame_set_indices(mi));

      null_mask = nn_nc_mask; % reassign null mask
      if ~isnan(x_spot_nn)
          nan_flag = isnan(nn_nc_mask(y_spot_nn,x_spot_nn));
      end

      % make sure size is reasonable 
      if sum(nn_nc_mask(:)) >= min_nucleus_area && sum(nn_nc_mask(:)) <= max_nucleus_area && ~nan_flag                

          nn_edge_dist_vec = nc_dist_frame(nn_nc_mask);
          nn_sep_vec = spot_dist_frame(nn_nc_mask);

          [spot_struct_protein(spotIndex).edge_null_x_vec(spotSubIndex), spot_struct_protein(spotIndex).edge_null_y_vec(spotSubIndex),...
           spot_struct_protein(spotIndex).edge_null_nc_vec(spotSubIndex), spot_struct_protein(spotIndex).edge_qc_flag_vec(spotSubIndex),~]...
              = find_control_sample(nn_edge_dist_vec, x_ref, y_ref, nn_sep_vec, spot_edge_dist,...
                   mi, tempParams.minSampleSep, null_mask,1);
          if spot_struct_protein(spotIndex).edge_qc_flag_vec(spotSubIndex)== 1
            spot_struct_protein(spotIndex).edge_qc_flag_vec(spotSubIndex) = 2;
          end
      end
  end 