function [new_edge_dist, new_serial_x, new_serial_y,qc_flag] = ...
                                drawSerializedControlSpot(samplingInfo,nucleus_mask,...
                                FrameVec, serial_null_x_vec, serial_null_y_vec)                                  

  % spot and nucleus distance info  
  nc_edge_dist_vec = samplingInfo.nc_dist_frame(nucleus_mask);
  
%   spot_struct_protein(spotIndex).spot_edge_dist_vec(spotSubIndex) = spot_edge_dist;        
  spot_sep_vec = samplingInfo.spot_dist_frame(nucleus_mask); 
  x_pos_vec = samplingInfo.x_ref(nucleus_mask);
  y_pos_vec = samplingInfo.y_ref(nucleus_mask);         
  CurrentFrame = samplingInfo.Frame;         
  
%   CurrentIndex = find(FrameVec == CurrentFrame);
  % generate sampling vector      
  sample_index_vec = find(spot_sep_vec >= samplingInfo.minSampleSep & ...
              nc_edge_dist_vec >= samplingInfo.minEdgeSep);
  
  % if this is the first sample for this spot, just find random
  % control snip. This will "seed" subsequent samples
  qc_flag = 1;
  if ~isempty(sample_index_vec)
      if all(isnan(serial_null_x_vec))          
    
          % Take a random sample filter for regions far enough away from locus    
      
          if length(sample_index_vec) > 1
            new_index = randsample(sample_index_vec,1);
          else
            new_index = sample_index_vec;
          end
          
          new_serial_x = x_pos_vec(new_index);
          new_serial_y = y_pos_vec(new_index);
          new_edge_dist = nc_edge_dist_vec(new_index);                             
            
    % otherwise, draw snip based on previous location
    else
        prevIndex = (find(~isnan(serial_null_x_vec),1,'last'));
        n_frames = double(CurrentFrame) - double(FrameVec(prevIndex)); % used to adjust jump weights
        old_x = double(serial_null_x_vec(prevIndex));
        old_y = double(serial_null_y_vec(prevIndex));    

        % calculate distance from previous location     
        drControl = double(sqrt((old_x-x_pos_vec(sample_index_vec)).^2+(old_y-y_pos_vec(sample_index_vec)).^2));   

        % calculate weights
        wt_vec = exp(-.5*((drControl/double(sqrt(n_frames)*(samplingInfo.driftTol))).^2));      

        % draw sample
        if any(wt_vec>0)
            if length(sample_index_vec)>1
                new_index = randsample(sample_index_vec,1,true,wt_vec);
            else
                new_index = sample_index_vec;
            end
            new_serial_x = x_pos_vec(new_index);
            new_serial_y = y_pos_vec(new_index);
            new_edge_dist = nc_edge_dist_vec(new_index);  
            
        else
            qc_flag = 0;
            new_index = randsample(sample_index_vec,1,true);
            new_serial_x = x_pos_vec(new_index);
            new_serial_y = y_pos_vec(new_index);
            new_edge_dist = nc_edge_dist_vec(new_index);  
        end
      end
  else
      warning('Unable to draw serial control spot. Check nucleus segmentation, and "minEdgeSep" and "minSampSep" parameters')
      qc_flag = 0;
      new_edge_dist = NaN;
      new_serial_x = NaN;
      new_serial_y = NaN;
  end
  
  % record info
%   serial_null_x_vec(samplingInfo.spotSubIndex) = new_serial_x;
% 	serial_null_y_vec(samplingInfo.spotSubIndex) = new_serial_y;
% 	serial_null_edge_dist_vec(samplingInfo.spotSubIndex) = new_edge_dist;