function [sample_x, sample_y, qc_flag, sample_edge_distance] = ...
    find_control_sample(...
    edge_dist_vec, samplingInfo, spot_sep_vec, spot_edge_dist ,nucleus_mask)
    
    % initialize variables  
    sample_x = NaN;
    sample_y = NaN;    
    sample_edge_distance = NaN;
    qc_flag = 0;   
    
    % get position vectors for nucleus mask
    x_pos_vec_spot = samplingInfo.x_ref(nucleus_mask);
    y_pos_vec_spot = samplingInfo.y_ref(nucleus_mask);
   
    sample_index_vec = 1:numel(spot_sep_vec);
    
    % find closest pixel that meets criteria
    candidate_pixel_filter = spot_sep_vec >= samplingInfo.minSampleSep & round(edge_dist_vec) == round(spot_edge_dist);
    sample_index_vec = sample_index_vec(candidate_pixel_filter);
    
    % if candidate found, then proceed. Else look to neighboring nuclei
    if ~isempty(sample_index_vec)
        sample_index = randsample(sample_index_vec,1);
        sample_x = x_pos_vec_spot(sample_index);
        sample_y = y_pos_vec_spot(sample_index);        
        sample_edge_distance = round(spot_edge_dist);
        qc_flag = 1;            
    end
%     elseif force_sample
%         new_filter = spot_sep_vec >= min_sample_sep;
%         distances = abs(spot_edge_dist-edge_dist_vec);
%         distances(~new_filter) = inf;
%         [~, sample_index] = min(distances);
%         sample_edge_distance = edge_dist_vec(sample_index);
%         sample_x = x_pos_vec_spot(sample_index);
%         sample_y = y_pos_vec_spot(sample_index);
%         sample_nucleus = index;
%         qc_flag = 3; 
%     end