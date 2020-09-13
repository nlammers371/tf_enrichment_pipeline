function [null_x, null_y, null_nc, qc_flag, null_e] = find_control_sample(...
    edge_dist_vec, x_ref, y_ref, spot_sep_vec, spot_edge_dist, index, min_sample_sep,null_mask,...
    force_sample)
    
    % initialize variables  
    null_x = NaN;
    null_y = NaN;
    null_nc = NaN;
    null_e = NaN;
    qc_flag = 0;   
    
    % get position vectors for nucleus mask
    x_pos_vec_spot = x_ref(null_mask);
    y_pos_vec_spot = y_ref(null_mask);
   
    sample_index_vec = 1:numel(spot_sep_vec);
    % find closest pixel that meets criteria
    cr_filter = spot_sep_vec >= min_sample_sep & round(edge_dist_vec) == round(spot_edge_dist);
    sample_index_vec = sample_index_vec(cr_filter);
    % if candidate found, then proceed. Else look to neighboring nuclei
    if ~isempty(sample_index_vec)
        sample_index = randsample(sample_index_vec,1);
        null_x = x_pos_vec_spot(sample_index);
        null_y = y_pos_vec_spot(sample_index);
        null_nc = index;
        null_e = round(spot_edge_dist);
        qc_flag = 1;               
    elseif force_sample
        new_filter = spot_sep_vec >= min_sample_sep;
        distances = abs(spot_edge_dist-edge_dist_vec);
        distances(~new_filter) = inf;
        [~, sample_index] = min(distances);
        null_e = edge_dist_vec(sample_index);
        null_x = x_pos_vec_spot(sample_index);
        null_y = y_pos_vec_spot(sample_index);
        null_nc = index;
        qc_flag = 1; 
    end