function [null_x, null_y, null_nc, qc_flag, null_mask,dist_mat_nn] = ...
    find_control_sample(dist_vec, x_ref, y_ref, x_spot, y_spot, spot_dist,...
    nc_x_vec, nc_y_vec, spot_x_vec, spot_y_vec, index, null_mask_spot, ...
    r_dist_mat, min_sample_sep, histone_image,id_array, nb_sz, nc_index_vec, ...
    area_bounds, centroid_flag)
    
    % dist_mat_nn defaults to zero matrix
    dist_mat_nn = zeros(size(x_ref));
    % extract area info
    min_area = area_bounds(1);
    max_area = area_bounds(2);
    % first try taking a sample from within same nucleus
    null_mask = null_mask_spot;
    % get position vectors for nucleus mask
    x_pos_vec_spot = x_ref(null_mask_spot);
    y_pos_vec_spot = y_ref(null_mask_spot);
    % calculate distance from spot
    x_sep_vec = x_pos_vec_spot - x_spot;        
    y_sep_vec = y_pos_vec_spot - y_spot;
    r_sep_vec = sqrt(x_sep_vec.^2 + y_sep_vec.^2);
    sample_index_vec = 1:numel(x_sep_vec);
    % find closest pixel that meets criteria
    cr_filter = r_sep_vec >= min_sample_sep & round(dist_vec) == round(spot_dist);
    sample_distances = r_sep_vec(cr_filter);
    sample_index_vec = sample_index_vec(cr_filter);
    % if candidate found, then proceed. Else look to neighboring nuclei
    if ~isempty(sample_distances)
%         [~, mi] = min(sample_distances);
        sample_index = randsample(sample_index_vec,1);
        null_x = x_pos_vec_spot(sample_index);
        null_y = y_pos_vec_spot(sample_index);
        null_nc = index;
        qc_flag = 1;               
    else
        % Find nearest neighbor nucleus
        r_vec = r_dist_mat(:,index);
        r_vec(index) = Inf;
        [~, mi] = min(r_vec);
        % segment nn nucleus
        xn_nn = round(nc_x_vec(mi));
        yn_nn = round(nc_y_vec(mi));   
        x_spot_nn = round(spot_x_vec(mi));
        y_spot_nn = round(spot_y_vec(mi));   

        nc_bw_final_nn = segment_nc_neighborhood(histone_image, xn_nn, yn_nn, x_spot_nn, ...
        y_spot_nn, id_array, nb_sz, nc_index_vec(mi));

        null_mask = nc_bw_final_nn; % reassign null mask
        if ~isnan(x_spot_nn)
            nan_flag = isnan(nc_bw_final_nn(y_spot_nn,x_spot_nn));
        end
        % make sure size is reasonable 
        if sum(nc_bw_final_nn(:)) >= min_area && sum(nc_bw_final_nn(:)) <= max_area &&...
                (isnan(x_spot_nn) || ~nan_flag)

            % look for potential control samples                         
            if centroid_flag
                x_centroid = round(mean(x_ref(nc_bw_final_nn)));
                y_centroid = round(mean(y_ref(nc_bw_final_nn)));
                centroid_mat = zeros(size(x_ref));
                centroid_mat(y_centroid,x_centroid) = 1;
                dist_mat_nn = bwdist(centroid_mat);
            else
                dist_mat_nn = bwdist(~nc_bw_final_nn);
            end
            dist_vec = dist_mat_nn(nc_bw_final_nn);
            x_pos_vec = x_ref(nc_bw_final_nn);
            y_pos_vec = y_ref(nc_bw_final_nn);
            % distance from particle
            if ~isnan(x_spot_nn)                
                x_sep_vec = x_pos_vec - x_spot_nn;               
                y_sep_vec = y_pos_vec - y_spot_nn;
                r_sep_vec = sqrt(x_sep_vec.^2 + y_sep_vec.^2);
            else                
                r_sep_vec = Inf(size(dist_vec));
            end            
            sample_index_vec = 1:numel(r_sep_vec);
            % find closest pixel that meets criteria
            cr_filter = r_sep_vec >= min_sample_sep & round(dist_vec) == round(spot_dist);       
            sample_index_vec = sample_index_vec(cr_filter);
            if ~isempty(sample_index_vec)
                sample_index = randsample(sample_index_vec,1);
                null_x = x_pos_vec(sample_index);
                null_y = y_pos_vec(sample_index);
                null_nc = mi;
                qc_flag = 2;
            else
                null_x = NaN;
                null_y = NaN;
                null_nc = NaN;
                qc_flag = 0;
            end
        else
            null_x = NaN;
            null_y = NaN;
            null_nc = NaN;
            qc_flag = 0;
        end 
    end