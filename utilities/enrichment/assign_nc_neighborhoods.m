function [id_array, yDim, xDim, y_ref, x_ref] = assign_nc_neighborhoods(histone_image, nc_x_vec, nc_y_vec, max_r, nc_index_vec)

    % dist ref arrays
    tic
    xDim = size(histone_image,2);
    yDim = size(histone_image,1);
    [x_ref, y_ref] = meshgrid(1:xDim, 1:yDim);
    % initialize binary array with nc centers
    dist_array = zeros(size(x_ref));
    nc_ind_list = sub2ind(size(x_ref),nc_y_vec,nc_x_vec);
    dist_array(nc_ind_list) = 1;
    % apply watershed
    dist_array = bwdist(dist_array);
    ws_array = watershed(dist_array);
    ws_array(dist_array>max_r) = NaN;    
    % reassign IDs
    id_array = NaN(size(ws_array));
    for j = 1:numel(nc_x_vec)
        init_id = ws_array(nc_y_vec(j),nc_x_vec(j));
        id_array(ws_array==init_id) = nc_index_vec(j);
    end
