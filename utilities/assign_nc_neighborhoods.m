function [id_array, yDim, xDim, y_ref, x_ref] = assign_nc_neighborhoods(histone_image, nc_x_vec, nc_y_vec, max_r, nc_index_vec)

    % dist ref arrays
    xDim = size(histone_image,2);
    yDim = size(histone_image,1);
    [x_ref, y_ref] = meshgrid(1:xDim, 1:yDim);
    % assign neighborhoods to each nucleus using a (reverse?) greedy allocation
    % process
    min_dist_array = Inf(size(x_ref));
    id_array = NaN(size(x_ref));
    for j = 1:numel(nc_x_vec)
        xn = nc_x_vec(j);
        yn = nc_y_vec(j);
        r_mat = sqrt((x_ref-xn).^2 + (y_ref-yn).^2);
        ids = find(r_mat<=max_r&r_mat<=min_dist_array);
        min_dist_array(ids) = r_mat(ids);
        id_array(ids) = nc_index_vec(j);
    end    