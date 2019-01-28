function [id_array, yDim, xDim, y_ref, x_ref] = assign_nc_neighborhoods(histone_image, nc_x_vec, nc_y_vec, max_r, nc_index_vec)

    % dist ref arrays
    xDim = size(histone_image,2);
    yDim = size(histone_image,1);
    [x_ref, y_ref] = meshgrid(1:xDim, 1:yDim);
    % assign neighborhoods to each nucleus using a (reverse?) greedy allocation
    % process
    min_dist_array = Inf(size(x_ref));
    id_array = NaN(size(x_ref));
    temp_mat = Inf(size(x_ref));
    for j = 1:numel(nc_x_vec)
        xn = round(nc_x_vec(j));
        yn = round(nc_y_vec(j));
        if xn < 2 || yn < 2 || xDim - xn < 2 || yDim - yn < 2
            continue
        end
        x_range = max(1,xn-round(max_r)):min(xDim,xn+round(max_r));
        y_range = max(1,yn-round(max_r)):min(yDim,yn+round(max_r));
       
        r_sub = sqrt((x_ref(y_range,x_range)-xn).^2 + (y_ref(y_range,x_range)-yn).^2);
      
        r_full = temp_mat;
        r_full(y_range,x_range) = r_sub;
        ids = find(r_full<=max_r&r_full<=min_dist_array);
        min_dist_array(ids) = r_full(ids);
        id_array(ids) = nc_index_vec(j);
    end    