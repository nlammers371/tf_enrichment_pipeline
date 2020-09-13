function nc_bw_final = segment_nc_neighborhood(histone_image, x_nucleus, y_nucleus, x_spot, ...
            y_spot, id_array, neighborhood_size, nc_ind)
        
    xDim = size(histone_image,2);
    yDim = size(histone_image,1);
    snip = histone_image(max(1,y_nucleus-neighborhood_size):min(yDim,y_nucleus+neighborhood_size),max(1,x_nucleus-neighborhood_size):min(xDim,x_nucleus+neighborhood_size));
    % generate binary histone image
    thresh = multithresh(snip(:));
    his_bin = histone_image > thresh;                
    his_lb = bwlabel(his_bin);
    % generate mask from neighborhood matrix
    id_mask = id_array==nc_ind;
    % make mask using binary histone image
    ID = his_lb(y_nucleus,x_nucleus);             
    nc_bw = his_lb==ID&(ID>0);          
    % enforce neighborhood boundaries
    nc_bw = nc_bw & id_mask;
    % enforce condition that spots must be inside nucleus
    if ~isnan(x_spot)
        nc_bw(y_spot,x_spot) = true;    
    end
    nc_bw_hull = bwconvhull(nc_bw);
    % prevent overlaps between nuclei
    nc_bw_final = nc_bw_hull & id_mask;
   