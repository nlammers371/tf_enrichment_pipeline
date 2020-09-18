function snip2D = sample_snip_3D(x_spot,y_spot,z_spot,z_ref,snip_size,z_sigma,data_stack)
    
    % generate rounded pos indices
    x_round = round(x_spot);
    y_round = round(y_spot);
    z_round = round(z_spot);
    
    % calculate vol dimensions
    xy_vol_dim = snip_size;
    z_vol_dim = ceil(2*z_sigma);
    
    % calculate dimensions
    xDim = size(z_ref,2);
    yDim = size(z_ref,1);
    zDim = size(z_ref,3);
   
    % volume protein sampling 
    x_range3 = max(1,x_round-xy_vol_dim):min(xDim,x_round+xy_vol_dim);
    y_range3 = max(1,y_round-xy_vol_dim):min(yDim,y_round+xy_vol_dim);
    z_range3 = max(1,z_round-z_vol_dim):min(zDim,z_round+z_vol_dim);
    x_range3_full = x_round-xy_vol_dim:x_round+xy_vol_dim;
    y_range3_full = y_round-xy_vol_dim:y_round+xy_vol_dim;
    z_range3_full = z_round-z_vol_dim:z_round+z_vol_dim;

    % generate protein sample and ref boxes
    pt_samp_box = NaN(numel(y_range3_full),numel(x_range3_full),numel(z_range3_full));
    pt_samp_box(ismember(y_range3_full,y_range3),ismember(x_range3_full,x_range3),...
        ismember(z_range3_full,z_range3)) = data_stack(y_range3,x_range3,z_range3);
    % z ref    
    z_ref_box(ismember(y_range3_full,y_range3),ismember(x_range3_full,x_range3),...
        ismember(z_range3_full,z_range3)) = z_ref(y_range3,x_range3,z_range3);
      
    % generate weight box
    wt_box = exp(-.5*(((z_spot-z_ref_box)./z_sigma).^2));
      
    % take weighted average
    snip2D = nansum(pt_samp_box.*wt_box,3) ./ nansum(~isnan(pt_samp_box).*wt_box,3);
    