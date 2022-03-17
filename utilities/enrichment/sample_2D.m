function pt_avg = sample_2D(x_spot,y_spot,x_ref,y_ref,xy_sigma,image_slice)
    
    % generate rounded pos indices
    x_round = round(x_spot);
    y_round = round(y_spot);
    
    % calculate vol dimensions
    xy_vol_dim = ceil(2*xy_sigma);
    
    % calculate dimensions
    xDim = size(x_ref,2);
    yDim = size(x_ref,1);
    
    % volume protein sampling 
    x_range3 = max(1,x_round-xy_vol_dim):min(xDim,x_round+xy_vol_dim);
    y_range3 = max(1,y_round-xy_vol_dim):min(yDim,y_round+xy_vol_dim);
    
    x_range3_full = x_round-xy_vol_dim:x_round+xy_vol_dim;
    y_range3_full = y_round-xy_vol_dim:y_round+xy_vol_dim;    

    % generate protein sample and ref boxes
    pt_samp_box = NaN(numel(y_range3_full),numel(x_range3_full));
    pt_samp_box(ismember(y_range3_full,y_range3),ismember(x_range3_full,x_range3)...
        ) = image_slice(y_range3,x_range3);
      
    % xyz ref
    y_ref_box = NaN(size(pt_samp_box));
    y_ref_box (ismember(y_range3_full,y_range3),ismember(x_range3_full,x_range3),...
        ismember(z_range3_full,z_range3)) = y_ref(y_range3,x_range3,z_range3);
    x_ref_box = NaN(size(pt_samp_box));
    x_ref_box(ismember(y_range3_full,y_range3),ismember(x_range3_full,x_range3),...
        ismember(z_range3_full,z_range3)) = x_ref(y_range3,x_range3,z_range3);

    % generate weight box
    wt_box = exp(-.5*(((y_spot-y_ref_box)./xy_sigma).^2+((x_spot-x_ref_box)./xy_sigma).^2 ...
        ));
    % take weighted average
    pt_avg = nansum(pt_samp_box(:).*wt_box(:)) / nansum(~isnan(pt_samp_box(:)).*wt_box(:));