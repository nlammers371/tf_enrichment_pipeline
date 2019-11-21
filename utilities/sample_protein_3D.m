function pt_avg = sample_protein_3D(x_spot,y_spot,z_spot,x_ref,y_ref,z_ref,xy_sigma,z_sigma,protein_stack)
    
    % generate rounded pos indices
    x_round = round(x_spot);
    y_round = round(y_spot);
    z_round = round(z_spot);
    % calculate vol dimensions
    xy_vol_dim = ceil(2*xy_sigma);
    z_vol_dim = ceil(2*z_sigma);
    % calculate dimensions
    xDim = size(x_ref,2);
    yDim = size(x_ref,1);
    zDim = size(x_ref,3);
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
        ismember(z_range3_full,z_range3)) = protein_stack(y_range3,x_range3,z_range3);
    % xyz ref
    y_ref_box = NaN(size(pt_samp_box));
    y_ref_box (ismember(y_range3_full,y_range3),ismember(x_range3_full,x_range3),...
        ismember(z_range3_full,z_range3)) = y_ref(y_range3,x_range3,z_range3);
    x_ref_box = NaN(size(pt_samp_box));
    x_ref_box(ismember(y_range3_full,y_range3),ismember(x_range3_full,x_range3),...
        ismember(z_range3_full,z_range3)) = x_ref(y_range3,x_range3,z_range3);
    z_ref_box = NaN(size(pt_samp_box));
    z_ref_box(ismember(y_range3_full,y_range3),ismember(x_range3_full,x_range3),...
        ismember(z_range3_full,z_range3)) = z_ref(y_range3,x_range3,z_range3);
    % generate weight box
    wt_box = exp(-.5*(((y_spot-y_ref_box)./xy_sigma).^2+((x_spot-x_ref_box)./xy_sigma).^2 + ...
        ((z_spot-z_ref_box)./z_sigma).^2));
    % take weighted average
    pt_avg = nansum(pt_samp_box(:).*wt_box(:)) / nansum(~isnan(pt_samp_box(:)).*wt_box(:));