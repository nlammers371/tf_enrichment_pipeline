function data_snip = sample_snip(x_spot,y_spot,pt_snippet_size,data_frame,spot_nc_mask)
    % get frame dimensions
    [yDim, xDim] = size(data_frame);
    % define range variables
    x_range = max(1,x_spot-pt_snippet_size):min(xDim,x_spot+pt_snippet_size);
    y_range = max(1,y_spot-pt_snippet_size):min(yDim,y_spot+pt_snippet_size);
    x_range_full = x_spot-pt_snippet_size:x_spot+pt_snippet_size;
    y_range_full = y_spot-pt_snippet_size:y_spot+pt_snippet_size;
    % define snip
    data_snip = NaN(numel(y_range_full),numel(x_range_full));
    data_snip(ismember(y_range_full,y_range),ismember(x_range_full,x_range)) = data_frame(y_range,x_range);
    % set everything outside nucleus mask to be NAN
    bound_snip = false(numel(y_range_full),numel(x_range_full));
    bound_snip(ismember(y_range_full,y_range),ismember(x_range_full,x_range)) = spot_nc_mask(y_range,x_range);  
    
    data_snip(~bound_snip) = NaN;