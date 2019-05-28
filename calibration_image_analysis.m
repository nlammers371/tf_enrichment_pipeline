% script to analyze 60-mer GFP data from injected embryos
clear
close all
% set file paths 
rawDynamicsData = 'E:\LocalEnrichment\Data\RawDynamicsData\';
date = '2019-05-21';
project = '60mer-eGFP_injectedEmbryos';
dataPath = [rawDynamicsData filesep date filesep project filesep 'AnalyzedImages' filesep];

% get list of TIFF files
file_list = dir([dataPath '*.tif']);

% determine list of unique z-stacks
trunc_name_cell = cell(1,numel(file_list));
for i = 1:numel(file_list)
    name = file_list(i).name;
    stopIndex = regexp(name,'t\d')+2;
    trunc_name_cell{i} = name(1:stopIndex);
end
stack_list = unique(trunc_name_cell);

% compile stacks
image_struct = struct;
for i = 1:numel(stack_list)
    stack_name = stack_list{i};
    indices = find(contains(trunc_name_cell,stack_name));
    im_stack = [];
    for j = indices
        im = imread([dataPath file_list(j).name]);
        im_stack = cat(3,im_stack,im);
    end
    image_struct(i).im_stack = im_stack;
    image_struct(i).source = stack_list{i};
end
% specify ranges to ignore during subsequent segmentation
% first four are easy
mask1 = true(size(im_stack,1),size(im_stack,2));
mask1(1:40,1:150) = false;
for i = 1:4
    image_struct(i).mask = mask1;
end
% now for first  set with visible embryo edge
im_project = mean(imgaussfilt(image_struct(5).im_stack,1),3);
mask = imbinarize(im_project);
mask = bwareaopen(mask,5);
mask_dist = bwdist(mask);
mask2 = true(size(mask1));
mask2(mask_dist<=50) = false;
for i = 5:8
    image_struct(i).mask = mask2;
end

% similar for second set with edge
im_project = mean(imgaussfilt(image_struct(9).im_stack,1),3);
mask = imbinarize(im_project);
mask = bwareaopen(mask,5);
mask_dist = bwdist(mask);
mask2 = true(size(mask1));
mask2(mask_dist<=50) = false;
for i = 9:12
    image_struct(i).mask = mask2;
end

%%% Segment and fit Z-projected stacks
for im = 1:numel(image_struct)
    close all
    show_figures = 0;
    disp('segmenting image...')
    % mean project
    mask = image_struct(im).mask;
    im_stack = image_struct(im).im_stack;
    for i = 1:size(im_stack,3)
        slice = im_stack(:,:,i);
        slice(~mask) = 0;
        im_stack(:,:,i) = slice;
    end
    im_gauss = imgaussfilt3(im_stack,1);
    im_project = max(im_gauss,[],3);
    if show_figures
        figure(1);
        imagesc(im_project)
    end
    % binarize using otsu's method
    im_project(~mask) = 0;
    im_bin = imbinarize(im_project);
    if show_figures
        figure(2);
        imagesc(im_bin)
    end
    % remove anamolously large regions as well as single pixels
    im_lb = bwlabel(im_bin);
    im_stats = regionprops(im_bin,'Area','Circularity','Centroid');
    area_vec = [im_stats.Area];
    keep_indices = find(area_vec > 2 & area_vec < 30);
    % keep_im_indices = find(ismember(im_lb(:),keep_indices));
    im_bin_size = zeros(size(im_bin));
    im_bin_size(ismember(im_lb,keep_indices)) = 1;
    im_stats_clean = im_stats(keep_indices);
    if show_figures
        figure(3);
        imagesc(im_bin_size)
    end
    % plot inferred spot centroids over top of original image
    centroid_array = vertcat(im_stats_clean.Centroid);
    if show_figures
        figure(4);
        imagesc(im_project);
        hold on
        scatter(centroid_array(:,1),centroid_array(:,2),20,'green')
    end
    % pull snips of local intensities for each putative spot for further
    % analysis
    xDim = size(im_bin,2);
    yDim = size(im_bin,1);
    snip_rad = 7;
    snip_stack = NaN(2*snip_rad+1,2*snip_rad+1,size(centroid_array,1));
    range_full = 1:2*snip_rad + 1;
    intensity_vec = NaN(1,size(centroid_array,1));
    for i = 20%1:size(centroid_array,1)
        xc = round(centroid_array(i,1));
        yc = round(centroid_array(i,2));
        x_range = xc-snip_rad:xc+snip_rad;
        y_range = yc-snip_rad:yc+snip_rad;
        ftx = x_range > 0 & x_range <=xDim;
        fty = y_range > 0 & y_range <=yDim;    
        snip_stack(fty,ftx,i) = im_mean_raw(y_range(fty),x_range(ftx));
    end

    % fit 2D Gaussians with offset to estimate total intensity of each putative
    % 60-mer
    disp('performing fits...')
    [mesh_x,mesh_y] = meshgrid(1:2*snip_rad+1,1:2*snip_rad+1);
    %fitting options
    lsqOptions=optimset('Display','none','maxfunevals',10000,'maxiter',10000);

    Gauss_fit_array_cst = NaN(size(snip_stack,3),5);
    Gauss_fit_array_full = NaN(size(snip_stack,3),7);
    Gauss_integral_vec_cst = NaN(1,size(snip_stack,3));
    Gauss_integral_vec_full = NaN(1,size(snip_stack,3));
    raw_integral_vec_cst = NaN(1,size(snip_stack,3));
    raw_integral_vec_full = NaN(1,size(snip_stack,3));

    % set parameter bounds
    lb = [.01 snip_rad-3 snip_rad-3 .1 .1 -1 0];
    ub = [10 snip_rad+3 snip_rad+3 5 5 1 1];
    for i = 1:size(snip_stack,3)
        snip2D = snip_stack(:,:,i);
        if sum(isnan(snip2D(:))) > .2 * numel(snip2D)
            continue
        end
        snip2D(isnan(snip2D)) = nanmean(snip2D(:));
        %%%%%%%%%%%%%%%%%
        % define constrained fit function
        Gauss2DConstrained = @(params) params(1).*...
        exp(-.5*(params(4)^-1.*(mesh_x-params(2)).^2  + params(4)^-1.*(mesh_y-params(3)).^2));    
        Gauss2DLossConstrained = @(params) Gauss2DConstrained(params) + params(5) - double(snip2D);    
        % perform fit
        initial_parameters = [1, snip_rad+1, snip_rad+1,1,.1];    
        [fits, ~, residual, ~, ~, ~, jacobian] = lsqnonlin(Gauss2DLossConstrained, ...
            initial_parameters,lb([1:4 7]),ub([1:4 7]), lsqOptions);    
        % record
        Gauss_fit_array_cst(i,:) = fits;
        sig_xy = sqrt(fits(4));
        Gauss_integral_vec_cst(i) = fits(1)*(2*pi*sig_xy^2);
        % raw intensity
        xy_rad = round(2*sig_xy);
        r_dist = sqrt((mesh_x-fits(2))^2+(mesh_y-fits(3))^2);
        raw_integral_vec_cst(i) = nansum(snip2D(r_dist<=xy_rad)) - fits(end)*sum(r_dist(:)<=xy_rad);

        %%%%%%%%%%%%%%
        % define "full" fit function
        Gauss2DFull = @(params) params(1).*...
        exp(-.5/(params(4)*params(5)*(1-params(6)^2)) * (params(5)*(mesh_x-params(2)).^2 + ...
            params(4)*(mesh_y-params(3)).^2 - 2*params(6)*sqrt(params(4)*params(5))*(mesh_x-params(2)).*(mesh_y-params(3))));    

        Gauss2DLossFull = @(params) Gauss2DFull(params) + params(7) - double(snip2D);
        % perfom fit
        initial_parameters = [1, snip_rad+1, snip_rad+1,1,1,.1,.1];        

        [fits, ~, residual, ~, ~, ~, jacobian] = lsqnonlin(Gauss2DLossFull, ...
            initial_parameters,lb,ub, lsqOptions); 
        % record
        Gauss_fit_array_full(i,:) = fits;
        sig = [fits(4) fits(6)*sqrt(fits(4)*fits(5)); fits(6)*sqrt(fits(4)*fits(5)) fits(5)];
        Gauss_integral_vec_full(i) = det(sig)^.5 * 2 * pi * fits(1);        
    end
    % save
    image_struct(im).Gauss_fit_array_full = Gauss_fit_array_full;
    image_struct(im).Gauss_fit_array_cst = Gauss_fit_array_cst;
    image_struct(im).Gauss_integral_vec_full = Gauss_integral_vec_full;
    image_struct(im).Gauss_integral_vec_cst = Gauss_integral_vec_cst;
    image_struct(im).raw_integral_vec_cst = raw_integral_vec_cst;
    image_struct(im).snip_stack = snip_stack;
    image_struct(im).im_stats = im_stats;
    image_struct(im).im_bin_size = im_bin_size;
end
%% Analyze results
cst_fit_array = vertcat(image_struct.Gauss_fit_array_cst);
xy_sig_vec = cst_fit_array(:,4);
cst_int_vec = [image_struct.Gauss_integral_vec_cst];
full_int_vec = [image_struct.Gauss_integral_vec_full];
histogram(cst_int_vec(xy_sig_vec<=1.5))
