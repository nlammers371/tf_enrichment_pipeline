% script to analyze 60-mer GFP data from injected embryos
clear
close all
% set file paths 
rawDynamicsData = 'E:\LocalEnrichment\Data\RawDynamicsData\';
date = '2019-05-21';
project = '60mer-eGFP_injectedEmbryos';
load("E:\Nick\Dropbox (Garcia Lab)\ProcessedEnrichmentData\Dl-Ven_snaBAC-mCh\psf_dims.mat");
dataPath = [rawDynamicsData filesep date filesep project filesep 'AnalyzedImages' filesep];
% specify fit parameters 
xy_snip_rad = 7;
z_snip_rad = 3;
min_blob_size = 5; % minimum number of pixels in blob
max_blob_size = 50; % this is a bit arbitrary at the moment
pixelSize = .1079;
zSize = .2;
ROIRadiusSpot = .2;
psf_dims_px = [psf_dims.xy_sigma psf_dims.xy_sigma psf_dims.z_sigma];
% generate approximate psf
% psf_filter = nonIsotropicGaussianPSF(psf_dims_px,3);
% vol_px = pi^1.5 * prod(psf_dims_px);
% D = ceil(size(psf_filter,1)/2);
% psf_filter = psf_filter(D-xy_snip_rad:D+xy_snip_rad,D-xy_snip_rad:D+xy_snip_rad,D-z_snip_rad:D+z_snip_rad);
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
% apply masks
for im = 1:numel(image_struct)
    mask = image_struct(im).mask;
    im_stack = image_struct(im).im_stack;
    for i = 1:size(im_stack,3)
        slice = im_stack(:,:,i);
        slice(~mask) = 0;
        im_stack(:,:,i) = slice;
    end
    image_struct(im).im_stack = im_stack;
    image_struct(im).dummy_flag = false;    
end
% generate set of random stacks to test validity of segmentation
full_stack = cat(3,image_struct.im_stack);
for im = numel(image_struct)+1:2*numel(image_struct)
    fake_vec = randsample(full_stack(:),numel(im_stack),true);
    image_struct(im).im_stack = reshape(fake_vec,size(im_stack));
    image_struct(im).dummy_flag = true;
    image_struct(im).mask = true(size(im_stack,1),size(im_stack,2));
    image_struct(im).source = 'fake';
end
%%% Segment and fit Z stacks
for im = 1:numel(image_struct)
    close all
    show_figures = 0;
    disp(['segmenting image ' num2str(im) ' (' num2str(numel(image_struct)) ')...'])
    % first remove small rafts of "hot" pixels  
    im_stack = image_struct(im).im_stack;
    im_stack_bin = bwareaopen(im_stack>0,min_blob_size);    
    im_stack_clean = double(im_stack).*double(im_stack_bin);
    % perform difference of gaussians filtering  
    im_gauss1 = imgaussfilt3(im_stack_clean,1);
    im_gauss5 = imgaussfilt3(im_stack_clean,5);
    im_gaussd = im_gauss1-im_gauss5;
    im_project = max(im_gaussd,[],3);
    if show_figures
        figure(1);
        imagesc(im_project)
    end
    % binarize using otsu's method
%     thresh = multithresh(im_gaussd,1);
    im_stack_bin = im_gaussd>.25; % set by eye currently
    % estimate bkg
    bkg_fluo = nanmean(im_stack(~im_stack_bin));
    % get image dimensions
    xDim = size(im_stack_bin,2);
    yDim = size(im_stack_bin,1);
    zDim = size(im_stack_bin,3);  
    % make labeled stack
    im_stack_lb = bwlabeln(im_stack_bin);
    im_stats = regionprops3(im_stack_lb,'Centroid','Solidity','Volume','VoxelList');
    % find number of unique z slices covered by each potential spot
    n_z_vec = NaN(1,numel(im_stats.VoxelList));
    for k = 1:numel(im_stats.VoxelList)
        voxel_list = im_stats.VoxelList{k};
        n_z_vec(k) = numel(unique(voxel_list(:,3)));
    end
    z_vec = im_stats.Centroid(:,3);
    vol_vec = [im_stats.Volume];
    sol_vec = [im_stats.Solidity];
    keep_indices = find(vol_vec >=min_blob_size & vol_vec <= max_blob_size & z_vec >1 & z_vec < zDim & sol_vec >=.8&n_z_vec'>1);
    im_stats_clean = im_stats(keep_indices,:);
    % plot inferred spot centroids over top of original image
    centroid_array = vertcat(im_stats_clean.Centroid);
    if show_figures
        figure(4);
        imagesc(im_project);
        hold on
        scatter(centroid_array(:,1),centroid_array(:,2),20,'green')
    end
    n_spots = size(centroid_array,1);
    % generate vector indicating fake vs. (putative) true spots
    if image_struct(im).dummy_flag
        spot_status_vec = -1*ones(1,n_spots);
    else
        dummy_centroids = [randsample(5:yDim-5,n_spots,true)',randsample(5:xDim-5,n_spots,true)',randsample(2:zDim-1,n_spots,true)'];
        centroid_array = vertcat(centroid_array,dummy_centroids);
        spot_status_vec = [ones(1,n_spots) zeros(1,n_spots)];
    end    
    % Attempt to fit 2D gaussians to all putative spots      
    snip_stack = zeros(2*xy_snip_rad+1,2*xy_snip_rad+1,2*z_snip_rad+1,size(centroid_array,1));    
    intensity_vec = NaN(1,size(centroid_array,1));
    for i = 1:size(centroid_array,1)
        xc = round(centroid_array(i,1));
        yc = round(centroid_array(i,2));
        zc = round(centroid_array(i,3));
        x_range = xc-xy_snip_rad:xc+xy_snip_rad;
        y_range = yc-xy_snip_rad:yc+xy_snip_rad;
        z_range = zc-z_snip_rad:zc+z_snip_rad;
        ftx = x_range > 0 & x_range <=xDim;
        fty = y_range > 0 & y_range <=yDim; 
        ftz = z_range > 0 & z_range <=zDim; 
        snip_stack(fty,ftx,ftz,i) = im_stack(y_range(fty),x_range(ftx),z_range(ftz));
    end

    % fit 3D Gaussians with offset to estimate total intensity of each putative
    % 60-mer
    disp(['parforming fits ' num2str(im) ' (' num2str(numel(image_struct)) ')...'])
    [mesh_x,mesh_y,mesh_z] = meshgrid(1:2*xy_snip_rad+1,1:2*xy_snip_rad+1,1:2*z_snip_rad+1);
    %fitting options
    lsqOptions=optimset('Display','none','maxfunevals',10000,'maxiter',10000);    
   
    Gauss_fit_array = NaN(size(snip_stack,3),8);    
    Gauss_integral_vec = NaN(1,size(snip_stack,3));    
    raw_integral_vec = NaN(1,size(snip_stack,3));    
    roi_integral_vec = NaN(1,size(snip_stack,3));    
    % set parameter bounds
    lb = [.01 xy_snip_rad-3 xy_snip_rad-3 z_snip_rad-1 .2 .2 .5 0];
    ub = [10 xy_snip_rad+3 xy_snip_rad+3 z_snip_rad+1 1.5 1.5 3 10];
    for i = 1:size(snip_stack,4)
        snip3D = snip_stack(:,:,:,i);   
        snipFilt = imgaussfilt3(snip3D,1);
        [~,mi] = max(snipFilt(:));
        [yc,xc,zc] = ind2sub(size(snip3D),mi);
        dx_mesh = abs(mesh_x-xc)*pixelSize;
        dy_mesh = abs(mesh_y-yc)*pixelSize;
        dz_mesh = abs(mesh_z-zc);
%         ft_mat = dx_mesh<=psf_dims_px(2)&dy_mesh<=psf_dims_px(1)&dz_mesh<=psf_dims_px(3);
        ft_mat = sqrt(dx_mesh(:,:,zc).^2 + dy_mesh(:,:,zc).^2)<=ROIRadiusSpot;
        roi_integral_vec(i) = sum(sum(ft_mat.*snip3D(:,:,zc))) - sum(ft_mat(:))*bkg_fluo;
        %%%%%%%%%%%%%%%%%
        % define constrained fit function
%         single3DGaussian = @(params) params(1)*...
%             exp(-.5*( ((params(2) - mesh_x)/params(5)).^2 + ((params(3)-mesh_y)/params(6)).^2 ...
%                  + ((params(4)-mesh_z)/params(7)).^2 ));
% 
%         objective_fun = @(params) single3DGaussian(params) + params(8) - snip3D;  
%         % perform fit
%         initial_parameters = [1, xy_snip_rad+1, xy_snip_rad+1 z_snip_rad+1,1,1,1,.1];    
%         [fits, ~, residual, ~, ~, ~, jacobian] = lsqnonlin(objective_fun, ...
%             initial_parameters,lb,ub, lsqOptions);    
        % record
%         Gauss_fit_array(i,:) = fits;       
%         if sqrt(fits(5)/fits(6))>.5  && sqrt(fits(5)/fits(6))<2 
%             Gauss_integral_vec(i) = fits(1)*(2*pi)^1.5 / sqrt(fits(4)*fits(5)*fits(6));
%             % raw intensity
%             gauss_pd = single3DGaussian(fits) ./ Gauss_integral_vec(i);
%             raw_integral_vec(i) = nansum(snip3D(gauss_pd>.01)) - fits(end)*sum(gauss_pd(:)>.01);   
%             % ROI intensity
%             z_plane = round(fits(4));
%             r_mat = sqrt(((fits(2)-mesh_x)*pixelSize).^2+((fits(3)-mesh_y)*pixelSize).^2+((fits(4)-mesh_z)*zSize).^2);
%             z_slice = snip3D(:,:,z_plane);
%             roi_integral_vec(i) = sum(snip3D(r_mat<=ROIRadiusSpot)) - fits(end)*sum(r_mat(:)<=ROIRadiusSpot);
%         end
%         Gauss_integral_vec(i) = sum(sum(sum(psf_filter.*snip3D)));
    end
    % save    
    image_struct(im).Gauss_fit_array = Gauss_fit_array;    
    image_struct(im).Gauss_integral_vec = Gauss_integral_vec;
    image_struct(im).raw_integral_vec = raw_integral_vec;
    image_struct(im).roi_integral_vec = roi_integral_vec;
    image_struct(im).snip_stack = snip_stack;
    image_struct(im).im_stats = im_stats_clean;    
    image_struct(im).spot_status_vec = spot_status_vec;
end

%%% Analyze results
Gauss_fit_array = vertcat(image_struct.Gauss_fit_array);
gauss_int_vec = [image_struct.Gauss_integral_vec];
roi_int_vec = [image_struct.roi_integral_vec];
raw_int_vec = [image_struct.raw_integral_vec];
spot_flag_vec = [image_struct.spot_status_vec];

qc_roi_fig = figure;
hold on
histogram(roi_int_vec(spot_flag_vec==1),'Normalization','probability')
xlabel('fluorescence (au)')
ylabel('share')
grid on
xlim([-5 30])

au_per_gfp = nanmedian(roi_int_vec(spot_flag_vec==1)) * 6 / 50 / 60

