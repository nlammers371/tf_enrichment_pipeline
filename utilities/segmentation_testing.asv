project = 'Bcd_GFP_hb_mCherry_Zoom3x_LowPower';
RawPath = 'E:\Nick\LivemRNA\Data\PreProcessedData\';
% Load trace data
DataPath = ['../../dat/' project '/'];
load([DataPath '/nucleus_struct.mat'],'nucleus_struct')
load([DataPath '/set_key.mat'],'set_key')
WritePath = ['../../dat/' project '/mf_images/'];
mkdir(WritePath)

% generate indexing vectors
frame_ref = [nucleus_struct.frames];
nc_x_ref = [nucleus_struct.xPos];
nc_y_ref = [nucleus_struct.yPos];
spot_x_ref = [nucleus_struct.xPosParticle];
spot_y_ref = [nucleus_struct.yPosParticle];

set_ref = [];
ind_ref = [];
sub_ind_ref = [];
pt_ref = [];
for i = 1:numel(nucleus_struct)
    ParticleID = nucleus_struct(i).ParticleID;    
    set_ref = [set_ref repelem(nucleus_struct(i).setID,numel(nucleus_struct(i).frames))];
    ind_ref = [ind_ref repelem(i,numel(nucleus_struct(i).frames))];
    sub_ind_ref = [sub_ind_ref 1:numel(nucleus_struct(i).frames)];
    if isnan(ParticleID)
        pt_ref = [pt_ref false(1,numel(nucleus_struct(i).frames))];
    else
        pt_ref = [pt_ref ~isnan(nucleus_struct(i).xPosParticle)];
    end
end

qc_ref = false(size(set_ref)); % if 1, indicates successful segmentation
nn_ref = NaN(size(set_ref));
meta_ind_ref = 1:numel(nc_x_ref);
%%% set snip size to use
set_vec = [nucleus_struct.setID];
set_index = unique(set_vec);
px_sizes = [];
for s = 1:numel(set_index)
    px = [nucleus_struct(set_vec==set_index(s)).PixelSize];
    px_sizes = [px_sizes px(1)];
end
% determine size of neighborhood to use
nb_sizes = round(6 ./ px_sizes);
% size of gaussian smoothing kernel 
sm_kernels = round(.2 / px_sizes);
% set min and max acceptable area
min_areas = round(pi*(2 ./ px_sizes).^2);
max_areas = round(pi*(4 ./ px_sizes).^2);

%%% make source key
src_cell = {1,numel(set_index)};
pt_id_vec = [nucleus_struct.ParticleID];
for s = 1:numel(set_index)
    ind = find([nucleus_struct.setID]==set_index(s)&~isnan(pt_id_vec),1);   
    src = nucleus_struct(ind).source_path;    
    src_cell{s} = src;     
end
%%% Generate reference array for set-frame combos
set_frame_array = unique([set_ref' frame_ref'],'row');
ref_frame_struct = struct;
%%% iterate

parfor i = 1:size(set_frame_array,1)
    setID = set_frame_array(i,1);
    frame = set_frame_array(i,2);       
    % get nucleus positions
    nc_x_vec = nc_x_ref(set_ref==setID&frame_ref==frame);
    nc_y_vec = nc_y_ref(set_ref==setID&frame_ref==frame);        
    spot_x_vec = spot_x_ref(set_ref==setID&frame_ref==frame);
    spot_y_vec = spot_y_ref(set_ref==setID&frame_ref==frame);        
    
    ind_vec = ind_ref(set_ref==setID&frame_ref==frame);    
    pt_vec = pt_ref(set_ref==setID&frame_ref==frame);    
    % get size params
    nb_sz = nb_sizes(set_index==setID);
    sm_kernel = sm_kernels(set_index==setID);
    min_area = min_areas(set_index==setID);
    max_area = max_areas(set_index==setID);
    max_r = round(sqrt(max_area/pi))';
    % load and invert MCP mCherry    
    src = set_key(set_key.setID==setID,:).prefix{1};    
    mcp_files = dir([RawPath src '/*' sprintf('%03d',frame) '*_ch02.tif']);
    mcp_stack = [];
    for im = 2:numel(mcp_files)-1
        image = double(mat2gray(imread([RawPath src '/' mcp_files(im).name])));
        mcp_stack(:,:,im-1) = image;
    end
    % Invert image. Exclude spot loci
    mcp_med = mat2gray(nanmedian(mcp_stack,3));
    mcp_med_inv = 1-mcp_med;
    
    spot_x_regions = repmat([spot_x_vec-1 spot_x_vec spot_x_vec+1],1,3);
    spot_y_regions = [repmat(spot_y_vec-1,1,3) repmat(spot_y_vec,1,3) repmat(spot_y_vec+1,1,3)];
    indices = sub2ind(size(mcp_med),spot_y_regions,spot_x_regions);
    indices = indices(~isnan(indices));
    mcp_med_inv(indices) = prctile(mcp_med_inv(:),99);
    
    % smooth and normalize
    mcp_sm = imgaussfilt(mcp_med_inv,sm_kernel);
    his_sm = mcp_sm / mean(mcp_sm(:));    
   
    % dist ref arrays
    [x_cp, y_cp] = meshgrid(1:xDim, 1:yDim);
    % assign neighborhoods to each nucleus using a (reverse?) greedy allocation
    % process
    min_dist_array = Inf(size(x_cp));
    id_array = NaN(size(x_cp));
    for j = 1:numel(nc_x_vec)
        xn = nc_x_vec(j);
        yn = nc_y_vec(j);
        r_mat = sqrt((x_cp-xn).^2 + (y_cp-yn).^2);
        ids = find(r_mat<=max_r&r_mat<=min_dist_array);
        min_dist_array(ids) = r_mat(ids);
        id_array(ids) = ind_vec(j);
    end    
    % for each nucleus center, segment neighborhood define boundaries
    % nuclei near edges will be excluded 
    new_frame = NaN(size(his_slice));
    se = strel('disk',1);
    for j = 1:numel(nc_x_vec)
        xn = round(nc_x_vec(j));
        yn = round(nc_y_vec(j));   
        xp = round(spot_x_vec(j));
        yp = round(spot_y_vec(j));   
        snip = his_sm(max(1,yn-nb_sz):min(yDim,yn+nb_sz),max(1,xn-nb_sz):min(xDim,xn+nb_sz));
        % generate binary histone image
        thresh = multithresh(snip(:));
        his_bin = his_sm > thresh;                
        his_lb = bwlabel(his_bin);
        % generate mask from neighborhood matrix
        id_mask = id_array==ind_vec(j);
        id_mask = imerode(id_mask,se)~=0;
        % make mask using binary histone image
        ID = his_lb(yn,xn);             
        nc_bw = his_lb==ID&(ID>0);          
        % enforce neighborhood boundaries
        nc_bw = nc_bw & id_mask;
        % enforce condition that spots must be inside nucleus
        if ~isnan(xp)
            nc_bw(yp,xp) = true;
        end
        nc_bw_hull = bwconvhull(nc_bw);
        % prevent overlaps between nuclei
        nc_bw_final = nc_bw_hull & id_mask;
        % make sure size is reasonable
        if sum(nc_bw_final(:)) > min_area && sum(nc_bw_final(:)) < max_area            
            new_frame(nc_bw_final) = ind_vec(j);
        end                
    end      
    nn_vec = NaN(size(nc_x_vec));
    % assign nearest neighor ID for each spot in the frame
    nz_ids = find(~isnan(new_frame));
    [y_ids, x_ids] = find(~isnan(new_frame));
    id_list = new_frame(nz_ids);        
    for j = 1:numel(nc_x_vec)
        if pt_vec(j)
            % center of current nucleus
            xn = nc_x_vec(j);
            yn = nc_y_vec(j);
            % nn candidates
            xc = x_ids;
            yc = y_ids;
            xc(id_list==ind_vec(j)) = Inf;
            yc(id_list==ind_vec(j)) = Inf;
            % find closest mask
            rn = sqrt((xc-xn).^2+(yc-yn).^2);
            [min_r, mi] = min(rn);
            if id_list(mi) ~= ind_vec(j) && min_r <= 3*max_r             
                nn_vec(j) = id_list(mi);
            end
        end
    end
    % record frame info
    ref_frame_struct(i).ref_frame = new_frame;
    ref_frame_struct(i).set = setID;
    ref_frame_struct(i).frame = frame;        
    ref_frame_struct(i).nn_vec = nn_vec;        
    ref_frame_struct(i).ind_vec = ind_vec;            
end
