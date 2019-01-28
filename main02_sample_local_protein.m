% main02_sample_local_protein(project, RawPath, keyword)
%
% DESCRIPTION
% Script to generate samples of local protein contentration at active gene
% loci and selected control locations
%
% ARGUMENTS
% project: master ID variable 
%
% RawPath: Full or relative path to PreProcessed folder (pr
%               equivalent)
%
%
% OUTPUT: ref_frame_struct: compiled data set

function nucleus_struct_protein = main02_sample_local_protein(project,RawPath,protein_channel,varargin)
tic
zeiss_flag = 0;
for i = 1:numel(varargin)
    if strcmpi(varargin{i}, 'zeiss')
        zeiss_flag = 1;
    end
end
% Load trace data
DataPath = ['../dat/' project '/'];
load([DataPath '/nucleus_struct.mat'],'nucleus_struct')
load([DataPath '/set_key.mat'],'set_key')
SnipPath = ['../dat/' project '/qc_images/'];
mkdir(SnipPath)
addpath('./utilities')
% get MCP channel
options = [1 2];
mcp_channel = options(options~=protein_channel);
% generate indexing vectors
frame_ref = [nucleus_struct.frames];
nc_x_ref = [nucleus_struct.xPos];
nc_y_ref = [nucleus_struct.yPos];
spot_x_ref = [nucleus_struct.xPosParticle]+zeiss_flag;
spot_y_ref = [nucleus_struct.yPosParticle]+zeiss_flag;
spot_z_ref = [nucleus_struct.brightestZs];

set_ref = [];
ind_ref = [];
sub_ind_ref = [];
pt_ref = [];
for i = 1:numel(nucleus_struct)
    ParticleID = nucleus_struct(i).ParticleID;    
    pt_ref = [pt_ref repelem(ParticleID, numel(nucleus_struct(i).frames))];
    set_ref = [set_ref repelem(nucleus_struct(i).setID,numel(nucleus_struct(i).frames))];
    ind_ref = [ind_ref repelem(i,numel(nucleus_struct(i).frames))];
    sub_ind_ref = [sub_ind_ref 1:numel(nucleus_struct(i).frames)];
end

%%% set snip size to use
set_vec = [nucleus_struct.setID];
set_index = unique(set_vec);
px_sizes = [];
for s = 1:numel(set_index)
    px = [nucleus_struct(set_vec==set_index(s)).PixelSize];
    px_sizes = [px_sizes px(1)];
end% determine size of neighborhood to use
nb_sizes = round(10 ./ px_sizes);
% size of gaussian smoothing kernel 
sm_kernels = round(.2 / px_sizes);
% set min and max acceptable area
min_areas = round(pi*(2 ./ px_sizes).^2);
max_areas = round(pi*(4 ./ px_sizes).^2);
% set snippet to be 3um in size
pt_snippet_size_vec = round(1.5 ./ px_sizes);
% set min separation between control and locus to 2um
min_sample_sep_vec = round(2 ./ px_sizes);

% Designate fields ot be added to nucleus structure
new_vec_fields = {'spot_protein_vec', 'null_protein_vec', 'spot_mcp_vec','null_mcp_vec',...
    'qc_flag_vec', 'fov_edge_flag_vec', 'null_x_vec', 'null_y_vec', 'null_nc_vec',...
    'spot_edge_dist_vec', 'spot_edge_dist_vec'};

new_snip_fields = {'spot_protein_snips', 'null_protein_snips','spot_mcp_snips','null_mcp_snips'};
% get most common snip size
default_snip_size = mode(pt_snippet_size_vec);
% Initialize fields
for i = 1:numel(nucleus_struct)
    ref = nucleus_struct(i).xPos;
    for j = 1:numel(new_vec_fields)
        nucleus_struct(i).(new_vec_fields{j}) = NaN(size(ref));
    end
    for j = 1:numel(new_snip_fields)
        nucleus_struct(i).(new_snip_fields{j}) = NaN(2*default_snip_size+1,2*default_snip_size+1,numel(ref));
    end
end
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
qc_structure = struct;
%%% iterate
for i = 1:size(set_frame_array,1)
    setID = set_frame_array(i,1);
    frame = set_frame_array(i,2);       
    % get nucleus and particle positions
    frame_set_filter = set_ref==setID&frame_ref==frame;
    nc_x_vec = nc_x_ref(frame_set_filter);
    nc_y_vec = nc_y_ref(frame_set_filter);        
    spot_x_vec = spot_x_ref(frame_set_filter);
    spot_y_vec = spot_y_ref(frame_set_filter);        
    spot_z_vec = spot_z_ref(frame_set_filter);  
    particle_id_vec = pt_ref(frame_set_filter);
    % find indices for spots
    if sum(~isnan(spot_x_vec)) == 0
        continue
    end
    % indexing variables
    nc_sub_index_vec = sub_ind_ref(frame_set_filter);     
    nc_index_vec = ind_ref(frame_set_filter);    
       
    %%%%% First load MCP frames and generate Nucleus Segmentation frame%%%%
    
    % get size params    
    pt_snippet_size = pt_snippet_size_vec(set_index==setID); % size of snippet
    min_sample_sep = min_sample_sep_vec(set_index==setID); % minimum distance between spot and control
    nb_sz = nb_sizes(set_index==setID); % size of nucleus neighborhood to use
    sm_kernel = sm_kernels(set_index==setID); %sigma for gaussian smoothing kerne;
    min_area = min_areas(set_index==setID); % lower bound on permitted nucleus size
    max_area = max_areas(set_index==setID); % upper bound
    max_r = round(sqrt(max_area/pi))'; % max nucleus neighborhood size
    % determine whether snips will need to be resampled
    snip_scale_factor =  default_snip_size / pt_snippet_size;
    % load and  MCP mCherry and protein stacks
    src = set_key(set_key.setID==setID,:).prefix{1};        
    [mcp_stack, protein_stack] = load_stacks(RawPath, src, frame, mcp_channel);
    
    % Invert image
    mcp_med = mat2gray(nanmedian(mcp_stack,3));
    mcp_med_inv = 1-mcp_med;
    % un-invert pixels around spot center
    spot_x_regions = repmat([spot_x_vec-1 spot_x_vec spot_x_vec+1],1,3);
    spot_y_regions = [repmat(spot_y_vec-1,1,3) repmat(spot_y_vec,1,3) repmat(spot_y_vec+1,1,3)];
    indices = sub2ind(size(mcp_med),spot_y_regions,spot_x_regions);
    indices = indices(~isnan(indices));
    mcp_med_inv(indices) = prctile(mcp_med_inv(:),99);
    
    % smooth and normalize
    mcp_sm = imgaussfilt(mcp_med_inv,sm_kernel);
    his_sm = mcp_sm / mean(mcp_sm(:));    
    
    [id_array, yDim, xDim, y_ref, x_ref] = assign_nc_neighborhoods(his_sm, nc_x_vec, nc_y_vec, max_r, nc_index_vec);
    % generate lookup table of inter-nucleus distances
    x_dist_mat = repmat(nc_x_vec,numel(nc_x_vec),1)-repmat(nc_x_vec',1,numel(nc_x_vec));
    y_dist_mat = repmat(nc_y_vec,numel(nc_y_vec),1)-repmat(nc_y_vec',1,numel(nc_y_vec));
    r_dist_mat = sqrt(x_dist_mat.^2 + y_dist_mat.^2);
    % for each spot, segment nearby nuclei and attempt to sample local
    % protein levels        
    % initialize arrays to store relevant info 
    for j = 1:numel(new_vec_fields)
        eval([new_vec_fields{j} ' = NaN(size(spot_x_vec));']);
    end
    for j = 1:numel(new_snip_fields)
        eval([new_snip_fields{j} ' = NaN(2*pt_snippet_size+1,2*pt_snippet_size+1,numel(spot_x_vec));']);
    end    
    
    % iterate through spots
    qc_mat = struct;
    parfor j = 1:numel(nc_x_vec)        
        % get location info
        xn = round(nc_x_vec(j));
        yn = round(nc_y_vec(j));   
        xp = round(spot_x_vec(j));
        yp = round(spot_y_vec(j));           
        zp = round(spot_z_vec(j));   
        if isnan(xp)
            continue
        end
        nc_bw_final = segment_nc_neighborhood(his_sm, xn, yn, xp, ...
            yp, id_array, nb_sz, nc_index_vec(j));
        % make sure size is reasonable and that spot is inside nucleus
        if sum(nc_bw_final(:)) < min_area || sum(nc_bw_final(:)) > max_area || ~nc_bw_final(yp,xp)   
            qc_flag_vec(j) = 0;
            continue
        end 
        
        %%%%%%%%%%%%%%%%%%%%%% Sample protein levels %%%%%%%%%%%%%%%%%%%%%%
        % get frames
        protein_frame = protein_stack(:,:,zp);
        mcp_frame = mcp_stack(:,:,zp);
        % take locus samples
        spot_protein_vec(j) = protein_frame(yp,xp);
        spot_mcp_vec(j) = mcp_frame(yp,xp);
        % pull snippets
        null_mask = nc_bw_final;
        temp_pt = protein_frame;
        temp_pt(~null_mask) = NaN;
        temp_mcp = mcp_frame;
        temp_mcp(~null_mask) = NaN;
        x_range = max(1,xp-pt_snippet_size):min(xDim,xp+pt_snippet_size);
        y_range = max(1,yp-pt_snippet_size):min(yDim,yp+pt_snippet_size);
        if numel(x_range) == 2*pt_snippet_size+1 && numel(y_range) == 2*pt_snippet_size+1
            spot_protein_snips(:,:,j) = temp_pt(y_range,x_range);
            spot_mcp_snips(:,:,j) = temp_mcp(y_range,x_range);
            fov_edge_flag_vec(j) = 0;
        else
            fov_edge_flag_vec(j) = 1;
        end
        % Now attempt to find control snip. First try same nucleus
        % distance from edge        
        dist_mat = bwdist(~nc_bw_final);        
        spot_edge_dist = dist_mat(yp,xp);
        spot_edge_dist_vec(j) = spot_edge_dist;
        edge_dist_vec = dist_mat(nc_bw_final);
        % distance from particle
        x_pos_vec = x_ref(nc_bw_final);
        x_sep_vec = x_pos_vec - xp;
        y_pos_vec = y_ref(nc_bw_final);
        y_sep_vec = y_pos_vec - yp;
        r_sep_vec = sqrt(x_sep_vec.^2 + y_sep_vec.^2);
        sample_index_vec = 1:numel(x_sep_vec);
        % find closest pixel that meets criteria
        cr_filter = r_sep_vec >= min_sample_sep & round(edge_dist_vec) == round(spot_edge_dist);
        sample_distances = r_sep_vec(cr_filter);
        sample_index_vec = sample_index_vec(cr_filter);
        % if candidate found, then proceed. Else look to neighboring nuclei
        if ~isempty(sample_distances)
            [~, mi] = min(sample_distances);
            null_x_vec(j) = x_pos_vec(sample_index_vec(mi));
            null_y_vec(j) = y_pos_vec(sample_index_vec(mi));
            null_nc_vec(j) = j;
            qc_flag_vec(j) = 1;        
        else
            % Find nearest neighbor nucleus
            r_vec = r_dist_mat(:,j);
            r_vec(j) = Inf;
            [~, mi] = min(r_vec);
            % segment nn nucleus
            xn_nn = round(nc_x_vec(mi));
            yn_nn = round(nc_y_vec(mi));   
            xp_nn = round(spot_x_vec(mi));
            yp_nn = round(spot_y_vec(mi));   
            
            nc_bw_final_nn = segment_nc_neighborhood(his_sm, xn_nn, yn_nn, xp_nn, ...
            yp_nn, id_array, nb_sz, nc_index_vec(mi));
        
            null_mask = nc_bw_final_nn; % reassign null mask
            if ~isnan(xp_nn)
                nan_flag = isnan(nc_bw_final_nn(yp_nn,xp_nn));
            end
            % make sure size is reasonable 
            if sum(nc_bw_final_nn(:)) >= min_area && sum(nc_bw_final_nn(:)) <= max_area &&...
                    (isnan(xp_nn) || ~nan_flag)
           
                % look for potential control samples                
                % distance from edge        
                dist_mat_nn = bwdist(~nc_bw_final_nn);
                edge_dist_vec = dist_mat_nn(nc_bw_final_nn);
                x_pos_vec = x_ref(nc_bw_final_nn);
                y_pos_vec = y_ref(nc_bw_final_nn);
                % distance from particle
                if ~isnan(xp_nn)                
                    x_sep_vec = x_pos_vec - xp_nn;               
                    y_sep_vec = y_pos_vec - yp_nn;
                    r_sep_vec = sqrt(x_sep_vec.^2 + y_sep_vec.^2);
                else                
                    r_sep_vec = Inf(size(edge_dist_vec));
                end            
                sample_index_vec = 1:numel(r_sep_vec);
                % find closest pixel that meets criteria
                cr_filter = r_sep_vec >= min_sample_sep & round(edge_dist_vec) == round(spot_edge_dist);       
                sample_index_vec = sample_index_vec(cr_filter);
                if ~isempty(sample_index_vec)
                    sample_index = randsample(sample_index_vec,1);
                    null_x_vec(j) = x_pos_vec(sample_index);
                    null_y_vec(j) = y_pos_vec(sample_index);
                    null_nc_vec(j) = mi;
                    qc_flag_vec(j) = 2;
                else
                    qc_flag_vec(j) = 0;
                end
            else
                qc_flag_vec(j) = 0;      
            end 
        end
        % phew. Now draw control samples (as appropriate)
        xc = NaN;
        yc = NaN;
        if qc_flag_vec(j) > 0
            xc = null_x_vec(j);
            yc = null_y_vec(j);            
            null_protein_vec(j) = protein_frame(yc,xc);
            null_mcp_vec(j) = mcp_frame(yc,xc);
            % draw snips
            temp_pt = protein_frame;
            temp_pt(~null_mask) = NaN;
            temp_mcp = mcp_frame;
            temp_mcp(~null_mask) = NaN;
            x_range = max(1,xc-pt_snippet_size):min(xDim,xc+pt_snippet_size);
            y_range = max(1,yc-pt_snippet_size):min(yDim,yc+pt_snippet_size);
            if numel(x_range) == 2*pt_snippet_size+1 && numel(y_range) == 2*pt_snippet_size+1
                null_protein_snips(:,:,j) = temp_pt(y_range,x_range);
                null_mcp_snips(:,:,j) = temp_mcp(y_range,x_range);
                fov_edge_flag_vec(j) = 0;
            else
                fov_edge_flag_vec(j) = 1;
            end
        end  
        % save qc data                 
        qc_mat(j).setID = setID;
        qc_mat(j).frame = frame;
        qc_mat(j).nc_index = nc_index_vec(j);
        qc_mat(j).nc_sub_index = nc_sub_index_vec(j);
        qc_mat(j).qc_flag = qc_flag_vec(j);            
        qc_mat(j).xp = xp;
        qc_mat(j).yp = yp;  
        qc_mat(j).xc = xc;
        qc_mat(j).yc = yc;
        qc_mat(j).ParticleID = particle_id_vec(j);
        sz = nb_sz;
%         rescale = 1;
        if qc_flag_vec(j) == 2
            xd = abs(xn - xc);
            yd = abs(yn - yc);
            sz = max([nb_sz,yd,xd]);
%             rescale = sz / nb_sz;
            dist_mat = dist_mat + dist_mat_nn;
        end
        y_range = max(1,yn-sz):min(yDim,yn+sz);
        x_range = max(1,xn-sz):min(xDim,xn+sz);
        qc_mat(j).x_center = median(x_range);
        qc_mat(j).y_center = median(y_range);
        qc_mat(j).mcp_snip = his_sm(y_range,x_range);
        qc_mat(j).dist_snip = dist_mat(y_range,x_range);        
    end 
    qc_structure(i).qc_mat = qc_mat;
    % map data back to nucleus_struct
    for j = 1:numel(nc_index_vec)
        nc_index = nc_index_vec(j);
        nc_sub_index = nc_sub_index_vec(j);
        for k = 1:numel(new_vec_fields)
            vec = eval(new_vec_fields{k});
            nucleus_struct(nc_index).(new_vec_fields{k})(nc_sub_index) = vec(j);
        end
        for k = 1:numel(new_snip_fields)
            snip = eval(new_snip_fields{k});
            if snip_scale_factor ~= 1
                snip = imresize(snip,snip_scale_factor);
            end
            nucleus_struct(nc_index).(new_snip_fields{k})(:,:,nc_sub_index) = snip(:,:,j);
        end
    end       
end

% save qc data
for i = 1:numel(qc_structure)
    qc_mat = qc_structure(i).qc_mat;
    for  j = 1:numel(qc_mat)
        qc_spot = qc_mat(j);
        if ~isfield(qc_spot,'ParticleID')
            continue
        end
        ParticleID = qc_spot.ParticleID;
        if isempty(ParticleID)
            continue
        end        
        frame = qc_spot.frame;      
        save_name = [SnipPath 'pt' num2str(1e4*ParticleID) '_frame' sprintf('%03d',frame) '.mat'];
        save(save_name,'qc_spot');
    end
end
toc
% save updated nucleus structure
nucleus_struct_protein = nucleus_struct;
save([DataPath 'nucleus_struct_protein.mat'],'nucleus_struct_protein') 