% main03_make_ref_images(project, zeiss_data,)
%
% DESCRIPTION
% Sample protein levels ant active loci and control regions
%
% ARGUMENTS
% project: master ID variable 
%
% RawPath: Full or relative path to PreProcessed folder (pr
%               equivalent)

%
% OUTPUT: ref_frame_struct: compiled data set

function [snip_struct, nucleus_struct_pt] = main03_sample_protein_levels(project,RawPath,zeiss_flag,pt_channel)

%%% Load trace data
DataPath = ['../../dat/' project '/'];

nucleus_struct = load([DataPath  'nucleus_struct.mat']);
nucleus_struct = nucleus_struct.nucleus_struct;
ref_frame_struct = load([DataPath  'ref_frame_struct.mat']);
ref_frame_struct = ref_frame_struct.ref_frame_struct;
%%% get path to ref images
RefPath = [DataPath '/mf_images/'];
nc_set_vec = [nucleus_struct.setID];
set_index = unique(nc_set_vec);
%%% set parameters for protein analysis 
pt_window = 1; % dim of region to integrate for pt level at locus
%%% set snip size to use
px_sizes = [];
for s = 1:numel(set_index)
    px = [nucleus_struct(nc_set_vec==set_index(s)).PixelSize];
    px_sizes = [px_sizes px(1)];
end
% set snippet to be 3um in size
pt_snippet_size = round(1.5 / mean(px_sizes));
% set min separation between samples
min_sample_sep = round(.5 / mean(px_sizes));

% vectors to index ref frame
ref_frame_sets = [ref_frame_struct.set]; 
ref_frame_frames = [ref_frame_struct.frame];

src_cell = {1,numel(set_index)};
pt_id_vec = [nucleus_struct.ParticleID];
for s = 1:numel(set_index)
    ind = find([nucleus_struct.setID]==set_index(s)&~isnan(pt_id_vec),1);   
    src = nucleus_struct(ind).source_path;    
    src_cell{s} = src;     
end

%%% initialize fields in nucleus structure
snip_names = {'fluo_snippet_spot','fluo_snippet_null','pt_snippet_spot',...
    'pt_snippet_null'};
vec_names = {'fluo_spot','fluo_null','pt_spot','pt_null','xPosNull_edge',...
    'yPosNull_edge','edge_flags','ctrl_flags','edgeDistNull','edgeDistSpot','CenterDistSpot'};
% add new fields to nucleus structure
for i = 1:numel(nucleus_struct)
    nc_time = nucleus_struct(i).time;    
    for j = 1:numel(vec_names)
        nucleus_struct(i).(vec_names{j}) = NaN(1,numel(nc_time));
    end    
end

%%% get list of x and y and t coordinates
pt_id_vec = [nucleus_struct.ParticleID];

% position and time vectors
ptX_vec = [nucleus_struct.xPosParticle];
ptY_vec = [nucleus_struct.yPosParticle];
ncX_vec = [nucleus_struct.xPos];
ncY_vec = [nucleus_struct.yPos];
ptZ_vec = [nucleus_struct.brightestZs];
ptFrame_vec = [nucleus_struct.frames];
ptTime_vec = [nucleus_struct.time];
AP_vec = [nucleus_struct.ap_vector];
ptID_vec = [];
for i = 1:numel(nucleus_struct)
    ptID_vec = [ptID_vec repelem(nucleus_struct(i).ParticleID, numel(nucleus_struct(i).brightestZs))];
end
% filter for time steps where spot is present
ptID_vec = ptID_vec(~isnan(ptX_vec));
set_vec = floor(ptID_vec);
ncX_vec = ncX_vec(~isnan(ptX_vec));
ncY_vec = ncY_vec(~isnan(ptX_vec));
ptY_vec = ptY_vec(~isnan(ptX_vec));
ptZ_vec = ptZ_vec(~isnan(ptX_vec));
ptFrame_vec = ptFrame_vec(~isnan(ptX_vec));
ptTime_vec = ptTime_vec(~isnan(ptX_vec));
AP_vec = AP_vec(~isnan(ptX_vec));
ptX_vec = ptX_vec(~isnan(ptX_vec));

%%% initialize temp arrays to store parfor results
fluo_snippet_spot_array = NaN(2*pt_snippet_size+1,2*pt_snippet_size+1,numel(AP_vec));
fluo_snippet_null_array = NaN(2*pt_snippet_size+1,2*pt_snippet_size+1,numel(AP_vec));

pt_snippet_spot_array = NaN(2*pt_snippet_size+1,2*pt_snippet_size+1,numel(AP_vec));
pt_snippet_null_array = NaN(2*pt_snippet_size+1,2*pt_snippet_size+1,numel(AP_vec));

%%% initialize vectors to store central intensity info
fluo_spot_vec = NaN(1,numel(AP_vec));
fluo_null_vec = NaN(1,numel(AP_vec));
pt_spot_vec = NaN(1,numel(AP_vec));
pt_null_vec = NaN(1,numel(AP_vec));

%%% aaaand vectors to store position info
xPosNull_edge_vec = NaN(1,numel(AP_vec));
yPosNull_edge_vec = NaN(1,numel(AP_vec));

edgeDistNull_vec = NaN(1,numel(AP_vec));
edgeDistSpot_vec = NaN(1,numel(AP_vec));
CenterDistSpot_vec = NaN(1,numel(AP_vec));
%%% iterate through projects
tic
ctrl_flags_vec = NaN(size(ptX_vec)); % if 1, full control
edge_flags_vec = NaN(size(ptX_vec)); % if 1, snippet truncated due to edge

id_null = 0;
id_ct = 0;
r_ct = 0;
parfor i = 1:numel(ptX_vec)
    % basic id info
    ParticleID = ptID_vec(i);  
    nc_ind = find(pt_id_vec==ParticleID);
    setID = floor(ParticleID);
    src = src_cell{set_index==floor(ParticleID)};    
    % get ref frame and space/time info
    frame = ptFrame_vec(i);     
    % load
    ref_name = [RefPath '/ref_frame_' sprintf('%03d',frame) '_set' sprintf('%03d',setID) '.mat'];
    rf = load(ref_name,'rf');
    
    ref_frame = rf.rf;       
    % all active particle positions (must avoid sampling near these)
    xp_vec = ptX_vec(set_vec==setID&ptFrame_vec==frame)-zeiss_flag;
    yp_vec = ptY_vec(set_vec==setID&ptFrame_vec==frame)-zeiss_flag;
    % particle of interest
    xp = ptX_vec(i)-zeiss_flag;
    yp = ptY_vec(i)-zeiss_flag;                
    
    % nucleus center positon    
    xn = ncX_vec(i);
    yn = ncY_vec(i);    
    % particle z position
    zp = ptZ_vec(i);  
    % distance from nucleus center
    CenterDistSpot_vec(i) = sqrt((xn-xp).^2+(yn-yp).^2);
    
    %%%%%%%%%% Perform Control Selection Using Nucleus Masks %%%%%%%%%%%%%%
    % initilaize control flag to 1 
    ctrl_edge_flag = 1;
    % get particle ID
    ID = ref_frame(yp,xp);    
    % get NN id   
    nn_vec = ref_frame_struct(ref_frame_sets==setID&ref_frame_frames==frame).nn_vec;
    ind_vec = ref_frame_struct(ref_frame_sets==setID&ref_frame_frames==frame).ind_vec;
    nn_ID = nn_vec(ind_vec==nc_ind);    
    nn_mask = ref_frame==nn_ID;
    D_nn = bwdist(~nn_mask);    
    % First try to select control snip from neighbor nucleus
    % Check that spot falls within ID'd region in ref frame, and that ID
    % corresponds to recorded value
    if ~isnan(ID) && ID == ind_vec(ind_vec==nc_ind)
        ref_mask = ref_frame==ID;
        % get distance of spot to nuclear edge
        D = bwdist(~ref_mask);
        spot_dist = D(yp,xp);
        edgeDistSpot_vec(i) = spot_dist;
        
        ids = find(round(D_nn)==round(spot_dist));        
        [y_ctrl, x_ctrl] = ind2sub(size(nn_mask),ids);
        r_ctrl = nanmin(sqrt((y_ctrl-yp_vec).^2+(x_ctrl-xp_vec).^2),[],2);
        s_indices = 1:numel(x_ctrl); 
        s_indices = s_indices(r_ctrl>=min_sample_sep);
        % if at least one pixel in NN nucleus fits criteria, take it
        if ~isempty(s_indices)
            control_index_edge = randsample(s_indices,1);       
            x_ctrl_edge = x_ctrl(control_index_edge);
            y_ctrl_edge = y_ctrl(control_index_edge);        
            edgeDistNull_vec(i) = D_nn(control_index_edge);
        else            
            r_ct = r_ct + 1;
            ctrl_edge_flag = 0;
        end
    else          
        id_ct = id_ct + 1;
        id_null = id_null + isnan(ID);
        ctrl_edge_flag = 0;
        control_index_edge = NaN;
    end       
    % load particle frame        
    mcp_name = [RawPath src '/' src '_' sprintf('%03d',frame) '_z' sprintf('%02d',zp) '_ch0' num2str(1*(2~=pt_channel) + 1) '.tif'];
    mcp_frame = imread(mcp_name);      
    mcp_frame_spot = mcp_frame;
    if ctrl_edge_flag
        mcp_frame_spot = NaN(size(ref_mask));
        mcp_frame_spot(ref_mask) = double(mcp_frame(ref_mask));
        mcp_frame_ctrl = NaN(size(nn_mask));
        mcp_frame_ctrl(nn_mask) = double(mcp_frame(nn_mask));        
    end
    % load protein frame 
    pt_name = [RawPath src '/' src '_' sprintf('%03d',frame) '_z' sprintf('%02d',zp) '_ch0' num2str(pt_channel) '.tif'];
    pt_frame = imread(pt_name);   
    pt_frame_spot = pt_frame;
    if ctrl_edge_flag
        pt_frame_spot = NaN(size(ref_mask));
        pt_frame_spot(ref_mask) = double(pt_frame(ref_mask));
        pt_frame_ctrl = NaN(size(ref_mask));
        pt_frame_ctrl(nn_mask) = double(pt_frame(nn_mask));             
    end
    
    full_size = pt_snippet_size*2 + 1;            
    %%%%%%%%%%%% Active Locus Intensities %%%%%%%%%%%%
    % spot MCP
    fluo_snippet_spot = mcp_frame_spot(max(1,yp-pt_snippet_size):min(size(mcp_frame,1),yp...            
        +pt_snippet_size),max(1,xp-pt_snippet_size):min(size(mcp_frame,2),xp+pt_snippet_size));      
    % calculate central intensities    
    fs = mcp_frame_spot(max(1,yp-pt_window):min(size(mcp_frame,1),yp...            
        +pt_window),max(1,xp-pt_window):min(size(mcp_frame,2),xp+pt_window));  
    fluo_spot_vec(i) = nanmean(fs(:));
    
    % spot Protein
    pt_snippet_spot = pt_frame_spot(max(1,yp-pt_snippet_size):min(size(mcp_frame,1),yp...            
        +pt_snippet_size),max(1,xp-pt_snippet_size):min(size(mcp_frame,2),xp+pt_snippet_size));    
    % protein central intensity        
    ps = pt_frame_spot(max(1,yp-pt_window):min(size(mcp_frame,1),yp...            
        +pt_window),max(1,xp-pt_window):min(size(mcp_frame,2),xp+pt_window));          
    pt_spot_vec(i) = nanmean(ps(:));        
    
    % record snippets
    spot_edge_flag = 0;
    if size(pt_snippet_spot,1) ==  full_size && size(pt_snippet_spot,2) == full_size
        fluo_snippet_spot_array(:,:,i) = fluo_snippet_spot;
        pt_snippet_spot_array(:,:,i) = pt_snippet_spot;
        spot_edge_flag = 1;
    end
    %%%%%% Control-- nucleus mask method %%%%%%%%
    null_edge_flag = 0;
    if ctrl_edge_flag
        % sample coordinates
        xc = x_ctrl_edge;
        yc = y_ctrl_edge; 
        % record control coordinates
        xPosNull_edge_vec(i) = xc;
        yPosNull_edge_vec(i) = yc;    
        % control MCP
        fluo_snippet_null = mcp_frame_ctrl(max(1,yc-pt_snippet_size):min(size(mcp_frame,1),yc...            
            +pt_snippet_size),max(1,xc-pt_snippet_size):min(size(mcp_frame,2),xc+pt_snippet_size));   
        % central intensity
        fn = mcp_frame_ctrl(max(1,yc-pt_window):min(size(mcp_frame,1),yc...            
            +pt_window),max(1,xc-pt_window):min(size(mcp_frame,2),xc+pt_window));   
        fluo_null_vec(i) = nanmean(fn(:));
        % control Protein                
        pt_snippet_null = pt_frame_ctrl(max(1,yc-pt_snippet_size):min(size(mcp_frame,1),yc...            
            +pt_snippet_size),max(1,xc-pt_snippet_size):min(size(mcp_frame,2),xc+pt_snippet_size));
        % central intensity
        pn = pt_frame_ctrl(max(1,yc-pt_window):min(size(mcp_frame,1),yc...            
            +pt_window),max(1,xc-pt_window):min(size(mcp_frame,2),xc+pt_window));
        pt_null_vec(i) = nanmean(pn(:));        
        % record snippets        
        if size(pt_snippet_null,1) == full_size && size(pt_snippet_null,2) == full_size
%             fluo_snippet_null(isnan(fluo_snippet_null)) = -1;
            fluo_snippet_null_array(:,:,i) = fluo_snippet_null;
%             pt_snippet_null(isnan(pt_snippet_null)) = -1;
            pt_snippet_null_array(:,:,i) = pt_snippet_null;   
            null_edge_flag = 1;
        end
    end            
    %%% save vectors
    ctrl_flags_vec(i) = ctrl_edge_flag;        
    edge_flags_vec(i) = spot_edge_flag&null_edge_flag;    
end
toc
disp(r_ct)
disp(id_ct)
disp(id_null)
% add new fields to nucleus structure
for i = 1:numel(pt_null_vec)
    ind = find(ptID_vec(i)==pt_id_vec);
    frame = ptFrame_vec(i);
    all_frames = nucleus_struct(ind).frames;        
    for j = 1:numel(vec_names)
        vn = vec_names{j};        
        vec = eval([vn '_vec']);
        nucleus_struct(ind).(vn)(all_frames==frame) = vec(i);
    end        
end
% filter nucleus_struct
nucleus_struct_pt = nucleus_struct(~isnan(pt_id_vec));

for j = 1:numel(snip_names)
    arr = eval([snip_names{j} '_array']);
    snip_struct.(snip_names{j}) = arr;
end
% add fields to snippet struct

snip_struct.ptID_vec = ptID_vec;
snip_struct.set_vec = set_vec;
snip_struct.ncX_vec = ncX_vec;
snip_struct.ncY_vec = ncY_vec;
snip_struct.ptY_vec = ptY_vec;
snip_struct.ptZ_vec = ptZ_vec;
snip_struct.ptFrame_vec = ptFrame_vec;
snip_struct.AP_vec = AP_vec;
snip_struct.ptX_vec = ptX_vec;
snip_struct.ptTime_vec = ptTime_vec;
snip_struct.ctrl_flags_edge = ctrl_flags_vec;

%%% save
save('-v7.3',[DataPath '/nucleus_struct_pt.mat'],'nucleus_struct_pt')
save('-v7.3',[DataPath '/snip_struct.mat'],'snip_struct')