clear
close all

% set paths
addpath('../utilities')
DropboxFolder =  'E:\Nick\LivemRNA\Dropbox (Personal)\';

% specify project
 project = 'Zld-GFP_snaBAC-mCh';
 
% get path to processed data
[ResultsRoot, ~, ~] =   header_function(DropboxFolder, project);
DataPath =  [DropboxFolder 'ProcessedEnrichmentData\' project '\'];
mkdir(DataPath);
sheet_path = [ResultsRoot 'DataStatus.xlsx'];
% path to raw tif stacks
RawDataPath = 'E:\LocalEnrichment\Data\PreProcessedData\';
ProteinChannel = 1; % should make this dynamic

% get prefix names
[~,sheet_names]=xlsfinfo(sheet_path);
sheet_index = find(ismember(sheet_names,project));
[~,~,sheet_cell] = xlsread(sheet_path,sheet_index);
name_col = sheet_cell(1:33,1); % hard coded for now
ready_ft = contains(name_col,'ReadyForEnrichment');
ready_cols = 1 + find([sheet_cell{ready_ft,2:end}]==1);
sheet_cell = sheet_cell(:,[1 ready_cols]);
% get list of project names
prefix_ft = contains(name_col,'Prefix');
prefix_cell_raw = sheet_cell(prefix_ft,2:end);
prefix_cell = {};
for i = 1:numel(prefix_cell_raw)
    if ~isempty(prefix_cell_raw{i})
        eval([prefix_cell_raw{i} ';'])
        prefix_cell = [prefix_cell{:} {Prefix}];
    end
end

%%% generate data structure with info about location of nuclei
tic
disp('initializing nuclear protein structure...')
nuclear_protein = [];
iter = 1;
for i = 1:numel(prefix_cell)
    load([ResultsRoot prefix_cell{i} '\' prefix_cell{i} '_lin.mat'])  
    load([ResultsRoot prefix_cell{i} '\FrameInfo.mat'])  
    % set nucleus rad
    PixelSize = FrameInfo(1).PixelSize;
    ncrad = ceil(4/PixelSize);
    % generate seperate entry for each tracked nucleus
    for j = 1:numel(schnitzcells)
        % initialize temp structure
        temp = struct;
        % find non-nan time points
        nan_flag = ~isnan(nanmean(schnitzcells(j).Fluo'));
        % get info
        if any(nan_flag)
            temp.xPos = schnitzcells(j).cenx(nan_flag);
            temp.yPos = schnitzcells(j).ceny(nan_flag);
            temp.ncrad = ncrad;
            [~, temp.zPos] = max(schnitzcells(j).Fluo(nan_flag,:)'); % let's sample brightest z slice

            temp.frames = schnitzcells(j).frames(nan_flag)';
            temp.APpos = schnitzcells(j).APpos(nan_flag);
            temp.Prefix = prefix_cell{i};
            temp.PrefixIDVec = repelem(i,numel(temp.yPos));
            temp.iter_id = repelem(iter,numel(temp.yPos));
            temp.iter_sub_id = 1:numel(temp.yPos);
            % initialize structure to store protein samples
            temp.nc_protein_stack = NaN(2*temp.ncrad+1,2*temp.ncrad+1,numel(temp.frames));
            % add to master structure
            nuclear_protein = [nuclear_protein temp];
            % increment
            iter = iter + 1;
        end
    end    
end
toc

% generate indexing vectors 
set_ref = [nuclear_protein.PrefixIDVec];
frame_ref = [nuclear_protein.frames];
x_ref = [nuclear_protein.xPos];
y_ref = [nuclear_protein.yPos];
z_ref = [nuclear_protein.zPos];
nc_id_ref = [nuclear_protein.iter_id];
nc_sub_id_ref = [nuclear_protein.iter_sub_id];
% get array of unique set-frame combinations
set_frame_array = unique([set_ref' frame_ref'],'row');
disp('extracting nuclear protein stacks...')
% initialize waitbar 
f = waitbar(0,'extracting nuclear protein');
  
for i = 1:size(set_frame_array,1)    
    % generate filter
    setID = set_frame_array(i,1);
    frame = set_frame_array(i,2);  
    frame_set_filter = set_ref == setID & frame_ref == frame;
    % extract position vectors
    nc_x_vec = x_ref(frame_set_filter);
    nc_y_vec = y_ref(frame_set_filter);
    nc_z_vec = z_ref(frame_set_filter);
    % indexing vectors
    nc_id_vec = nc_id_ref(frame_set_filter);
    nc_sub_id_vec = nc_sub_id_ref(frame_set_filter);    
    % get basic frame info 
    Prefix = prefix_cell{setID};
    load([ResultsRoot Prefix '\FrameInfo.mat'])
    yDim = FrameInfo(1).LinesPerFrame;
    xDim = FrameInfo(1).PixelsPerLine;
    zDim = FrameInfo(1).NumberSlices;
    
    % load stacks        
    files = dir([RawDataPath Prefix '/*_' sprintf('%03d',frame) '*_ch0' num2str(ProteinChannel) '.tif']);
    protein_stack = NaN(yDim,xDim,zDim);
    for im = 2:zDim + 1      
        protein_stack(:,:,im-1) = double(imread([RawDataPath Prefix '/' files(im).name]));
    end    
    
    % iterate through nuclei
    for j = 1:numel(nc_x_vec)
        % position
        ncx = nc_x_vec(j);
        ncy = nc_y_vec(j);
        ncz = nc_z_vec(j)-1;
        % index
        nci = nc_id_vec(j);
        ncsi = nc_sub_id_vec(j);
        ncrad = nuclear_protein(nci).ncrad;
        % indexing vectors
        x_range_full = ncx - ncrad:ncx + ncrad;
        x_ft = x_range_full > 0 & x_range_full <= xDim;
        y_range_full = ncy - ncrad:ncy + ncrad;
        y_ft = y_range_full > 0 & y_range_full <= yDim;
        % store protein
        nuclear_protein(nci).nc_protein_stack(y_ft,x_ft,ncsi) = protein_stack(y_range_full(y_ft),x_range_full(x_ft),ncz);
    end
    waitbar(i/size(set_frame_array,1),f)
end
clear f
disp('saving...')
save([DataPath 'nuclear_protein.mat'],'nuclear_protein')
disp('done')