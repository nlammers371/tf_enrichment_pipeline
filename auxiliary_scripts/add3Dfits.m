function add3Dfits(project,DropboxFolder,varargin)
addpath('./utilities')
% set defaults
[RawResultsRoot, ~, ~] =   header_function(DropboxFolder, project);
LivemRNAPath = 'E:\Nick\LivemRNA\mRNADynamics';

for i = 1:numel(varargin)
    if ischar(varargin{i}) && i~= numel(varargin)        
        eval([varargin{i} '=varargin{i+1};']);        
    end
end

addpath(LivemRNAPath);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Set Path Specs, ID Vars %%%%%%%%%%%%%%%%%%%%%%%%

% find sheet
sheet_path = [RawResultsRoot 'DataStatus.xlsx'];
[~,sheet_names]=xlsfinfo(sheet_path);
sheet_index = find(ismember(sheet_names,project));
if isempty(sheet_index)
    error('no tab matching "project" string found in DataStatus')
end
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

disp('performing 3D fits to data...')
% Loop through filenames    
for i = 1:length(prefix_cell) 
    Prefix = prefix_cell{i};
    spot_path = [RawResultsRoot Prefix '/Spots.mat'];
    % check that Spots file exists
    if exist(spot_path)
        disp(['loading spots mat for  ' Prefix '...']);
        tic
        load(spot_path);
        toc
    else
        warning(['Spots structure not found for ' Prefix '. Skipping...']);
        continue
    end
    % call fitting function
    disp(['conducting fits for  ' Prefix '...']);
    tic
    fit3DGaussiansToAllSpots(Prefix, 1, 'segmentSpots',Spots);    
    toc
    disp(['finished fits for  ' Prefix ' (' num2str(i) ' of ' num2str(numel(prefix_cell))' '...']);    
end
