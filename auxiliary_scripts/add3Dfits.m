function add3Dfits(project,DropboxFolder,varargin)
addpath('./utilities')
% set defaults
[RawResultsRoot, ~, ~] =   header_function(DropboxFolder, project);
LivemRNAPath = 'E:\Nick\LivemRNA\mRNADynamics';
minDT = 737735;
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
    cp_path = [RawResultsRoot Prefix '/CompiledParticles.mat'];
    sp_token_path = [RawResultsRoot Prefix '/Spots3DToken.mat'];
    cp_token_path = [RawResultsRoot Prefix '/CompiledParticlesToken.mat.mat'];
    % check that Spots file exists
    if exist(spot_path)
        % first check to see if 3D fits have been performed recently
        valid_fits = false;
        if exist(sp_token_path)
            load(sp_token_path)
            valid_fits = Spots3DToken > minDT;
        end   
        if valid_fits
            disp(['found recent fits for ' Prefix ' skipping...']);
        else
            disp(['loading spots mat for  ' Prefix '...']);
            tic
            load(spot_path);
            toc
            % call fitting function
            disp(['conducting fits for  ' Prefix '...']);
            tic
            fit3DGaussiansToAllSpots(Prefix, 1, 'segmentSpots',Spots);    
            toc
            disp(['finished fits for  ' Prefix ' (' num2str(i) ' of ' num2str(numel(prefix_cell))' '...']);    
        end              
        % now re-run CompileParticles 
        if exist(cp_path)
            % check for cp token
            valid_cp_token = false;
            if exist(cp_token_path)
                load(cp_token_path)
                valid_cp_token = (CompiledParticlesToken > Spots3DToken) & valid_fits;
            end
            if valid_cp_token
                disp('recently compiled particles set found. Skipping...')
            else
                disp('re-running CompileParticles...')
                tic
                CompileParticles(Prefix,'ApproveAll','SkipAll')
                toc
            end
        else
            disp('no CompiledParticles set found. Please run necessary pipeline scripts')
        end
    else
        warning(['Spots structure not found for ' Prefix '. Skipping...']);        
    end    
end
