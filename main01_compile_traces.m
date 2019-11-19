% main01_compile_data(project,DropboxFolder,varargin)
%
% DESCRIPTION
% Funcion to compile relevant outputs from image analysis pipeline across
% multiple experiments

%
% ARGUMENTS
% project: master ID variable (should match a tab name in the Data Status
% sheet)
% DropboxFolder: full file path to folder containing compiled imaging
% results

% OPTIONS
% dropboxFolder: Pass this option, followed by the path to data folder 
%                where you wish to save
%                pipeline-generated data sets and figures. If this
%                var is not specified, output will be saved one level
%                above git repo in the folder structure
% first_nc: script defaults to taking only nc14 traces. If you want
%           earlier traces, pass 'first_nc', followed by desired nuclear cycle
%           number
%
% OUTPUT: nucleus_struct: compiled data set contain key nucleus and
% particle attributes

function nucleus_struct = main01_compile_traces(project,DropboxFolder,varargin)
addpath('./utilities')
% set defaults
firstNC = 14;
minDP = 15;
pctSparsity = 1;
two_spot_flag = contains(project, '2spot');
min_time = 6*60; % take no fluorescence data prior to this point
TresInterp = 20; 
[RawResultsRoot, DataPath, ~] =   header_function(DropboxFolder, project);
for i = 1:numel(varargin)
    if ischar(varargin{i})
        if ismember(varargin{i},{'includeVec','firstNC','expType','minDP'})
           eval([varargin{i} '=varargin{i+1};']);
        end
    end
end

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
    
% make filepath
mkdir(DataPath);
% assign save names
nucleus_name = [DataPath 'nucleus_struct.mat']; % names for compiled elipse struct

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Obtain Relevant Filepaths %%%%%%%%%%%%%%%%%%%%%%%
% initialize filename vectors
cp_filenames = {}; % particles
cn_filenames = {}; % protein data
nc_filenames = {}; % nuclei
fov_filenames = {}; % fov info
for d = 1:numel(prefix_cell)
    thisdir = prefix_cell{d};            
    % append file paths
    cn_filenames = [cn_filenames {[RawResultsRoot thisdir '/CompiledNuclei.mat']}];
    cp_filenames = [cp_filenames {[RawResultsRoot thisdir '/CompiledParticles.mat']}];        
    nc_filenames = [nc_filenames {[RawResultsRoot thisdir '/' thisdir '_lin.mat']}];           
    fov_filenames = [fov_filenames {[RawResultsRoot thisdir '/FrameInfo.mat']}];        
end

% generate set key data structure
set_key = array2table((1:numel(prefix_cell))','VariableNames',{'setID'});
set_key.prefix = prefix_cell';
save([DataPath 'set_key.mat'],'set_key')

% Generate master structure with info on all nuclei and traces in
% constituent sets

disp('compiling data...')
nucleus_struct = [];
% Loop through filenames    
for i = 1:length(cp_filenames) 
    % read in raw files
    try
        load(nc_filenames{i}) % Ellipse Info
        load(fov_filenames{i}) % FrameInfo Info                
    catch
        warning(['failed to load one or more files for prefix: ' prefix_cell{i} '. Skipping...' ])
        continue
    end     
    
    % extract data structures
    processed_data = load(cp_filenames{i}); % processed particles    
    % Extract compiled particles structure         
    cp_particles = processed_data.CompiledParticles;    
    if iscell(cp_particles)
        cp_particles = cp_particles{1};
    end   
%     threeD_flag = isfield(cp_particles,'Fluo3DRaw');
    threeD_flag = isfield(cp_particles,'FluoGauss3D');
    % check for 3D fit data
    % set identifier
    setID = i;            
    % pull trace, time, and frame variables
    time_raw = processed_data.ElapsedTime*60; % time vector   
    traces_raw = processed_data.AllTracesVector; % array with a column for each trace 
    % check to see if traces are stored in cell array
    if iscell(traces_raw)
        traces_raw = traces_raw{1};
    end
    frames_raw = 1:length(time_raw); % Frame list   
    % find index for first frame
    first_frame = max([1, processed_data.(['nc' num2str(firstNC)])]); % Default to start of nc14 for the moment    
    % filter trace mat and time
    traces_clean = traces_raw(first_frame:end,:);
    time_clean = time_raw(first_frame:end);    
    time_clean = time_clean - min(time_clean); % Normalize to start of nc    
    frames_clean = frames_raw(first_frame:end);        
    % get basic frame info
    yDim = FrameInfo(1).LinesPerFrame;
    xDim = FrameInfo(1).PixelsPerLine;
    zDim = FrameInfo(1).NumberSlices;
    zStep = FrameInfo(1).ZStep;
    % compile schnitz info
    s_cells = struct;
    e_pass = 1;    
    for e = 1:length(schnitzcells)
        e_frames = schnitzcells(e).frames;
        nc_filter = ismember(e_frames,frames_clean);
        nc_frames = e_frames(nc_filter);        
        if length(nc_frames) >= 1 % skip nuclei not desired nc range   
            s_cells(e_pass).xDim = xDim;
            s_cells(e_pass).yDim = yDim;
            s_cells(e_pass).zDim = zDim;
            s_cells(e_pass).zStep = zStep;
            s_cells(e_pass).PixelSize = FrameInfo(1).PixelSize;    
            
            % Initialize particle fields...will be set to particle values 
            % for nuclei with matching particle
            s_cells(e_pass).ParticleID = NaN;
            s_cells(e_pass).xPosParticle = NaN(1,sum(nc_filter));
            s_cells(e_pass).yPosParticle = NaN(1,sum(nc_filter));
            s_cells(e_pass).zPosParticle = NaN(1,sum(nc_filter)); 
            if threeD_flag
                s_cells(e_pass).xPosParticle3D = NaN(1,sum(nc_filter));
                s_cells(e_pass).yPosParticle3D = NaN(1,sum(nc_filter));
                s_cells(e_pass).zPosParticle3D = NaN(1,sum(nc_filter)); 
                s_cells(e_pass).fluo3D = NaN(1,sum(nc_filter));
            end
            s_cells(e_pass).fluo = NaN(1,sum(nc_filter)); 
            s_cells(e_pass).fluoOffset = NaN(1,sum(nc_filter));  

            s_cells(e_pass).qc_flag = NaN;
            s_cells(e_pass).N = NaN;
            s_cells(e_pass).sparsity = NaN;
            
            % add core nucleus info
            x = schnitzcells(e).cenx;            
            y = schnitzcells(e).ceny;                           
            s_cells(e_pass).xPos = x(nc_filter);
            s_cells(e_pass).yPos = y(nc_filter); 
            s_cells(e_pass).frames = nc_frames';            
            s_cells(e_pass).Nucleus = e;     
            s_cells(e_pass).ncID = eval([num2str(setID) '.' sprintf('%04d',e)]);
            s_cells(e_pass).xMean = mean(x(nc_filter));
            s_cells(e_pass).yMean = mean(y(nc_filter));
            s_cells(e_pass).ncStart = firstNC;
            s_cells(e_pass).minDP = minDP;
            s_cells(e_pass).min_time = min_time;
            s_cells(e_pass).pctSparsity = pctSparsity;
            
            % time and set info
            s_cells(e_pass).time = time_clean(ismember(frames_clean,nc_frames));
            if numel(s_cells(e_pass).time) > 1
                if s_cells(e_pass).time(2) == 0
                    error('wtf')
                end
            end
            s_cells(e_pass).setID = setID;
            fn = cp_filenames{i}; % Get filename to store in struct  
            fn = fn(1:strfind(fn,'/')-1);
            s_cells(e_pass).source_path = fn;                                     
            s_cells(e_pass).ap_flag = 0;
            
            % sister fields
            s_cells(e_pass).sister_Index = NaN;
            s_cells(e_pass).sister_ParticleID = NaN; 
            e_pass = e_pass + 1;
        end
    end
    % Index vector to cross-ref w/ particles            
    e_index = [s_cells.Nucleus];          
   
    % iterate through traces
    for j = 1:size(traces_clean,2)  
        % Raw fluo trace
        raw_trace = traces_clean(:,j); 
        % Get nucleus ID
        schnitz = cp_particles(j).schnitz;        
        trace_start = find(time_clean>=min_time&~isnan(raw_trace'),1);
        trace_stop = find(~isnan(raw_trace),1,'last');        
        %Creat versions with all intervening frames present (missing frames
        %appear as NaNs)
        trace_full = raw_trace(trace_start:trace_stop)';                                   
        frames_full = frames_clean(trace_start:trace_stop); 
        % skip particles not in nc range
        if sum(~isnan(trace_full)) == 0
            continue
        end                              
        % Find intersection btw full frame range and CP frames        
        raw_pt_frames = cp_particles(j).Frame;
        % perform qc tests    
        nDP = sum(~isnan(trace_full));
        sparsity = prctile(diff(find(~isnan(trace_full))),pctSparsity);
        qc_flag = ~(nDP < minDP || sparsity > 1);     
        % trace-nucleus mapping may be many-to-1
        nc_ind = find(e_index==schnitz);  
        if length(nc_ind) ~= 1
            warning('Problem with Particle-Nucleus Crossref')
            continue
        end 
        % Identifier variable                        
        particle = cp_particles(j).OriginalParticle;                        
        ParticleID = eval([num2str(setID) '.' sprintf('%04d',particle)]);
        
        % check to see if a different particle has already been assigned
        % if so, create a new entry
        if ~isnan(s_cells(nc_ind).ParticleID)
            two_spot_flag = true;
            s_cells = [s_cells s_cells(nc_ind)];
            % update cross-reference variables
            s_cells(nc_ind).sister_ParticleID = ParticleID;
            s_cells(end).sister_ParticleID = s_cells(nc_ind).sister_ParticleID;
            s_cells(nc_ind).sister_Index = numel(s_cells);
            s_cells(end).sister_Index = nc_ind;
            % redefine index variables
            nc_ind = numel(s_cells);
            % reset particle fields
            n_entries = numel(s_cells(nc_ind).xPosParticle);
            s_cells(nc_ind).ParticleID = NaN;
            s_cells(nc_ind).xPosParticle = NaN(1,n_entries);
            s_cells(nc_ind).yPosParticle = NaN(1,n_entries);
            s_cells(nc_ind).zPosParticle = NaN(1,n_entries);
            if threeD_flag
                s_cells(nc_ind).xPosParticle3D = NaN(1,sum(nc_filter));
                s_cells(nc_ind).yPosParticle3D = NaN(1,sum(nc_filter));
                s_cells(nc_ind).zPosParticle3D = NaN(1,sum(nc_filter));
                s_cells(nc_ind).fluo3D = NaN(1,n_entries);
            end
            s_cells(nc_ind).fluo = NaN(1,n_entries);    
            s_cells(nc_ind).fluoOffset = NaN(1,n_entries);  
        end
        s_cells(nc_ind).ParticleID = ParticleID;
        % find overlap between nucleus and trace
        nc_frames = s_cells(nc_ind).frames;         
        spot_filter = ismember(nc_frames,frames_full);            
        s_cells(nc_ind).spot_frames = spot_filter;
        if sum(spot_filter) < numel(frames_full)            
            error('Inconsistent particle and nucleus frames')
        end
        % record fluorescence info             
        s_cells(nc_ind).fluo(spot_filter) = trace_full;  
        % add offset info
        s_cells(nc_ind).fluoOffset(ismember(nc_frames,raw_pt_frames)) = ...
            cp_particles(j).Off(ismember(raw_pt_frames,nc_frames));
        % x, y, and z info                                
        s_cells(nc_ind).xPosParticle(ismember(nc_frames,raw_pt_frames)) = ...
            cp_particles(j).xPos(ismember(raw_pt_frames,nc_frames));
        s_cells(nc_ind).yPosParticle(ismember(nc_frames,raw_pt_frames)) = ...
            cp_particles(j).yPos(ismember(raw_pt_frames,nc_frames));   
        s_cells(nc_ind).zPosParticle(ismember(nc_frames,raw_pt_frames)) = ...
            cp_particles(j).zPos(ismember(raw_pt_frames,nc_frames));
        
        % 3D info
        if threeD_flag
            s_cells(nc_ind).xPosParticle3D(ismember(nc_frames,raw_pt_frames)) = ...
                cp_particles(j).xPosGauss3D(ismember(raw_pt_frames,nc_frames));
            
            s_cells(nc_ind).yPosParticle3D(ismember(nc_frames,raw_pt_frames)) = ...
                cp_particles(j).yPosGauss3D(ismember(raw_pt_frames,nc_frames));  
            
            s_cells(nc_ind).zPosParticle3D(ismember(nc_frames,raw_pt_frames)) = ...
                cp_particles(j).zPosGauss3D(ismember(raw_pt_frames,nc_frames));
            
            s_cells(nc_ind).fluo3D(ismember(nc_frames,raw_pt_frames)) = ...
                cp_particles(j).Fluo3DRaw(ismember(raw_pt_frames,nc_frames));
        end
        % add qc info
        s_cells(nc_ind).N = nDP;
        s_cells(nc_ind).sparsity = sparsity;        
        s_cells(nc_ind).qc_flag = qc_flag;         
    end      
    nucleus_struct = [nucleus_struct  s_cells];        
end
% add fields related to 2 spot analyses
for i = 1:numel(nucleus_struct)
    nucleus_struct(i).two_spot_flag = two_spot_flag;
    nucleus_struct(i).threeD_flag = threeD_flag;
    nucleus_struct(i).target_locus_flag = NaN;
    nucleus_struct(i).control_locus_flag = NaN;
end

disp('interpolating data...')
% generate interpolation fields
if threeD_flag
    interp_fields = {'fluo','fluo3D'};
else
    interp_fields = {'fluo'};
end
interpGrid = 0:TresInterp:60*60;
for i = 1:numel(nucleus_struct)
    time_vec = nucleus_struct(i).time;
    fluo_vec = nucleus_struct(i).fluo;
    start_i = find(~isnan(fluo_vec),1);
    stop_i = find(~isnan(fluo_vec),1,'last');    
    time_vec = time_vec(start_i:stop_i);       
    if nucleus_struct(i).qc_flag == 1
        time_interp = interpGrid(interpGrid>=time_vec(1)&interpGrid<=time_vec(end));
        for  j = 1:numel(interp_fields)
            vec = nucleus_struct(i).(interp_fields{j})(start_i:stop_i);
            %Look for clusters of 6 or more NaNs
            kernel = [1,1,1,1,1];
            vec_nans = isnan(vec);
            vec_conv = conv(kernel,vec_nans);
            vec_conv = vec_conv(3:end-2);
            z_ids = find(vec_conv==numel(kernel));
            z_ids = unique([z_ids-1 z_ids z_ids+1]); % get set of z_ids    
            vec(z_ids) = 0; % set clusters to zeros    
            vec(vec<0) = 0; % deal with negative values    
            % find single dp "blips". These will be replaced via interpolation
            % interpolate remaining NaNs    
            query_points = time_vec(isnan(vec));
            interp_t = time_vec(~isnan(vec));
            interp_f = vec(~isnan(vec)); 
            new_f = interp1(interp_t,interp_f,query_points);      
            vec(ismember(time_vec,query_points)) = new_f;        
            % Interpolate to standardize spacing    
            trace_interp = interp1(time_vec,vec,time_interp);                                    
            nucleus_struct(i).([interp_fields{j} '_interp']) = interp1(time_vec,vec,time_interp);
        end
    else
        time_interp = nucleus_struct(i).time;
        for  j = 1:numel(interp_fields)
            vec = nucleus_struct(i).(interp_fields{j});
            nucleus_struct(i).([interp_fields{j} '_interp']) = vec(1);
        end
    end           
    nucleus_struct(i).time_interp = time_interp;
    nucleus_struct(i).TresInterp = TresInterp;
end
% call function to calculate average psf difs
% save
save(nucleus_name ,'nucleus_struct') 

disp('calculating psf dims...')
calculate_average_psf(project,DropboxFolder);

disp('done.')
