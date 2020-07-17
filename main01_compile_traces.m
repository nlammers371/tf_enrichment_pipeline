% main01_compile_data(DataStatusTab,DropboxFolder,varargin)
%
% DESCRIPTION
% Funcion to compile relevant outputs from image analysis pipeline across
% multiple experiments
%
%
% ARGUMENTS
% DataStatusTab: master ID variable (should match a tab name in the Data Status
% sheet)
% DropboxFolder: full file path to folder containing compiled imaging
% results
%
% OPTIONS
% firstNC: script defaults to taking only nc14 traces. If you want
%           earlier traces, pass 'firstNC', followed by desired nuclear cycle
%           number
%
% OUTPUT
% nucleus_struct: compiled data set contain key nucleus and
% particle attributes

function nucleusStruct = main01_compile_traces(dataStatusTab,dropboxFolder,varargin)
addpath('./utilities')

% Set defaults
firstNC = 14;
minDP = 15;
pctSparsity = 50;
twoSpotFlag = contains(dataStatusTab, '2spot');
minTime = 0*60; % take no fluorescence data prior to this point
tresInterp = 20; 
calculatePSF = false;
projectName = dataStatusTab;

% Process input parameters
for i = 1:numel(varargin)
    if ischar(varargin{i}) && i < numel(varargin) && mod(i,2)==1
        eval([varargin{i} '=varargin{i+1};']);
    end
end

[rawResultsRoot, ~, ~] =   header_function(dropboxFolder, dataStatusTab);
[~, dataPath, ~] =   header_function(dropboxFolder, projectName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Set Path Specs, ID Vars %%%%%%%%%%%%%%%%%%%%%%%%

% Find the DataStatus.xlsx file and grab only those datasets marked as
% approved by 'ReadyForEnrichmentAnalysis' flag
approvedFlag = 'ReadyForEnrichmentAnalysis';
readyPrefixes = getProjectPrefixes(projectName,'customApproved',approvedFlag);
prefixes = readyPrefixes;
    
% Make the output filepath
mkdir(dataPath);
% Assign save names
nucleusName = [dataPath 'nucleus_struct.mat']; % names for compiled elipse struct

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Obtain Relevant Filepaths %%%%%%%%%%%%%%%%%%%%%%%

% Get all the file paths for the relevant pipeline outputs
numExperiments = numel(prefixes);
compiledParticlesFilenames = cell(1,numExperiments); % particles
compiledNucleiFilenames = cell(1,numExperiments); % nuclear protein
nucleiLineagesFilenames = cell(1,numExperiments); % nuclei lineages
frameInfoFilenames = cell(1,numExperiments); % movie frame info
for i = 1:numExperiments
    currDir = prefixes{i};            
    compiledNucleiFilenames{i} = [rawResultsRoot, currDir, filesep, 'CompiledNuclei.mat'];
    compiledParticlesFilenames = [rawResultsRoot, currDir, filesep, 'CompiledParticles.mat'];        
    nucleiLineagesFilenames = [rawResultsRoot, currDir, filesep, currDir, '_lin.mat'];           
    frameInfoFilenames = [rawResultsRoot, currDir, filesep, 'FrameInfo.mat'];        
end

% Generate set key data structure
setKey = array2table((1:numel(prefixes))','VariableNames',{'setID'});
setKey.prefix = prefixes';
save([dataPath 'set_key.mat'],'setKey')

% Generate master structure with info on all nuclei and traces in
% constituent sets
disp('Compiling data...')
nucleusStruct = [];
% Loop through filenames    
for i = 1:numExperiments
    % Try to read in raw files
    %nuclear lineages  
    try
        load(nucleiLineagesFilenames{i})              
    catch
        warning(['Failed to load nuclei lineages (''_lin.mat'') for prefix: ' prefixes{i} '. Skipping...' ])
        continue
    end
    %FrameInfo 
    try
        load(frameInfoFilenames{i})
    catch
        warning(['failed to load FrameInfo.mat files for prefix: ' prefixes{i} '. Skipping...' ])
        continue
    end    
    
    % MT 2020-07-15 I don't think this is ever used ... maybe remove?
%     processedNucleiData = load(compiledNucleiFilenames{i}); % processed nuclei 
    
    % Extract compiled particles structure
    processedSpotData = load(compiledParticlesFilenames{i});
    compiledParticles = processedSpotData.CompiledParticles;    
    if iscell(compiledParticles)
        compiledParticles = compiledParticles{1};
    end   
    
    % Check if certain info is present for this experiment    
    threeDFlag = isfield(compiledParticles,'Fluo3DRaw');    %3D spot fitting data        
    apFlag = isfield(compiledParticles, 'APpos');   %AP position info
        
    % Set the data set identifier
    setID = i;            
    % pull trace, time, and frame variables
    time_raw = processedSpotData.ElapsedTime*60; % time vector   
    traces_raw = processedSpotData.AllTracesVector; % array with a column for each trace 
    % check to see if traces are stored in cell array
    if iscell(traces_raw)
        traces_raw = traces_raw{1};
    end
    %%%%% Basic data characteristics %%%%%%
    frames_raw = 1:length(time_raw); % Frame list 
    % find index for first frame
    first_frame = max([1, processedSpotData.(['nc' num2str(firstNC)])]); 
    % filter trace mat and time
    traces_clean = traces_raw(first_frame:end,:);
    time_clean = time_raw(first_frame:end);    
    time_clean = time_clean - min(time_clean); % Normalize to start of nc
    % temporary fix for problematic time vector
    if strcmp(compiledParticlesFilenames{i},'E:\Nick\LivemRNA\Dropbox\LocalEnrichmentResults\2020-02-07-Dl_Ven_snaBAC_MCPmCh_F-F-F_Leica_Zoom2_7uW14uW_06/CompiledParticles.mat') %wtf
        time_clean = linspace(0,max(time_clean),numel(time_clean));
    end
    frames_clean = frames_raw(first_frame:end);        
    % get basic frame info
    yDim = FrameInfo(1).LinesPerFrame;
    xDim = FrameInfo(1).PixelsPerLine;
    zDim = FrameInfo(1).NumberSlices;
    zStep = FrameInfo(1).ZStep;
    
    %%%%% Compile nucleus schnitz info %%%%%%
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
            
            % Initialize particle fields--will be set to particle values 
            % for nuclei with matching particle
            s_cells(e_pass).ParticleID = NaN;
            s_cells(e_pass).xPosParticle = NaN(1,sum(nc_filter));
            s_cells(e_pass).yPosParticle = NaN(1,sum(nc_filter));
            s_cells(e_pass).zPosParticle = NaN(1,sum(nc_filter)); 
            s_cells(e_pass).APPosParticle = NaN(1,sum(nc_filter)); 
            s_cells(e_pass).spot_frames = NaN(1,sum(nc_filter)); 
            if threeDFlag
                s_cells(e_pass).xPosParticle3D = NaN(1,sum(nc_filter));
                s_cells(e_pass).yPosParticle3D = NaN(1,sum(nc_filter));
                s_cells(e_pass).zPosParticle3D = NaN(1,sum(nc_filter)); 
                s_cells(e_pass).fluo3D = NaN(1,sum(nc_filter));
            end
            s_cells(e_pass).fluo = NaN(1,sum(nc_filter)); 
            s_cells(e_pass).fluoOffset = NaN(1,sum(nc_filter));  
            
            % QC-related flags
            s_cells(e_pass).qc_flag = NaN;
            s_cells(e_pass).N = NaN;
            s_cells(e_pass).sparsity = NaN;
            
            % Add core nucleus info
            x = schnitzcells(e).cenx;            
            y = schnitzcells(e).ceny;                           
            s_cells(e_pass).xPos = x(nc_filter);
            s_cells(e_pass).yPos = y(nc_filter);
            % protein 
            s_cells(e_pass).raw_nc_protein = nanmax(schnitzcells(e).Fluo(nc_filter,:),[],2);
            s_cells(e_pass).frames = nc_frames';            
            s_cells(e_pass).Nucleus = e;     
            s_cells(e_pass).ncID = eval([num2str(setID) '.' sprintf('%04d',e)]);
            s_cells(e_pass).ncStart = firstNC;
            s_cells(e_pass).minDP = minDP;
            s_cells(e_pass).min_time = min_time;
            s_cells(e_pass).pctSparsity = pctSparsity;
            
            % Add time and set info
            s_cells(e_pass).time = time_clean(ismember(frames_clean,nc_frames));
            if numel(unique(s_cells(e_pass).time))~=numel(s_cells(e_pass).time)
                error('non-unique time values. Check FrameInfo')
            end
            s_cells(e_pass).setID = setID;
            fn = compiledParticlesFilenames{i}; % Get filename to store in struct  
            fn = fn(1:strfind(fn,'/')-1);
            s_cells(e_pass).source_path = fn;                                     
            % AP flag not relevant atm
            s_cells(e_pass).ap_flag = 0;
            
            % sister spot fields
            s_cells(e_pass).sister_Index = NaN;
            s_cells(e_pass).sister_ParticleID = NaN; 
            
            % increment
            e_pass = e_pass + 1;
        end
    end
    % Index vector to cross-ref w/ particles            
    e_index = [s_cells.Nucleus];          
   
    % Iterate through traces
    for j = 1:size(traces_clean,2)  
        % Raw fluo trace
        raw_trace = traces_clean(:,j); 
        % Get nucleus ID
        schnitz = compiledParticles(j).schnitz;    
        % Find first and last expression frames, requiring that spots
        % cannot apear earlier than some fixed time (min_time)
        trace_start = find(time_clean>=min_time&~isnan(raw_trace'),1);
        trace_stop = find(~isnan(raw_trace),1,'last');        
        %Creat versions with all intervening frames present (missing frames
        %appear as NaNs)
        trace_full = raw_trace(trace_start:trace_stop)';                                   
        trace_frames_full = frames_clean(trace_start:trace_stop); 
        % skip particles not in specified nc range
        if sum(~isnan(trace_full)) == 0
            continue
        end   
        
        % Find intersection btw full frame range and CP frames        
        raw_pt_frames = compiledParticles(j).Frame;
        raw_pt_frames = raw_pt_frames(ismember(raw_pt_frames,trace_frames_full));
        % Perform qc tests    
        nDP = sum(~isnan(trace_full));
        sparsity = prctile(diff(find(~isnan(trace_full))),pctSparsity);
        qc_flag = nDP >= minDP && sparsity == 1;     
        % trace-nucleus mapping may be many-to-1
        nc_ind = find(e_index==schnitz);  
        if length(nc_ind) ~= 1
            warning('Problem with Particle-Nucleus Crossref')
            continue
        end 
        
        % Identifier variable                        
        particle = compiledParticles(j).OriginalParticle;                        
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
            if threeDFlag
                s_cells(nc_ind).xPosParticle3D = NaN(1,sum(nc_filter));
                s_cells(nc_ind).yPosParticle3D = NaN(1,sum(nc_filter));
                s_cells(nc_ind).zPosParticle3D = NaN(1,sum(nc_filter));
                s_cells(nc_ind).fluo3D = NaN(1,n_entries);
            end
            s_cells(nc_ind).fluo = NaN(1,n_entries);    
            s_cells(nc_ind).fluoOffset = NaN(1,n_entries);  
        end
        % Record particle identifier
        s_cells(nc_ind).ParticleID = ParticleID;
        % find overlap between nucleus and trace
        nc_frames = s_cells(nc_ind).frames;         
        spot_filter = ismember(nc_frames,trace_frames_full);            
        s_cells(nc_ind).spot_frames = spot_filter;
        if sum(spot_filter) < numel(trace_frames_full)            
            error('Inconsistent particle and nucleus frames')
        end
        % record fluorescence info             
        s_cells(nc_ind).fluo(spot_filter) = trace_full; % note that
        % make filters
        nc_sp_ft1 = ismember(nc_frames,raw_pt_frames);
        nc_sp_ft2 = ismember(raw_pt_frames,nc_frames);
        % add offset info
        s_cells(nc_ind).fluoOffset(nc_sp_ft1) = compiledParticles(j).Off(nc_sp_ft2);
        % x, y, and z info                                
        s_cells(nc_ind).xPosParticle(nc_sp_ft1) = compiledParticles(j).xPos(nc_sp_ft2);
        s_cells(nc_ind).yPosParticle(nc_sp_ft1) = compiledParticles(j).yPos(nc_sp_ft2);   
        s_cells(nc_ind).zPosParticle(nc_sp_ft1) = compiledParticles(j).zPos(nc_sp_ft2);
        if apFlag
            s_cells(nc_ind).APPosParticle(nc_sp_ft1) = compiledParticles(j).APposParticle(nc_sp_ft2)*100;
        end
        % 3D info
        if threeDFlag
            s_cells(nc_ind).xPosParticle3D(nc_sp_ft1) = compiledParticles(j).xPosGauss3D(nc_sp_ft2);            
            s_cells(nc_ind).yPosParticle3D(nc_sp_ft1) = compiledParticles(j).yPosGauss3D(nc_sp_ft2);              
            s_cells(nc_ind).zPosParticle3D(nc_sp_ft1) = compiledParticles(j).zPosGauss3D(nc_sp_ft2);            
            s_cells(nc_ind).fluo3D(nc_sp_ft1) = compiledParticles(j).Fluo3DRaw(nc_sp_ft2);
        end
        % add qc info
        s_cells(nc_ind).N = nDP;
        s_cells(nc_ind).sparsity = sparsity;        
        s_cells(nc_ind).qc_flag = qc_flag;  
    end      
    nucleusStruct = [nucleusStruct  s_cells];        
end
% add additional fields
for i = 1:numel(nucleusStruct)
    nucleusStruct(i).two_spot_flag = two_spot_flag;
    nucleusStruct(i).threeD_flag = threeDFlag;
    nucleusStruct(i).target_locus_flag = NaN;
    nucleusStruct(i).control_locus_flag = NaN;
end    

disp('interpolating data...')
% generate interpolation fields
if threeDFlag
    interp_fields = {'fluo','fluo3D'};
else
    interp_fields = {'fluo'};
end
interpGrid = 0:TresInterp:60*60;
for i = 1:numel(nucleusStruct)
    time_vec = nucleusStruct(i).time;
    fluo_vec = nucleusStruct(i).fluo;
    start_i = find(~isnan(fluo_vec),1);
    stop_i = find(~isnan(fluo_vec),1,'last');    
    time_vec = time_vec(start_i:stop_i);       
    if numel(time_vec)>1%nucleus_struct(i).qc_flag == 1
        time_interp = interpGrid(interpGrid>=time_vec(1)&interpGrid<=time_vec(end));
        for  j = 1:numel(interp_fields)
            vec = nucleusStruct(i).(interp_fields{j})(start_i:stop_i);
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
            if ~isempty(query_points)
                new_f = interp1(interp_t,interp_f,query_points);  
                vec(ismember(time_vec,query_points)) = new_f;   
            else
                vec = interp_f;
            end                 
            % Interpolate to standardize spacing   
            try
                nucleusStruct(i).([interp_fields{j} '_interp']) = interp1(time_vec,vec,time_interp);
            catch
                error('why?')
            end
        end
    else
        time_interp = NaN;
        for  j = 1:numel(interp_fields)
            nucleusStruct(i).([interp_fields{j} '_interp']) = NaN;
        end
    end
    nucleusStruct(i).time_interp = time_interp;
    nucleusStruct(i).TresInterp = TresInterp;
end
% call function to calculate average psf difs
% save
save(nucleusName ,'nucleus_struct') 

if calculatePSF
    disp('calculating psf dims...')
    calculate_average_psf(projectName,DropboxFolder);
end

disp('done.')
