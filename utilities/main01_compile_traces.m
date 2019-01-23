% main01_compile_traces(project, FolderPath, keyword, include_vec)
%
% DESCRIPTION
% Funcion to compile relevant outputs from image analysis pipeline across
% multiple experi

%
% ARGUMENTS
% project: master ID variable 
%
% FolderPath: Full or relative path to DynamicsResults folder (pr
%               equivalent)
%
% keyword: String contained in all folder names for projects one wishes to 
%          compile. For instance: 'Eve2MS2', will full all projects
%          containing this string
%
% include_vec:  Vector specifying IDs within list of projects
%               pulled by 'keyword' to keep. If passed as empty vector, all
%               matching projects will be taken
%
% OUTPUT: nucleus_struct: compiled data set

function nucleus_struct = main01_compile_traces(project,FolderPath,keyword,include_vec)
 
%--------------------------Set Path Specs, ID Vars------------------------%    
%%% folders
data_path = ['../../dat/' project '/']; % data mat directory
%%% make filepath
mkdir(data_path);
%%% assign save names
nucleus_name = [data_path 'nucleus_struct.mat']; % names for compiled elipse struct

%--------------------------Obtain Relevant Filepaths----------------------%
% obtain and store set names
dirinfo = dir([FolderPath '*' keyword '*']);
dirinfo(~[dirinfo.isdir]) = []; %remove non-directories
if isempty(include_vec)
    include_vec = 1:numel(dirinfo);
end
% initialize filename vectors
cp_filenames = {}; % particles
cn_filenames = {}; % protein data
ap_filenames = {}; % ap info
nc_filenames = {}; % nuclei
fov_filenames = {}; % fov info
sp_filenames = {}; % spot info
pa_filenames = {}; % particle info that contains z info
for d = 1:numel(include_vec)
    thisdir = dirinfo(include_vec(d)).name;            
    % append file paths
    cn_filenames = [cn_filenames {[thisdir '/CompiledNuclei.mat']}];
    cp_filenames = [cp_filenames {[thisdir '/CompiledParticles.mat']}];    
    ap_filenames = [ap_filenames {[thisdir '/APDetection.mat']}];    
    nc_filenames = [nc_filenames {[thisdir '/' thisdir '_lin.mat']}];           
    fov_filenames = [fov_filenames {[thisdir '/FrameInfo.mat']}];
    sp_filenames = [sp_filenames {[thisdir '/Spots.mat']}];
    pa_filenames = [pa_filenames {[thisdir '/Particles.mat']}];
end
%%% generate set key
set_key = array2table(include_vec','VariableNames',{'setID'});
dirs = {dirinfo(include_vec).name};
set_key.prefix = dirs';
save([data_path 'set_key.mat'],'set_key')

%%% Generate master structure with info on all nuclei and traces in
%%% constituent sets
nucleus_struct = [];
% Loop through filenames    
for i = 1:length(cp_filenames) 
    % read in raw files
    load([FolderPath ap_filenames{i}]) % AP Info   
    load([FolderPath nc_filenames{i}]) % Ellipse Info
    load([FolderPath fov_filenames{i}]) % FrameInfo Info
    load([FolderPath sp_filenames{i}]) % Spots info
    load([FolderPath pa_filenames{i}]); % Raw Particles
    raw_data = load([FolderPath cp_filenames{i}]); % Particles    
    protein_data = load([FolderPath cn_filenames{i}]); % Protein in Nuclei     
    % set identifier
    setID = include_vec(i);            
    % pull trace and nuclei variables
    time_raw = raw_data.ElapsedTime*60; % time vector            
    traces_raw = raw_data.AllTracesVector; % Array with a column for each trace    
    if iscell(traces_raw)
        traces_raw = traces_raw{1};
    end
    frames_raw = 1:length(time_raw); % Frame list    
    first_frame = raw_data.nc14; % Default to start of nc14 for the moment    
    % filter trace mat and time
    traces_clean = traces_raw(first_frame:end,:);
    time_clean = time_raw(first_frame:end);    
    time_clean = time_clean - min(time_clean); % Normalize to start of nc14    
    frames_clean = frames_raw(first_frame:end);    
    % ectract protein data
    cp_protein = protein_data.CompiledNuclei;   
    if iscell(cp_protein)
        cp_protein= cp_protein{1};
    end
    % compile schnitz info
    s_cells = struct;
    e_pass = 1;    
    for e = 1:length(schnitzcells)
        e_frames = schnitzcells(e).frames;
        nc14_filter = ismember(e_frames,frames_clean);
        nc14_frames = e_frames(nc14_filter);
        if length(nc14_frames) >= 1% skip nuclei not in nc14
            x = schnitzcells(e).cenx;            
            y = schnitzcells(e).ceny;                           
            s_cells(e_pass).xPos = x(nc14_filter);
            s_cells(e_pass).yPos = y(nc14_filter);             
            %Will be set to particle position for nuclei with matching
            %particle
            s_cells(e_pass).ParticleID = NaN;
            s_cells(e_pass).xPosParticle = NaN(1,sum(nc14_filter));
            s_cells(e_pass).yPosParticle = NaN(1,sum(nc14_filter));
            s_cells(e_pass).brightestZs = NaN(1,sum(nc14_filter));
            s_cells(e_pass).spot_z_cell = cell(1,sum(nc14_filter));
            % add dummy fluorescence field
            s_cells(e_pass).fluo = NaN(1,sum(nc14_filter));
            s_cells(e_pass).frames = nc14_frames';            
            s_cells(e_pass).Nucleus = e;                
            s_cells(e_pass).ap_vector = schnitzcells(e).APpos(nc14_filter);
            s_cells(e_pass).ncID = eval([num2str(setID) '.' sprintf('%04d',e)]);
            %Will be set to mean particle position for nuclei with matiching
            %particles
            s_cells(e_pass).xMean = mean(x);
            s_cells(e_pass).yMean = mean(y);
            % time and set info
            s_cells(e_pass).time = time_clean(ismember(frames_clean,nc14_frames));                        
            s_cells(e_pass).setID = setID;
            fn = cp_filenames{i}; % Get filename to store in struct  
            fn = fn(1:strfind(fn,'/')-1);
            s_cells(e_pass).source_path = fn;                        
            s_cells(e_pass).PixelSize = FrameInfo(1).PixelSize;                        
            % add protein info              
            prot_nuc = cp_protein([cp_protein.schnitz] == e);                     
            pt_vec = prot_nuc.FluoMax(nc14_filter);
            s_cells(e_pass).protein = pt_vec';
            e_pass = e_pass + 1;
        end
    end
    
    %%% Now add particle info
    
    % Index vector to cross-ref w/ particles            
    e_index = [s_cells.Nucleus]; 
    
    % Extract compiled particles structure       
    cp_particles = raw_data.CompiledParticles; 
    if iscell(cp_particles)
        cp_particles = cp_particles{1};
    end    
    % iterate through traces
    for j = 1:size(traces_clean,2)  
        % Raw fluo trace
        raw_trace = traces_clean(:,j); 
        % Get nucleus ID
        schnitz = cp_particles(j).schnitz;        
        % skip particles not in nc14
        if sum(~isnan(raw_trace)) == 0
            continue
        end
        trace_start = find(~isnan(raw_trace),1);
        trace_stop = find(~isnan(raw_trace),1,'last');        
        %Creat versions with all intervening frames present (missing frames
        %appear as NaNs
        trace_full = raw_trace(trace_start:trace_stop)';               
        time_full = time_clean(trace_start:trace_stop);                
        frames_full = frames_clean(trace_start:trace_stop);                       
        % Check for mapping issues btw AllTracesVector and
        % CompiledParticles        
        raw_pt_frames = cp_particles(j).Frame;
        pt_frames = raw_pt_frames(ismember(raw_pt_frames,frames_full));        
        
        % Initialize arrays to store Z info
        bZ = NaN(1, length(frames_full));
        z_cell = cell(1, length(frames_full));
        % Get spot Z location info
        part_i = cp_particles(j).OriginalParticle;        
        PT_Frames = Particles(part_i).Frame;
        for id = 1:length(pt_frames)            
            cor_idx = find(pt_frames(id) == frames_full);
            rel_id = raw_pt_frames==pt_frames(id);
            if Particles(part_i).Index(PT_Frames==pt_frames(id)) > length(Spots(pt_frames(id)) ...
                    .Fits)
                error('Mismatch between Particles and Spots dimensions')                                
            else
                bZ(cor_idx) = Spots(pt_frames(id)) ...
                        .Fits(Particles(part_i).Index(rel_id)).brightestZ;
                z_cell{cor_idx} = Spots(Particles(part_i).Frame(rel_id)) ...
                        .Fits(Particles(part_i).Index(rel_id)).z;
            end
        end
        % find corresponding nucleus index        
        nc_ind = find(e_index==schnitz);        
        if length(nc_ind) ~= 1
            error('Error: Problem with Particle-Nucleus Crossref')
        end                                         
        % find overlap between nucleus and trace
        nc_frames = s_cells(nc_ind).frames;         
        spot_filter = ismember(nc_frames,frames_full);            
        s_cells(nc_ind).spot_frames = spot_filter;
        if sum(s_cells(nc_ind).spot_frames) < numel(time_full)
            error('inconsistent particle and nucleus frames')
        end
        % record fluorescence info             
        s_cells(nc_ind).fluo(s_cells(nc_ind).spot_frames) = trace_full;               
        % x and y info                                
        s_cells(nc_ind).xPosParticle(ismember(nc_frames,pt_frames)) = ...
            cp_particles(j).xPos(ismember(pt_frames,nc_frames));
        s_cells(nc_ind).yPosParticle(ismember(nc_frames,pt_frames)) = ...
            cp_particles(j).yPos(ismember(pt_frames,nc_frames));
        % add z info                       
        s_cells(nc_ind).brightestZs(spot_filter) = bZ;   
        s_cells(nc_ind).spot_z_cell(spot_filter) = z_cell;
        % Identifier variables                        
        particle = cp_particles(j).OriginalParticle;                        
        s_cells(nc_ind).ParticleID = eval([num2str(setID) '.' sprintf('%04d',particle)]);        
    end      
    nucleus_struct = [nucleus_struct  s_cells];        
end
% save
save(nucleus_name ,'nucleus_struct') 

