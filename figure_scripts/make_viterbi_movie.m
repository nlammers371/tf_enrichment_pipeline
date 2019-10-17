% Script to generate viterbi movie illustrating promoter dynamics

function make_viterbi_movie(project,varargin)

% set defaults
CurrentNC = 14;
K = 2; %number of states to display 

RawDataRoot = 'E:\LocalEnrichment\Data\PreProcessedData\';
setID = 1;
% project = 'Dl-Ven_snaBAC-mCh';
wInference = 7;
KInference = 3;

for i = 1:(numel(varargin)-1)  
    if ischar(varargin{i}) && ~strcmpi(varargin{i},'dropboxFolder')        
        eval([varargin{i} '=varargin{i+1};']);                
    elseif strcmpi(varargin{i},'dropboxFolder')
        dataPath = [varargin{i+1} '\ProcessedEnrichmentData\' project '/'];
    end
end

%Embryo ID and Folder location
PrefixDropboxRoot = 'E:\Nick\LivemRNA\Dropbox (Personal)\';
PrefixDropboxFolder = [PrefixDropboxRoot 'LocalEnrichmentResults\' Prefix '\'];
ProcessedDataRoot = [PrefixDropboxRoot 'ProcessedEnrichmentData\' project '\'];
FISHPath = [RawDataRoot Prefix];

%Load the data
load([PrefixDropboxFolder,'\CompiledParticles.mat'])
load([PrefixDropboxFolder,'\' Prefix,'_lin.mat'])
load([PrefixDropboxFolder,'\Ellipses.mat'])
load([PrefixDropboxFolder,'\FrameInfo.mat'])
load([ProcessedDataRoot 'nucleus_struct.mat'])
load([ProcessedDataRoot 'hmm_inference\w' num2str(wInference) '_K' num2str(KInference) '\soft_fit_struct.mat'])

% get image dims
xDim = FrameInfo(1).PixelsPerLine;
yDim = FrameInfo(1).LinesPerFrame;
MaxRadius = round(5/FrameInfo(1).PixelSize);
se = strel('disk',1);
% generate ref frame
[px, py] = meshgrid(1:xDim,1:yDim);
if iscell(MeanVectorAP)
    MeanVectorAP = MeanVectorAP{1};
end
FrameRange = nc14:length(MeanVectorAP);
TimeVec = (ElapsedTime - ElapsedTime(nc14))*60;
% set write path 
writePath = [PrefixDropboxRoot 'LocalEnrichmentFigures/hmm_movies/' Prefix '/'];
mkdir(writePath);

% Set write directory
if K == 2
    subdir = 'two_state_movie';   
elseif K == 3
    subdir = 'three_state_movie';    
end
StackPath = [writePath '\' subdir '\'];
mkdir(StackPath)

%%% Make some useful ref cells
tr_particle_vec = [nucleus_struct.ParticleID];
ss_particle_vec = soft_fit_struct.particle_index;
ss_set_indices = find(floor(ss_particle_vec)==setID);
nc_set_indices = find(floor(tr_particle_vec)==setID);
if K == 2
    color_cell = {[213 108 85],[122 169 116],[122 169 116]};
elseif K == 3
    color_cell = {[115 143 193]/yDim, [122 169 116]/yDim, [213 108 85]/yDim};
end

% generate lists of extant particles for each frame 
extant_particle_cell = cell(1,length(FrameRange));
viterbi_state_cell = cell(1,length(FrameRange));
        
for i = ss_set_indices
    ParticleID = ss_particle_vec(i);
    % generate inclusinve frame vector
    frames = nucleus_struct(tr_particle_vec==ParticleID).frames; 
    time_interp = nucleus_struct(tr_particle_vec==ParticleID).time_interp; 
    p_frames = frames(~isnan(nucleus_struct(tr_particle_vec==ParticleID).xPosParticle));
    p_frames = min(p_frames):max(p_frames);
    p_times = TimeVec(p_frames);
   
    % calculate and store effective states
    pt_ft = tr_particle_vec==ParticleID;
    Nucleus = floor((nucleus_struct(pt_ft).ncID - setID)*1e4);
    z_mat_raw = exp(soft_fit_struct.p_z_log_soft{ss_particle_vec==ParticleID});     
    z_mat_interp = interp1(time_interp',z_mat_raw',p_times,'nearest','extrap')';        
    if K == 2 && KInference == 3
        z_vec_mid = NaN(2,size(z_mat_interp,2));
        z_vec_mid(1,:) = z_mat_interp(1,:);
        z_vec_mid(2,:) = sum(z_mat_interp(2:3,:));
    else
        z_vec_mid = z_mat_interp;
    end
    [~, z_vec_int] = max(z_vec_mid);
    
    % update cell structures
    for p = 1:length(p_frames)
        f_filter = FrameRange==p_frames(p);
        if sum(f_filter) > 0 
            % extanmt particle list
            extant_list = extant_particle_cell{f_filter};
            extant_list = [extant_list ParticleID];
            extant_particle_cell{f_filter} = extant_list;
            % promoter state list
            state_list = viterbi_state_cell{f_filter};
            state_list = [state_list z_vec_int(p)];
            viterbi_state_cell{f_filter} = state_list;
        end
    end
end

%Iterate Through Frames
for CurrentFrame = FrameRange
    
    % In initialize arrays
    NucleusStateMat = zeros(yDim,xDim,3); 
    NucleusDistMat = zeros(yDim,xDim); 
    
    % Loop through ALL nuclei (inactive and active) and assign pixels
    % to a nucleus
    for s = 1:numel(nc_set_indices)
        ind = nc_set_indices(s);        
        frame_vec = nucleus_struct(ind).frames;
        x =  nucleus_struct(ind).xPos(frame_vec==CurrentFrame);
        y =  nucleus_struct(ind).yPos(frame_vec==CurrentFrame);        
        NucleusDistMat(y,x) = floor((nucleus_struct(ind).ncID - setID)*1e4);
    end
    [d_mat, idx] = bwdist(NucleusDistMat);
    NucleusIDMat = reshape(NucleusDistMat(idx),yDim,xDim);
    NucleusIDMat(d_mat > MaxRadius) = 0;
    % reference frame index cell
    extant_particles = extant_particle_cell{CurrentFrame-nc14 + 1};
    particle_states = viterbi_state_cell{CurrentFrame-nc14 + 1};
    ParticlesToShow = NaN(size(extant_particles));
    for i = 1:length(extant_particles)
        ParticleID = extant_particles(i);
        z_state = particle_states(i);
        pt_ft = tr_particle_vec==ParticleID;
        Nucleus = floor((nucleus_struct(pt_ft).ncID - setID)*1e4);            
        % Generate patch
        filter = NucleusIDMat==Nucleus;
        color_vec = fliplr(color_cell{z_state});
        for k = 1:3
            nc_slice = NucleusStateMat(:,:,k);
            nc_temp = zeros(size(nc_slice));
            nc_temp(filter) = ones;
            nc_temp = imerode(nc_temp,se);
            nc_slice(nc_temp==1) = color_vec(k);
            NucleusStateMat(:,:,k) = nc_slice;
        end
        ParticlesToShow(i) = Nucleus;
    end

    %Make a maximum projection of the mRNA channel
    n_char = length(num2str(CurrentFrame));
    z_string = '000';
    index_string = [z_string(1:3-n_char) num2str(CurrentFrame)];
    DmRNA=dir([FISHPath,'\',Prefix,'_',index_string,'_z*_ch02.tif']);
    %Do not load the first and last frame as they are black
    ImageTemp=[];
    for m=2:(length(DmRNA)-1)
        ImageTemp(:,:,m-1)=imread([FISHPath,'\',DmRNA(m).name]);
    end
    mRNAImage= imadjust(mat2gray(max(ImageTemp,[],3)));           
    mRNAImage = mat2gray(mRNAImage.* reshape([162 60 150],1,1,3));
    %Load the corresponding histone image
    DProtein=dir([FISHPath,'\',Prefix,'_',index_string,'_z*_ch01.tif']);
    %Do not load the first and last frame as they are black
    ImageTemp=[];
    for m=2:(length(DProtein)-1)
        ImageTemp(:,:,m-1)=imread([FISHPath,'\',DProtein(m).name]);
    end
    HistoneImage= mat2gray(median(ImageTemp,3));         
    HistoneImage = mat2gray(HistoneImage .* reshape([12 118 60],1,1,3));
       
    patch_stack = mat2gray(cat(3,NucleusStateMat(:,:,3),NucleusStateMat(:,:,2),NucleusStateMat(:,:,1)));
    
    ImOverlay= mat2gray(.75*HistoneImage + .5*mRNAImage + patch_stack);
    
    %%% Make figure
    OverlayFig = figure('Visible','off');
    clf
    imshow(fliplr(ImOverlay))   
    
    write_time = round(ElapsedTime(CurrentFrame)-ElapsedTime(FrameRange(1)));
    write_time_prev = round(ElapsedTime(CurrentFrame-1)-ElapsedTime(FrameRange(1)));
    
    xlim([0+MaxRadius,xDim-MaxRadius])    
    ylim([0+MaxRadius,yDim-MaxRadius])
    

    text(15+MaxRadius,yDim-30-MaxRadius,[sprintf('%02d',write_time),...
        ' min'],'Color','k','FontSize',10,'BackgroundColor',[1,1,1,.5])
    drawnow
    
    saveas(gcf,[StackPath '\nc',num2str(CurrentNC),...
        '-',sprintf('%03d',CurrentFrame), '.tif']);   
    close all
           
end
