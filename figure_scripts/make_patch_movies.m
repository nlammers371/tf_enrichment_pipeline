function make_patch_movies(Prefix, PreProcPath, PostProcPath, nc, varargin)
%
% Script to generate movies with overlaid colored patches indicating mean activity
% or EVER ON status
%
%author: Nick Lammers
%
%
%TO DO: 
% 1. Add instantaneous fluorescence visualization
% 2. Scale spot channel and opacity better 
%%% Load Data & Create Write Paths

close all;
au_to_polII = 13;
alpha = 1.4;
w = 7;
meanRate = 1;
duration = 0;
on_only = 0; % if 1,  only indicate whether nucleus turns on. 0 Mean rate patches
visible = 0; % if 1, shows each frame
saveFigs = 1;
last_time = NaN;
maxDuration = 30*3;
for i = 1:length(varargin)
    if strcmpi(varargin{i}, 'on_only')
        on_only = 1;
        meanRate = 0;
    elseif strcmpi(varargin{i}, 'duration')
        duration = 1;
        meanRate = 0;
    elseif strcmpi(varargin{i}, 'visible')
        visible = 1;
    elseif strcmpi(varargin{i}, 'noSave')
        saveFigs = 0;
    elseif strcmpi(varargin{i}, 'last_time')
        last_time = varargin{i+1};
    end
end

% set write path 
currentFolder = pwd;
slashes = strfind(currentFolder,'\');
fName = currentFolder(slashes(end)+1:end);
figPath = ['../../../../figures/' fName '/movies/'];
mkdir(figPath);
% [~,~,PrefixDropboxFolder,~, PreProcPath,...
%     ~, Prefix, ~,~,~,~, ~] = readMovieDatabase(Prefix);
if duration
    fill_color = [115 143 193]/256 * 1.1; % patch color
elseif on_only
    fill_color = [234 194 100]/256; % patch color
elseif meanRate
    fill_color = [85 169 116]/256; % patch color
end

% make write path
PostProcPath = [PostProcPath '/' Prefix '/'];
PreProcPath = [PreProcPath '/' Prefix '/'];

if on_only
    writePath = [figPath '/FractionOnFrames/'];
    save_prefix = 'fraction_on';
elseif meanRate
    writePath = [figPath '/MeanRateFrames/'];
    save_prefix = 'mean_rate';
elseif duration
    writePath = [figPath '/DurationFrames/'];
    save_prefix = 'duration';    
end
mkdir(writePath);


%Load the data
load([PostProcPath,'\CompiledParticles.mat']);
load([PostProcPath,'\' Prefix,'_lin.mat'], 'schnitzcells');
load([PostProcPath,'\Ellipses.mat'], 'Ellipses');
load([PostProcPath,'\FrameInfo.mat'], 'FrameInfo');

%%% calculate additional (derivative) movie parameters
PixelSize = FrameInfo(1).PixelSize;
MaxRadius = 5 / PixelSize; % um
xDim = FrameInfo(1).PixelsPerLine;
yDim = FrameInfo(1).LinesPerFrame;
zSlices = FrameInfo(1).NumberSlices;

[px, py] = meshgrid(1:xDim,1:yDim);
channel = 1; %this can be modified for multiple channels if needed

if iscell(AllTracesVector)
    AllTracesVector = AllTracesVector{channel};
end
if iscell(CompiledParticles)
    CompiledParticles = CompiledParticles{channel};
end

MaxmRNAFluo = prctile(nanmean(AllTracesVector),99);    
first_frame = eval(['nc' num2str(nc)]);

% make colorbar (save separately)
gradient = repmat(linspace(0,1,128)',1,3);
cmap = gradient .* fill_color;
fig = figure;
colormap(cmap);
if meanRate
    alpha_factor =  w / (w - .5*1.4);
    scale_factor = 1 / au_to_polII * alpha_factor / (w * 20) * 60;
    caxis([0 round(MaxmRNAFluo*scale_factor / 5)*5])
elseif duration
    caxis([0 maxDuration/3])
end
colorbar
axis off
saveas(fig,[writePath 'colorbar.pdf'])



if ~isnan(last_time)
    last_frame = find(floor(ElapsedTime-ElapsedTime(first_frame)) == last_time,1);
elseif nc == 14
    last_frame = numel(ElapsedTime);
else
    try
        last_frame = eval(['nc' num2str(nc+1)]);
    catch
        last_frame = numel(ElapsedTime);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate Movies %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if visible
    OverlayFig = figure();
else
    OverlayFig = figure('Visible','off');
end

OverlayAxes = axes(OverlayFig);


FrameRange=first_frame:last_frame;

%Iterate Through Frames
for Dummy = FrameRange    
    CurrentFrame = max(FrameRange) - Dummy + first_frame;
   
    %Track pixel assignments
    NucleusStateMat = zeros(yDim,xDim,3);
    NucleusIDMat = zeros(yDim,xDim);
    NucleusDistMat = ones(yDim,xDim)*MaxRadius;
    %Loop through ALL nuclei (inactive and active) and assign patches
    for s = 1:length(schnitzcells)
        MaxFrame = max(schnitzcells(s).frames);
        MinFrame = min(schnitzcells(s).frames);  
        % fix issue of patfches dropping off late
        if any(ismember(MaxFrame,120:175)) && CurrentFrame > MaxFrame
            cf = MaxFrame;
        else
            cf = CurrentFrame;
        end
        CurrentEllipse= schnitzcells(s).cellno(...
                        schnitzcells(s).frames==...
                        cf);  
        CurrentEllipse= schnitzcells(s).cellno(...
                        schnitzcells(s).frames==...
                        CurrentFrame);  
        x =  Ellipses{CurrentFrame}(CurrentEllipse,1)+1;
        y =  Ellipses{CurrentFrame}(CurrentEllipse,2)+1;
        if any(ismember(MinFrame,120:175)) 
            continue
        elseif isempty(x)
            continue
        end 
        distances = ((px-x).^2 + (py-y).^2).^.5;
        candidate_indices = NucleusDistMat > distances; 
        %Record Fluorescence        
        NucleusIDMat(candidate_indices) = s;
        NucleusDistMat(candidate_indices) = distances(candidate_indices);
    end
    %Loop Through Particles to see which remain on and how long these have
    %been active
    ParticlesToShow = [];
    for i = 1:length(CompiledParticles)        
        cp_frames = CompiledParticles(i).Frame;        
        all_frames = min(cp_frames):max(cp_frames);    
        %the patch turns on when the spot is first spotted and then stays
        %on for the rest of the movie. 
        extant_frames = all_frames((all_frames <= CurrentFrame)&...
            (all_frames>first_frame)); 
        
        if ~isempty(extant_frames)               
            filter = NucleusIDMat==CompiledParticles(i).Nucleus;
            frame_filter = ismember(cp_frames,extant_frames);
            
            %this calculation could use some clarification

            metric = nanmean(CompiledParticles(i).Fluo(frame_filter)) ;
            metric = min(MaxmRNAFluo,metric) / MaxmRNAFluo;

            if on_only
                metric = .85;
            elseif duration                
                metric = min(1,numel(extant_frames) / maxDuration);
            end
            for k = 1:3
                slice = NucleusStateMat(:,:,k);                        
                slice(filter) = fill_color(k)*metric;
                NucleusStateMat(:,:,k) = slice;
            end                    
            ParticlesToShow = [ParticlesToShow CompiledParticles(i).Nucleus];            
        end
    end
    %Now Draw Nucleus Borders for active nuclei
    NucleusBorderMat = zeros(size(NucleusIDMat));
    window = 1; %radius of convolution window
    for i = ParticlesToShow
        %get coordinates of nucleus patch
        if sum(sum(NucleusIDMat==i)) > 0
            x_vec = reshape(px(NucleusIDMat==i),[],1);
            y_vec = reshape(py(NucleusIDMat==i),[],1);
            for j = 1:length(x_vec)
                neighborhoodScore = sum(sum(NucleusIDMat(max(1,y_vec(j)-window):min(y_vec(j)+window,yDim),...
                         max(1,x_vec(j)-window):min(x_vec(j) + window,xDim))));
                if neighborhoodScore~= i*(2*window+1)^2
                    NucleusBorderMat(y_vec(j),x_vec(j)) = 1;
                end
            end
        end        
    end
    %Prevent overlap between fluroescence mask and borders
    for k = 1:3
        slice = NucleusStateMat(:,:,k);
        slice(NucleusBorderMat>0) = 0;
        NucleusStateMat(:,:,k) = slice;
    end

    %Make a maximum projection of the mRNA channel
    D=dir([PreProcPath,filesep,Prefix,'_',num2str(CurrentFrame,'%03d'),'_z*.tif']);
    %Do not load the first and last frame as they are black
    ImageTemp=zeros(yDim, xDim, zSlices) ;
    for i=2:(length(D)-1)
        ImageTemp(:,:,i-1)=imread([PreProcPath,D(i).name]);
    end
    mRNAImageRaw=max(ImageTemp,[],3);    
    %Load the corresponding histone image
    HistoneImage=imread([PreProcPath,filesep,...
        Prefix,'-His_',num2str(CurrentFrame,'%03d'),'.tif']);        

    %Overlay all channels
    MCPshading = NucleusStateMat(:,:,2);  %this is the patch. should be dynamically scaled. 
    scale = max(MaxmRNAFluo) / 2;
    MCPChannel =  mRNAImageRaw/scale + MCPshading; %this is the spot on top of the patch. 150 is the number that made the spot 
    %visible over the patch. 
    MCPChannel(MCPChannel>1) = 1;
         
    
    HistoneChannel=  mat2gray(HistoneImage)/2 + NucleusStateMat(:,:,1);%
    HistoneChannel(HistoneChannel>1) = 1;
    StateChannel = NucleusStateMat(:,:,3);    
    
    ImOverlay=cat(3,HistoneChannel,MCPChannel,StateChannel);
    imshow(fliplr(ImOverlay),'DisplayRange',[], 'Parent', OverlayAxes);
    ylim([16,yDim-16])
    if ceil(ElapsedTime(CurrentFrame)-ElapsedTime(FrameRange(1))) == 40        
        saveas(gcf,[writePath save_prefix '_patch_overlay_40min.tif']); 
    end
           
    text(OverlayAxes, 446,225,[num2str(round(ElapsedTime(CurrentFrame)-ElapsedTime(FrameRange(1))),'%02d'),...
        ' min'],'Color','k','FontSize',10,'BackgroundColor',[1,1,1,.5])	    %not sure about these numbers
    drawnow
    if saveFigs
        saveas(OverlayFig,[writePath '\nc',num2str(nc),...
            '-',num2str(CurrentFrame,'%02d'),'.tif']);   
    end
    
end
