% Script to investigate spatio-temporal dynamics of protein distributions
% within nuclei. Ultimate goal is to find a way to infer position of active
% loci using protein channel alone
clear
close all

rawPath = 'E:\LocalEnrichment\Data\PreProcessedData\';
dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
project = 'Dl_Venus_snaBAC_mCherry_Leica_hp';
varargin = {'dropboxFolder',dropboxFolder};

proteinChannel = 1;
zeissFlag = 0;
ROIRadiusSpot = .2; % radus (um) of region used to query and compare TF concentrations
mfTolerance = ROIRadiusSpot;
minSampleSep = 1; %um
dataPath = ['../dat/' project '/'];
figPath = ['../fig/' project '/'];
for i = 1:numel(varargin)  
    if ischar(varargin{i})
        if ismember(varargin{i},{'zeissFlag','dropboxFolder','ROIRadiusSpot','ROIRadiusControl','minSampleSep'})       
            eval([varargin{i} '=varargin{i+1};']);
        end
        if strcmpi(varargin{i},'dropboxFolder')
            dataPath = [varargin{i+1} '\ProcessedEnrichmentData\' project '/'];
            figPath = [varargin{i+1} '\LocalEnrichmentFigures\' project '/'];
        end
    end
end

% Load trace data
load([dataPath '/nucleus_struct.mat'],'nucleus_struct')
load([dataPath '/set_key.mat'],'set_key')
snipPath = [figPath '/nc_intensity_images/'];
mkdir(snipPath)
addpath('./utilities')
% get MCP channel
options = [1 2];
mcp_channel = options(options~=proteinChannel);


%%% set snip size to use
set_vec = [nucleus_struct.setID];
set_index = unique(set_vec);
px_sizes = [];
for s = 1:numel(set_index)
    px = [nucleus_struct(set_vec==set_index(s)).PixelSize];
    px_sizes = [px_sizes px(1)];
end% determine size of neighborhood to use
nb_sizes = round(3.5 ./ px_sizes);
% calculate ROI size in pixels for spot and control
roi_rad_spot_pix = round(ROIRadiusSpot ./ px_sizes);
% initialize structure to store results
nc_snip_structure = struct;
%%% iterate
for i = 1:numel(nucleus_struct)
    setID = nucleus_struct(i).setID;    
    ncID = nucleus_struct(i).ncID;    
    
    % get size params    
    nb_sz = nb_sizes(set_index==setID); % size of nucleus neighborhood to use
    int_kernel = roi_rad_spot_pix(set_index==setID); %sigma for gaussian smoothing kerne;
    roi_spot = roi_rad_spot_pix(set_index==setID);   
    % load and  MCP mCherry and protein stacks
    src = set_key(set_key.setID==setID,:).prefix{1};        
        
    xPosParticle = nucleus_struct(i).xPosParticle;
    pt_ft = ~isnan(xPosParticle);
    frameVec = nucleus_struct(i).frames(pt_ft);    
    frameVecFull = nucleus_struct(i).frames;
    frameVecFull = frameVecFull(find(pt_ft,1):find(pt_ft,1,'last'));
    
    % for now limit ourselves to nearly complete traces
    if sum(pt_ft) / numel(frameVecFull) < .8 || numel(frameVec) < 20
        continue
    end
    
    
    timeVecFull = nucleus_struct(i).time;
    timeVec = timeVecFull(pt_ft);
    timeVecFull = timeVecFull(find(pt_ft,1):find(pt_ft,1,'last'));
    
    fluoVec = nucleus_struct(i).fluo(find(pt_ft,1):find(pt_ft,1,'last'));
    
    xPosParticle = xPosParticle(pt_ft);
    yPosParticle = nucleus_struct(i).yPosParticle(pt_ft);
    zPosParticle = nucleus_struct(i).brightestZs(pt_ft);
    xPosNC = nucleus_struct(i).xPos(pt_ft);
    yPosNC = nucleus_struct(i).yPos(pt_ft);       
    
    if min([xPosNC,yPosNC]) < nb_sz + 1 || max([xPosNC,yPosNC]) > 856 - nb_sz -1
        continue
    end
        
    % interpolate
    xPosParticleFull = imgaussfilt(interp1(timeVec,xPosParticle,timeVecFull),1);      
    yPosParticleFull = imgaussfilt(interp1(timeVec,yPosParticle,timeVecFull),1);
    zPosParticleFull = round(imgaussfilt(interp1(timeVec,zPosParticle,timeVecFull),1));
%     fluoVecFull = round(imgaussfilt(interp1(timeVec,fluoVec,timeVecFull),1));
    xPosNCFull = round(imgaussfilt(interp1(timeVec,xPosNC,timeVecFull),1));
    yPosNCFull = round(imgaussfilt(interp1(timeVec,yPosNC,timeVecFull),1));
    % make filepath 
    ncPath = [snipPath '/nc_' num2str(ncID) '/'];
    mkdir(ncPath)
    
    fNorm = fluoVec / nanmax(fluoVec) * 50;
    fNorm(fNorm<0) = NaN;
    nc_snip_stack = NaN(2*nb_sz+1,2*nb_sz+1,numel(zPosParticleFull));    
    for j = 1:numel(frameVecFull)
        x_nucleus = xPosNCFull(j);
        y_nucleus = yPosNCFull(j);
        
        frame = frameVecFull(j);
        fileList = dir([rawPath src '/*_' sprintf('%03d',frame)  '*z' sprintf('%02d',zPosParticleFull(j)) '_ch0' num2str(proteinChannel) '.tif']);
        fileName = [fileList(1).folder '/' fileList(1).name];
        protein_frame = double(imread(fileName));
        protein_snip = protein_frame(y_nucleus-nb_sz:y_nucleus+nb_sz,x_nucleus-nb_sz:x_nucleus+nb_sz);
        nc_snip_stack(:,:,j) = protein_snip;
    end
    
    nc_snip_stack_sm = imgaussfilt3(nc_snip_stack,roi_spot/1.5);
    
    max_vec = imgaussfilt(reshape(max(max(nc_snip_stack_sm,[],2)),[],1),5);
    min_vec = imgaussfilt(reshape(min(min(nc_snip_stack_sm,[],2)),[],1),5);        
    
    % plot intensity fluctuations over time
    for j = 1:size(nc_snip_stack_sm,3)
        pt_field_fig = figure('Visible','off');
        imagesc(nc_snip_stack_sm(:,:,j))
        caxis([min_vec(j) max_vec(j)])
        hold on
        scatter(xPosParticleFull(j)-xPosNCFull(j) + nb_sz + 1,yPosParticleFull(j)-yPosNCFull(j)+ nb_sz + 1,...
            fNorm(j),'MarkerFaceColor','red','MarkerFaceAlpha',.3)
        h = colorbar;
        ylabel(h,'Dl concentration (au)')
        saveas(pt_field_fig,[ncPath 'frame_' sprintf('%03d',frameVecFull(j)) '.tif'])
        close all
    end    
                         
end