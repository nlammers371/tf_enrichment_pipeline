% Script to investigate spatio-temporal dynamics of protein distributions
% within nuclei. Ultimate goal is to find a way to infer position of active
% loci using protein channel alone
clear
close all

rawPath = 'E:\LocalEnrichment\Data\PreProcessedData\';
dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
project = 'Dl_Venus_snaBAC_MCPmCherry_Leica_Zoom2_7uW14uW';
varargin = {'dropboxFolder',dropboxFolder};

proteinChannel = 1;
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
protein_dynamics_struct = struct;
nc_iter = 0;

plot_flags = false(size(nucleus_struct));
for i = 1:numel(nucleus_struct)
    fluo = nucleus_struct(i).fluo;
    duration = find(~isnan(fluo),1,'last')-find(~isnan(fluo),1);
    if ~isempty(duration)
        plot_flags(i) = sum(~isnan(fluo)) > 40 & (sum(~isnan(fluo)) / duration)>.8;
    end
end
% draw sample set
rng(123)
plot_indices = randsample(find(plot_flags),50,false);

%%% iterate
for i = plot_indices
    setID = nucleus_struct(i).setID;    
    ncID = nucleus_struct(i).ncID;    
    
    % get size params    
    nb_sz = nb_sizes(set_index==setID); % size of nucleus neighborhood to use
    int_kernel = roi_rad_spot_pix(set_index==setID); %sigma for gaussian smoothing kerne;
    roi_spot = roi_rad_spot_pix(set_index==setID);   
    % load and  MCP mCherry and protein stacks
    src = set_key(set_key.setID==setID,:).prefix{1};        
        
    xPosParticle = nucleus_struct(i).xPosParticle;
    nan_ft = ~isnan(xPosParticle);
    pt_ft = 1:find(nan_ft,1,'last');
    frameVec = nucleus_struct(i).frames(nan_ft);    
    frameVecFull = nucleus_struct(i).frames(pt_ft);    
    
    % for now limit ourselves to nearly complete traces
    if sum(nan_ft) / numel(frameVec) < .8 || numel(frameVec) < 20
        continue
    end    
    
    timeVec = nucleus_struct(i).time(pt_ft);    
    
    fluoVec = nucleus_struct(i).fluo(find(pt_ft,1):find(pt_ft,1,'last'));
    
    xPosParticle = xPosParticle(pt_ft);
    yPosParticle = nucleus_struct(i).yPosParticle(pt_ft);
    zPosParticle = nucleus_struct(i).zPosParticle(pt_ft);
    xPosNC = nucleus_struct(i).xPos(pt_ft);
    yPosNC = nucleus_struct(i).yPos(pt_ft);       
    
    if min([xPosNC,yPosNC]) < nb_sz + 1 || max([xPosNC,yPosNC]) > 856 - nb_sz -1
        continue
    end
    nc_iter = nc_iter+1;    
%     % interpolate
%     xPosParticle = xPostParticle;%interp1(timeVec,xPosParticle,timeVecFull);      
%     yPosParticle = yPosParticle;%interp1(timeVec,yPosParticle,timeVecFull);
%     zPosParticle = zPosParticle;%round(interp1(timeVec,zPosParticle,timeVecFull));
% %     fluoVecFull = round(imgaussfilt(interp1(timeVec,fluoVec,timeVecFull),1));
%     xPosNCFull = xPosNC; %round(interp1(timeVec,xPosNC,timeVecFull),1);
%     yPosNCFull = yPosNC;% round(interp1(timeVec,yPosNC,timeVecFull),1);
    % make filepath 
    ncPath = [snipPath '/nc_' num2str(1e4*ncID) '/'];
    mkdir(ncPath)
    % record fields 
    protein_dynamics_struct(nc_iter).xPosParticle = xPosParticle;      
    protein_dynamics_struct(nc_iter).yPosParticle = yPosParticle;      
    protein_dynamics_struct(nc_iter).zPosParticle = zPosParticle;      
    protein_dynamics_struct(nc_iter).xPosNC = xPosNC;
    protein_dynamics_struct(nc_iter).yPosNC = yPosNC;
    protein_dynamics_struct(nc_iter).timeVec = timeVec;
    
    fNorm = fluoVec / nanmax(fluoVec) * 60;
    fNorm(fNorm<0) = NaN;
    nc_snip_stack = NaN(2*nb_sz+1,2*nb_sz+1,numel(zPosParticle));  
    particleFrames = frameVecFull;
    particleFrames(~nan_ft) = NaN;
    for j = 1:numel(frameVecFull)
        x_nucleus = xPosNC(j);
        y_nucleus = yPosNC(j);
        frame = frameVecFull(j);
        [~,mi] = nanmin(abs(particleFrames-frame));        
        fileList = dir([rawPath src '/*_' sprintf('%03d',frame)  '*z' sprintf('%02d',zPosParticle(mi)) '_ch0' num2str(proteinChannel) '.tif']);
        fileName = [fileList(1).folder '/' fileList(1).name];
        protein_frame = double(imread(fileName));
        protein_snip = protein_frame(y_nucleus-nb_sz:y_nucleus+nb_sz,x_nucleus-nb_sz:x_nucleus+nb_sz);
        nc_snip_stack(:,:,j) = protein_snip;
    end
    protein_dynamics_struct(nc_iter).nc_snip_stack = nc_snip_stack;   
    % find percentile rank of local concentration over time
    prctile_vec = NaN(1,numel(frameVecFull));
    nc_snip_stack_sm = NaN(size(nc_snip_stack));
    for j = 1:numel(frameVecFull)
        nc_snip_stack_sm(:,:,j) = imgaussfilt(nc_snip_stack(:,:,j),roi_spot/1.5);
        % approximate mast for now 
        dist_mat = zeros(size(nc_snip_stack_sm(:,:,j)));
        dist_mat(nb_sz+1,nb_sz+1) = 1;
        dist_mat = bwdist(dist_mat);
        mask = dist_mat <=25;
        yp = round(yPosParticle(j)-yPosNC(j)) + nb_sz + 1;
        xp = round(xPosParticle(j)-xPosNC(j)) + nb_sz + 1;
        masked = nc_snip_stack_sm(:,:,j);
        masked(~mask) = NaN;
        try
            spot_val = masked(yp,xp);
        catch
            continue
        end
        prctile_vec(j) = sum(masked(:)<spot_val) / sum(~isnan(masked(:)));
    end
    protein_dynamics_struct(nc_iter).nc_snip_stack_sm = nc_snip_stack_sm;
    protein_dynamics_struct(nc_iter).prctile_vec = prctile_vec;

    max_vec = imgaussfilt(reshape(max(max(nc_snip_stack_sm,[],2)),[],1),5);
    min_vec = imgaussfilt(reshape(min(min(nc_snip_stack_sm,[],2)),[],1),5);            
    % plot intensity fluctuations over time
    for j = 1:size(nc_snip_stack_sm,3)
        pt_field_fig = figure('Visible','off');
        imagesc(nc_snip_stack_sm(:,:,j))
        caxis([min_vec(j) max_vec(j)])
        hold on
        scatter(xPosParticle(j)-xPosNC(j) + nb_sz + 1,yPosParticle(j)-yPosNC(j)+ nb_sz + 1,...
            20+fNorm(j),'MarkerFaceColor','red','MarkerFaceAlpha',.5)
        h = colorbar;
        ylabel(h,'Dl concentration (au)')
        saveas(pt_field_fig,[ncPath 'frame_' sprintf('%03d',frameVecFull(j)) '.tif'])
        close all
    end                             
end
save([dataPath 'protein_dynamics_struct.mat'],'protein_dynamics_struct')