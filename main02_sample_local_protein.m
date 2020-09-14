% main02_sample_local_protein(projectName, varargin)
%
% DESCRIPTION
% Script to generate samples of local protein contentration at active gene
% loci and selected control locations
%
% ARGUMENTS
% projectName: master ID variable (should match a tab name in the 
%              DataStatus.xlsx spreadsheet)
%
% OPTIONS
% script allows any default variable to be set using format:
%       "VariableNameString", VariableValue
%
% OUTPUT: spot_struct_protein: compiled data set with protein samples

function spot_struct_protein = main02_sample_local_protein(projectName,varargin)

%% %%%%%%%%%%%%%%%%%%%%%%% Set Defaults %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('utilities'))
ROIRadiusSpot_um = .2; % radus (um) of region used to query and compare TF concentrations
minSampleSep_um = 1.5; %um
minEdgeSep_um = .25; %um
snippet_size_um = 1.5;
min_nucleus_radius_um = 2;
max_nucleus_radius_um = 4;

segmentNuclei = 0;
use3DSpotInfo = 1;

% PSF info for 3D sampling
use_psf_fit_dims = false; % NL: currently no supported
xy_sigma_um = 0.25;% um 
xy_sigma_nuclear_um = 1.5;
z_sigma_um = 0.6; % um
ignoreQC = false;
write_snip_flag = false; %NL: what does this do?

%% %%%%%%%%%%%%%%%%%%%%%%% Check for optional inputs %%%%%%%%%%%%%%%%%%%%%%
for i = 1:(numel(varargin)-1)  
    if i ~= numel(varargin)        
        eval([varargin{i} '=varargin{i+1};']);                
    end    
end

%% %%%%%%%%%%%%%%%%%%%%%%% Save key sampling parameters %%%%%%%%%%%%%%%%%%%
proteinSamplingInfo = struct;
proteinSamplingInfo.ROIRadiusSpot = ROIRadiusSpot_um;
proteinSamplingInfo.minSampleSep_um = minSampleSep_um;
proteinSamplingInfo.minEdgeSep_um = minEdgeSep_um; 
proteinSamplingInfo.xy_sigma_nuclear_um = xy_sigma_nuclear_um;
proteinSamplingInfo.ROIRadiusSpot_um = ROIRadiusSpot_um;
proteinSamplingInfo.snippet_size_um = snippet_size_um;
proteinSamplingInfo.xy_sigma_um = xy_sigma_um;
proteinSamplingInfo.z_sigma_um = z_sigma_um;
proteinSamplingInfo.min_nucleus_radius_um = min_nucleus_radius_um;
proteinSamplingInfo.max_nucleus_radius_um = max_nucleus_radius_um;

%% %%%%%%%%%%%%%%%%%%%%%%% Get project info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[liveProject, ~, nucleusName, hasAPInfo, has3DSpotInfo] = headerFunction(projectName);
use3DSpotInfo = use3DSpotInfo&&has3DSpotInfo;
proteinSamplingInfo.use3DSpotInfo = use3DSpotInfo;

%% %%%%%%%%%%%%%%%%%%%%%%% Load data and clean trace data %%%%%%%%%%%%%%%%%
load(nucleusName,'nucleus_struct')

if use_psf_fit_dims
    load([liveProject.dataPath '/psf_dims.mat'],'psf_dims')
end

% make paths 
snipPath = [liveProject.dataPath '/qc_images/'];
refPath = [liveProject.dataPath '/refFrames/'];
mkdir(refPath)
mkdir(snipPath)

% remove frames where no particle was observed
spot_struct_protein = truncateParticleFields(nucleus_struct,use3DSpotInfo,hasAPInfo);

%% %%%%%%%%%%%%%%%%%%%%%%% Generate indexing vectors  %%%%%%%%%%%%%%%%%%%%%
refVecStruct = generateReferenceVectors(spot_struct_protein,refPath,use3DSpotInfo,ignoreQC);

qc_structure = struct;


%% %%%%%%%%%%%%%%%%%%%%%%% Initialize enrichment-related fields  %%%%%%%%%%

[spot_struct_protein, new_vec_fields] = initializeProteinFields(spot_struct_protein, use3DSpotInfo);

new_snip_fields = {'spot_protein_snips', 'edge_null_protein_snips',...
    'spot_mcp_snips','edge_null_mcp_snips'};


%% %%%%%%%%%%%%%%%%%%%%%%% Nucleus segmentation  %%%%%%%%%%%%%%%%%%%%%%%%%%

% first check to see if segmentation files exist
[segmentNuclei, segmentIndices] = ...
          handleSegmentationOptions(refVecStruct,segmentNuclei);

if segmentNuclei
    disp('segmenting nuclei...')    
    nuclearSegmentation(liveProject, refVecStruct, segmentIndices, NumWorkers);      
end

%% %%%%%%%%%%%%%%%%%%%%%%% Local Protein Sampling %%%%%%%%%%%%%%%%%%%%%%%%%

disp('taking protein samples...')
for i = 1:size(set_frame_array,1)    
    tic
    currentSetID = set_frame_array(i,1);
    currentFrame = set_frame_array(i,2);  
    
    Prefix = liveProject.includedExperimentNames{currentSetID}; 
    currExperiment = liveProject.includedExperiments{currentSetID};                  
    proteinChannel = currExperiment.inputChannels;
    mcpChannel = currExperiment.spotChannels;
    
    %% %%%%%%%%%%%%%%%%%%%%% Set size parameters  %%%%%%%%%%%%%%%%%%%%%%%%%
    currExperiment = LiveExperiment(Prefix);
    PixelSize = currExperiment.pixelSize_nm / 1e3; 
    zStep = currExperiment.zStep_um;
    
    if length(mcpChannel) > 1 || length(mcpChannel) > 1
      error('This pipeline does not currently support multiple spot channels')
    elseif length(proteinChannel) > 1
      error('This pipeline does not currently support multiple protein channels')
    end
    
    % Generate reference vectors
    xDim = currExperiment.xDim;
    yDim = currExperiment.yDim;
    zDim = currExperiment.zDim;
    [x_ref,y_ref,z_ref] = meshgrid(1:xDim,1:yDim,1:zDim);
    
    % calculate basic parameters for sampling
    sampParamNames = fieldnames(proteinSamplingInfo);
    tempParams = struct;
    for s = 1:length(sampParamNames)
      paramName = sampParamNames{s};
      tempParams.(paramName(1:end-3)) = proteinSamplingInfo.(paramName) / PixelSize;
    end
    
    min_nucleus_area = pi*tempParams.min_nucleus_radius^2;
    max_nucleus_area = pi*tempParams.max_nucleus_radius^2;
    
    % calculate characteristic drift to use for simulated spot
    driftTol = calculateVirtualSpotDrift(refVecStruct,PixelSize);
    
    % load spot and nucleus reference frames
    nc_ref_name = [refPath 'nc_ref_frame_set' sprintf('%02d',currentSetID) '_frame' sprintf('%03d',currentFrame) '.mat'];
    load(nc_ref_name,'nc_ref_frame');
    nc_dist_frame = bwdist(~nc_ref_frame);    
    spot_ref_name = [refPath 'spot_roi_frame_set' sprintf('%02d',currentSetID) '_frame' sprintf('%03d',currentFrame) '.mat'];
    load(spot_ref_name,'spot_dist_frame');
%     spot_roi_frame = bwlabel(spot_dist_frame <= roi_rad_spot_pix);
    
    % get nucleus
    frame_set_filter = refVecStruct.setID_ref==currentSetID&refVecStruct.frame_ref==currentFrame;
    frame_set_indices = find(frame_set_filter);                

    % load stacks    
    proteinPath = [currExperiment.preFolder  Prefix '_' sprintf('%03d',currentFrame) '_ch0' num2str(proteinChannel) '.tif'];
    protein_stack = imreadStack2(proteinPath, yDim, xDim, zDim+2);  
    protein_stack = protein_stack(:,:,2:end); %NL: need to make this dynamic
    
    mcpPath = [currExperiment.preFolder  Prefix '_' sprintf('%03d',currentFrame) '_ch0' num2str(mcpChannel) '.tif'];
    mcp_stack = imreadStack2(mcpPath, yDim, xDim, zDim+2);   
    mcp_stack = mcp_stack(:,:,2:end);
    
    % generate lookup table of inter-nucleus distances
    nc_x_vec = round(refVecStruct.nc_x_ref(frame_set_filter));
    nc_y_vec = round(refVecStruct.nc_j_ref(frame_set_filter));  
    x_dist_mat = repmat(nc_x_vec,numel(nc_x_vec),1)-repmat(nc_x_vec',1,numel(nc_x_vec));
    y_dist_mat = repmat(nc_y_vec,numel(nc_y_vec),1)-repmat(nc_y_vec',1,numel(nc_y_vec));
    r_dist_mat = sqrt(double(x_dist_mat).^2 + double(y_dist_mat).^2);            
    
    % initialize arrays to store relevant info 
%     for j = 1:numel(new_snip_fields)
%         eval([new_snip_fields{j} ' = NaN(2*pt_snippet_size+1,2*pt_snippet_size+1,numel(spot_x_vec));']);
%     end    
    
    % iterate through spots
    qc_mat = struct;
    j_pass = 1;
    for j = frame_set_indices%1:numel(nc_x_vec)
        % get indexing info         
        spotIndex = refVecStruct.particle_index_ref(j);
        spotSubIndex = refVecStruct.particle_subindex_ref(j);
        
        % get location info
        x_nucleus = round(refVecStruct.nc_x_ref(j));
        y_nucleus = round(refVecStruct.nc_j_ref(j));          
            
        x_spot = refVecStruct.spot_x_ref(j);
        x_index = min([ xDim max([1 round(x_spot)])]);
        y_spot = refVecStruct.spot_y_ref(j);
        y_index = min([ yDim max([1 round(y_spot)])]);
        z_spot = refVecStruct.spot_z_ref(j)-1.0; % adjust for z padding
        
        if ~refVecStruct.particle_qcFlag_ref(j)
            continue
        end
                
        % extract mask 
        spot_nc_mask = nc_ref_frame == refVecStruct.master_nucleusID_ref(j);             
            
        %%%%%%%%%%%%%%%%%%%%%% Sample protein levels %%%%%%%%%%%%%%%%%%%%%%

        % sample protein near locus
        spot_struct_protein(spotIndex).spot_protein_vec(spotSubIndex) = sample_protein_3D(...
              x_spot,y_spot,z_spot,x_ref,y_ref,z_ref,tempParams.xy_sigma,tempParams.z_sigma,protein_stack);
        spot_struct_protein(spotIndex).spot_mcp_vec(spotSubIndex) = sample_protein_3D(...
              x_spot,y_spot,z_spot,x_ref,y_ref,z_ref,tempParams.xy_sigma,tempParams.z_sigma,mcp_stack); 
        
        % make sure size is reasonable and that spot is inside nucleus
        if sum(spot_nc_mask(:)) < min_nucleus_area || sum(spot_nc_mask(:)) > max_nucleus_area || ~spot_nc_mask(y_index,x_index) %|| int_it==0  
            spot_struct_protein(spotIndex).edge_qc_flag_vec(spotSubIndex) = -1;            
            spot_struct_protein(spotIndex).serial_qc_flag_vec(spotSubIndex) = -1;
            continue
        end 
                      
        % sample snippets
        spot_protein_snips(:,:,j) = sample_snip(x_spot,y_spot,snippet_size,protein_slice,spot_nc_mask);
        spot_mcp_snips(:,:,j) = sample_snip(x_spot,y_spot,snippet_size,mcp_slice,spot_nc_mask);                 

        % Take average across all pixels within 1.5um of nuclues center           
        spot_struct_protein(spotIndex).nucleus_protein_vec(spotSubIndex) = sample_protein_3D(x_nucleus,y_nucleus,z_spot,...
            x_ref,y_ref,z_ref,tempParams.xy_sigma_nuclear,tempParams.z_sigma,protein_stack);    
        
          
        %% %%%%%%%%%%%% Draw edge control spot %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [spot_struct_protein, null_mask] = findEdgeControlWrapper(spot_struct_protein,refVecStruct,x_index,y_index,spot_nc_mask,...
          spotIndex,spotSubIndex,r_dist_mat,frame_set_indices,j_pass);
        
        % Draw control samples (as appropriate)    
        if spot_struct_protein(spotIndex).edge_qc_flag_vec(spotSubIndex) > 0  
            edge_control_x = spot_struct_protein(spotIndex).edge_null_x_vec(spotSubIndex);
            edge_control_y = spot_struct_protein(spotIndex).edge_null_y_vec(spotSubIndex);     
            
            spot_struct_protein(spotIndex).edge_null_protein_vec(spotSubIndex) = sample_protein_3D(...
              edge_control_x,edge_control_y,z_spot,x_ref,y_ref,z_ref,tempParams.xy_sigma,tempParams.z_sigma,protein_stack);
            spot_struct_protein(spotIndex).edge_null_mcp_vec(spotSubIndex) = sample_protein_3D(...
              edge_control_x,edge_control_y,z_spot,x_ref,y_ref,z_ref,tempParams.xy_sigma,tempParams.z_sigma,protein_stack);
         
            % draw snips    
            edge_null_protein_snips(:,:,j) = sample_snip(edge_control_x,edge_control_y,snippet_size,protein_slice,nc_ref_frame>0);
            edge_null_mcp_snips(:,:,j) = sample_snip(edge_control_x,edge_control_y,snippet_size,mcp_slice,nc_ref_frame>0);                                        
        end                  
        
        %% %%%%%%%%%%%% Draw serialized control spot %%%%%%%%%%%%%%%%%%%%%%%           
        [serial_control_x, serial_control_y, serial_edge_dist] = ...
          drawSerializedControlSpot(spot_struct_protein,spotIndex,spotSubIndex,...
                   spot_sep_vec, tempParams, nc_edge_dist_vec);
        
        % sample protein 
        spot_struct_protein(spotIndex).serial_null_protein_vec(spotSubIndex) = sample_protein_3D(...
              serial_control_x,serial_control_y,z_spot,x_ref,y_ref,z_ref,tempParams.xy_sigma,tempParams.z_sigma,protein_stack);
        spot_struct_protein(spotIndex).serial_null_mcp_vec(spotSubIndex) = sample_protein_3D(...
              serial_control_x,serial_control_y,z_spot,x_ref,y_ref,z_ref,tempParams.xy_sigma,tempParams.z_sigma,mcp_stack);
        % record
        spot_struct_protein(spotIndex).serial_null_x_vec(spotSubIndex) = serial_control_x;
        spot_struct_protein(spotIndex).serial_null_y_vec(spotSubIndex) = serial_control_y;
        spot_struct_protein(spotIndex).serial_null_edge_dist_vec(spotSubIndex) = serial_edge_dist;
        
        %% %%%%%%%%%%%% Check for sister spot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        x_spot_sister = NaN;
        y_spot_sister = NaN;
        nucleus_indices = find(refVecStruct.master_nucleusID_ref(j)==refVecStruct.master_nucleusID_ref(frame_set_filter));
        if length(nucleus_indices) == 2            
            sister_index = nucleus_indices(nucleus_indices~=j_pass);
            x_spot_sister = spot_x_vec(sister_index);
            y_spot_sister = spot_y_vec(sister_index);
        end            
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        % save qc data                 
        qc_mat(numel(nc_x_vec)-j+1).setID = currentSetID; 
        qc_mat(numel(nc_x_vec)-j+1).frame = currentFrame;
        qc_mat(numel(nc_x_vec)-j+1).nc_index = nc_lin_index_vec(j);
        qc_mat(numel(nc_x_vec)-j+1).nc_sub_index = nc_sub_index_vec(j);
        qc_mat(numel(nc_x_vec)-j+1).qcFlag = edge_qc_flag_vec(j);            
        qc_mat(numel(nc_x_vec)-j+1).xp = x_spot;
        qc_mat(numel(nc_x_vec)-j+1).yp = y_spot;  
        qc_mat(numel(nc_x_vec)-j+1).xp_sister = x_spot_sister;
        qc_mat(numel(nc_x_vec)-j+1).yp_sister = y_spot_sister;  
        qc_mat(numel(nc_x_vec)-j+1).xc_edge = edge_null_x_vec(j);
        qc_mat(numel(nc_x_vec)-j+1).yc_edge = edge_null_y_vec(j);
        qc_mat(numel(nc_x_vec)-j+1).xc_serial = serial_null_x_vec(j);
        qc_mat(numel(nc_x_vec)-j+1).yc_serial = serial_null_y_vec(j);       
        qc_mat(numel(nc_x_vec)-j+1).particleID = particle_id_vec(j);
        qc_mat(numel(nc_x_vec)-j+1).serial_qc_flag = serial_qc_flag_vec(j);
        qc_mat(numel(nc_x_vec)-j+1).edge_qc_flag = edge_qc_flag_vec(j);        
        sz = nb_size;
        edge_dist_mat = nc_dist_frame;
        edge_dist_mat(~spot_nc_mask&~null_mask) = 0;
        if edge_qc_flag_vec(j) == 2         
            sz = max([nb_size,abs(x_nucleus - edge_null_x_vec(j)),abs(y_nucleus - edge_null_y_vec(j))...
                abs(x_nucleus - serial_null_x_vec(j)),abs(y_nucleus - serial_null_y_vec(j))]);            
        end
        y_range = max(1,y_nucleus-sz):min(yDim,y_nucleus+sz);
        x_range = max(1,x_nucleus-sz):min(xDim,x_nucleus+sz);
        qc_mat(numel(nc_x_vec)-j+1).x_origin = x_range(1);
        qc_mat(numel(nc_x_vec)-j+1).y_origin = y_range(1);
        qc_mat(numel(nc_x_vec)-j+1).mcp_snip = mcp_slice(y_range,x_range);
        qc_mat(numel(nc_x_vec)-j+1).protein_snip = protein_slice(y_range,x_range);
        qc_mat(numel(nc_x_vec)-j+1).edge_dist_snip = edge_dist_mat(y_range,x_range);   
        
        % increment
        j_pass = j_pass + 1;
    end 
    qc_structure(i).qc_mat = fliplr(qc_mat);  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save snip data   
    % initialize struct to store snip data
    snip_data = struct;        
    
    % store key ID variables
    snip_data.frame = currentFrame;
    snip_data.setID = currentSetID;
    snip_data.particle_id_vec = particle_id_vec;
    snip_data.spot_edge_dist_vec = spot_edge_dist_vec;
    % indexing vectors    
    snip_data.nc_sub_index_vec = nc_sub_index_vec; 
    snip_data.nc_lin_index_vec = nc_lin_index_vec; 
    snip_data.nc_master_vec = nc_master_vec;    
    % specify name
     % read snip file    
    snip_name = ['snip_data_F' sprintf('%03d',currentFrame) '_S' sprintf('%02d',currentSetID)]; 
    if write_snip_flag            
        blank = struct;
        save([dataPath 'snip_data.mat'],'blank','-v7.3')    
        write_snip_flag = true;
    end
    snip_file = matfile([dataPath 'snip_data.mat'],'Writable',true);    
    snip_file.(snip_name)= snip_data;        
    clear snip_file;    
    % report time
    t = round(toc);
    disp([num2str(i) ' of ' num2str(size(set_frame_array,1)) ' frames completed (' num2str(t) ' sec)'])         
end
disp('saving qc frames...')
% save qc data
tic
particle_index = unique([spot_struct_protein.particleID]);
particle_index = particle_index(~isnan(particle_index));
qc_particles = randsample(particle_index,min([100,numel(particle_index)]),false);
particle_index_full = [];
particle_frames_full = [];
for i = 1:numel(qc_structure)
    qc_mat = qc_structure(i).qc_mat;
    for  j = 1:numel(qc_mat)
        qc_spot = qc_mat(j);
        if ~isfield(qc_spot,'ParticleID')
            continue
        end
        ParticleID = qc_spot.particleID;
        if isempty(ParticleID) || ~ismember(ParticleID,qc_particles)
            continue
        end        
        currentFrame = qc_spot.frame;      
        particle_index_full = [particle_index_full ParticleID];
        particle_frames_full = [particle_frames_full currentFrame];        
        save_name = [snipPath 'pt' num2str(1e4*ParticleID) '_frame' sprintf('%03d',currentFrame) '.mat'];
        save(save_name,'qc_spot');
    end
end
[particle_index_full, si] = sort(particle_index_full);
particle_frames_full = particle_frames_full(si);

qc_ref_struct.particle_frames_full = particle_frames_full;
qc_ref_struct.particle_index_full = particle_index_full;
toc
% save updated nucleus structure
disp('saving nucleus structure...')
spot_struct_protein = spot_struct_protein;
save([dataPath 'qc_ref_struct.mat'],'qc_ref_struct')
save([dataPath 'spot_struct_protein.mat'],'spot_struct_protein','-v7.3') 