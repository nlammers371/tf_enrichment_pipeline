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
close all force

ROIRadiusSpot_um = .2; % radus (um) of region used to query and compare TF concentrations
minSampleSep_um = 1.5; %um
minEdgeSep_um = .25; %um
snippet_size_um = 1.5;
min_nucleus_radius_um = 2;
max_nucleus_radius_um = 4;

segmentNuclei = 0;
use3DSpotInfo = 1;
NumWorkers = 24;
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
[refVecStruct, set_frame_array] = generateReferenceVectors(spot_struct_protein,refPath,use3DSpotInfo,ignoreQC);

qc_structure = struct;


%% %%%%%%%%%%%%%%%%%%%%%%% Initialize enrichment-related fields  %%%%%%%%%%

[spot_struct_protein, ~] = initializeProteinFields(spot_struct_protein, use3DSpotInfo);

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

h = waitbar(0,'Sampling local protein...');
for i = 1:size(set_frame_array,1)    
    waitbar(i/size(set_frame_array,1),h)
    
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
    tempParams.z_sigma = proteinSamplingInfo.z_sigma_um / zStep;
    min_nucleus_area = pi*tempParams.min_nucleus_radius^2;
    max_nucleus_area = pi*tempParams.max_nucleus_radius^2;
    tempParams.snippet_size = round(tempParams.snippet_size);
    
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
    nc_y_vec = round(refVecStruct.nc_y_ref(frame_set_filter));  
    x_dist_mat = repmat(nc_x_vec,numel(nc_x_vec),1)-repmat(nc_x_vec',1,numel(nc_x_vec));
    y_dist_mat = repmat(nc_y_vec,numel(nc_y_vec),1)-repmat(nc_y_vec',1,numel(nc_y_vec));
    r_dist_mat = sqrt(double(x_dist_mat).^2 + double(y_dist_mat).^2);            
    
    % initialize temporary arrays to store snip info 
    for j = 1:numel(new_snip_fields)
        eval([new_snip_fields{j} ' = NaN(2*tempParams.snippet_size+1,2*tempParams.snippet_size+1,numel(nc_x_vec));']);
    end    
    
    % iterate through spots
    qc_mat = struct;
    j_pass = 1;
    for j = frame_set_indices%1:numel(nc_x_vec)
        % get indexing info         
        spotIndex = refVecStruct.particle_index_ref(j);
        spotSubIndex = refVecStruct.particle_subindex_ref(j);
        
        % get location info
        x_nucleus = round(refVecStruct.nc_x_ref(j));
        y_nucleus = round(refVecStruct.nc_y_ref(j));          
            
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
            
        %% %%%%%%%%%%%%%%%%%%% Sample protein levels %%%%%%%%%%%%%%%%%%%%%%

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
        spot_protein_snips(:,:,j_pass) = sample_snip_3D(x_spot,y_spot,z_spot,z_ref,tempParams.snippet_size,tempParams.z_sigma,protein_stack);
        spot_mcp_snips(:,:,j_pass) = sample_snip_3D(x_spot,y_spot,z_spot,z_ref,tempParams.snippet_size,tempParams.z_sigma,mcp_stack); 

        % Take average across all pixels within 1.5um of nuclues center           
        spot_struct_protein(spotIndex).nucleus_protein_vec(spotSubIndex) = sample_protein_3D(x_nucleus,y_nucleus,z_spot,...
            x_ref,y_ref,z_ref,tempParams.xy_sigma_nuclear,tempParams.z_sigma,protein_stack);    
        
          
        %% %%%%%%%%%%%% Draw edge control spot %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        spot_struct_protein = findEdgeControlWrapper(spot_struct_protein,refVecStruct,x_index,y_index,spot_nc_mask,...
                                                      spotIndex,spotSubIndex,r_dist_mat,frame_set_indices,nc_dist_frame,...
                                                      spot_dist_frame,x_ref,y_ref,tempParams,j_pass);
        
        % Draw control samples (as appropriate)    
        if spot_struct_protein(spotIndex).edge_qc_flag_vec(spotSubIndex) > 0  
            edge_control_x = spot_struct_protein(spotIndex).edge_null_x_vec(spotSubIndex);
            edge_control_y = spot_struct_protein(spotIndex).edge_null_y_vec(spotSubIndex);     
            
            spot_struct_protein(spotIndex).edge_null_protein_vec(spotSubIndex) = sample_protein_3D(...
              edge_control_x,edge_control_y,z_spot,x_ref,y_ref,z_ref,tempParams.xy_sigma,tempParams.z_sigma,protein_stack);
            spot_struct_protein(spotIndex).edge_null_mcp_vec(spotSubIndex) = sample_protein_3D(...
              edge_control_x,edge_control_y,z_spot,x_ref,y_ref,z_ref,tempParams.xy_sigma,tempParams.z_sigma,protein_stack);
         
            % draw snips    
            edge_null_protein_snips(:,:,j_pass) = sample_snip_3D(edge_control_x,edge_control_y,z_spot,z_ref,tempParams.snippet_size,tempParams.z_sigma,protein_stack);
            edge_null_mcp_snips(:,:,j_pass) = sample_snip_3D(edge_control_x,edge_control_y,z_spot,z_ref,tempParams.snippet_size,tempParams.z_sigma,mcp_stack);
        end                  
        
        %% %%%%%%%%%%%% Draw serialized control spot %%%%%%%%%%%%%%%%%%%%%%%           
        [serial_control_x, serial_control_y, serial_edge_dist] = ...
          drawSerializedControlSpot(spot_struct_protein,spotIndex,spotSubIndex,...
                                  tempParams,nc_dist_frame,spot_dist_frame,spot_nc_mask,driftTol);
        
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
%         x_spot_sister = NaN;
%         y_spot_sister = NaN;
%         nucleus_indices = find(refVecStruct.master_nucleusID_ref(j)==refVecStruct.master_nucleusID_ref(frame_set_filter));
%         if length(nucleus_indices) == 2            
%             sister_index = nucleus_indices(nucleus_indices~=j_pass);
%             x_spot_sister = spot_x_vec(sister_index);
%             y_spot_sister = spot_y_vec(sister_index);
%         end            
               
        %% %%%%%%%%%%%%% save qc data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 
        qc_mat = updateQCmat;
        
        % increment
        j_pass = j_pass + 1;
    end 
    qc_structure(i).qc_mat = fliplr(qc_mat);  
        
    %% %%%%%%%%%%%%%%%%%%%% save snip data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    % initialize struct to store snip data
    snip_data = struct;        
    
    % store key ID variables
    snip_data.frame = currentFrame;
    snip_data.setID = currentSetID;
    snip_data.particle_id_vec = refVecStruct.particleID_ref(frame_set_filter);
%     snip_data.spot_edge_dist_vec = spot_edge_dist_vec;
    % indexing vectors    
    snip_data.nc_sub_index_vec = refVecStruct.particle_subindex_ref(frame_set_filter);
    snip_data.nc_lin_index_vec = refVecStruct.particle_index_ref(frame_set_filter); 
    snip_data.nc_master_vec = refVecStruct.master_nucleusID_ref(frame_set_filter);    
    
    % specify name      
    snip_name = ['snip_data_F' sprintf('%03d',currentFrame) '_S' sprintf('%02d',currentSetID)]; 
%     if write_snip_flag            
%         blank = struct;
%         save([dataPath 'snip_data.mat'],'blank','-v7.3')    
%         write_snip_flag = false;
%     end
    snip_file = matfile([liveProject.dataPath 'snip_data.mat'],'Writable',true);    
    snip_file.(snip_name)= snip_data;        
    clear snip_file;    
        
end


% disp('saving qc frames...')
% % save qc data
% tic
% particle_index = unique([spot_struct_protein.particleID]);
% particle_index = particle_index(~isnan(particle_index));
% qc_particles = randsample(particle_index,min([100,numel(particle_index)]),false);
% particle_index_full = [];
% particle_frames_full = [];
% for i = 1:numel(qc_structure)
%     qc_mat = qc_structure(i).qc_mat;
%     for  j = 1:numel(qc_mat)
%         qc_spot = qc_mat(j);
%         if ~isfield(qc_spot,'ParticleID')
%             continue
%         end
%         ParticleID = qc_spot.particleID;
%         if isempty(ParticleID) || ~ismember(ParticleID,qc_particles)
%             continue
%         end        
%         currentFrame = qc_spot.frame;      
%         particle_index_full = [particle_index_full ParticleID];
%         particle_frames_full = [particle_frames_full currentFrame];        
%         save_name = [snipPath 'pt' num2str(1e4*ParticleID) '_frame' sprintf('%03d',currentFrame) '.mat'];
%         save(save_name,'qc_spot');
%     end
% end
% [particle_index_full, si] = sort(particle_index_full);
% particle_frames_full = particle_frames_full(si);
% 
% qc_ref_struct.particle_frames_full = particle_frames_full;
% qc_ref_struct.particle_index_full = particle_index_full;
% toc
% save updated nucleus structure
disp('saving nucleus structure...')

% save([dataPath 'qc_ref_struct.mat'],'qc_ref_struct')
save([liveProject.dataPath 'spot_struct_protein.mat'],'spot_struct_protein','-v7.3') 
save([liveProject.dataPath 'proteinSamplingInfo.mat'],'proteinSamplingInfo');