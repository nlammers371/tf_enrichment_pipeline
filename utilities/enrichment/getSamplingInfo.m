function samplingInfo = getSamplingInfo(samplingInfo,liveProject,...
              proteinSamplingInfo,RefStruct,NewSetFlag)                

    % get path to reference frames   
    refPath = [liveProject.dataPath '/refFrames/'];        
    currExperiment = liveProject.includedExperiments{samplingInfo.SetID}; 
    Prefix = liveProject.includedExperimentNames{samplingInfo.SetID}; 
    
             
    samplingInfo.proteinChannel = currExperiment.inputChannels;
    samplingInfo.mcpChannel = currExperiment.spotChannels;        
    
    % generate path info for protein stacks
    samplingInfo.proteinPath = [currExperiment.preFolder  Prefix '_' sprintf('%03d',samplingInfo.Frame) '_ch0' num2str(samplingInfo.proteinChannel) '.tif'];
    samplingInfo.mcpPath = [currExperiment.preFolder  Prefix '_' sprintf('%03d',samplingInfo.Frame) '_ch0' num2str(samplingInfo.mcpChannel) '.tif'];
    
    %% %%%%%%%%%%%%%%%%%%%%% Set size parameters  %%%%%%%%%%%%%%%%%%%%%%%%%
    if NewSetFlag      
      PixelSize = currExperiment.pixelSize_nm / 1e3;  
      zStep = currExperiment.zStep_um;    

      if length(samplingInfo.mcpChannel) > 1 
        error('This pipeline does not currently support multiple spot channels')
      elseif length(samplingInfo.proteinChannel) > 1
        error('This pipeline does not currently support multiple protein channels')
      end

      % Generate reference vectors
      samplingInfo.xDim = currExperiment.xDim;
      samplingInfo.yDim = currExperiment.yDim;
      samplingInfo.zDim = currExperiment.zDim;
      [samplingInfo.x_ref,samplingInfo.y_ref,samplingInfo.z_ref] = meshgrid(1:samplingInfo.xDim,1:samplingInfo.yDim,1:samplingInfo.zDim);
    
      % calculate basic parameters for sampling
      sampParamNames = fieldnames(proteinSamplingInfo);   
      for s = 1:length(sampParamNames)
        paramName = sampParamNames{s};
        samplingInfo.(paramName(1:end-3)) = proteinSamplingInfo.(paramName) / PixelSize;
      end
      samplingInfo.z_sigma = proteinSamplingInfo.z_sigma_um / zStep;
      samplingInfo.min_nucleus_area = pi*samplingInfo.min_nucleus_radius^2;
      samplingInfo.max_nucleus_area = pi*samplingInfo.max_nucleus_radius^2;
      samplingInfo.snippet_size = round(samplingInfo.snippet_size);

      % calculate characteristic drift to use for simulated spot
      samplingInfo.driftTol = calculateVirtualSpotDrift(RefStruct,samplingInfo.SetID);
    end
    
    % load spot and nucleus reference frames     
    nc_ref_name = [refPath 'nc_ref_frame_set' sprintf('%02d',samplingInfo.SetID) '_frame' sprintf('%03d',samplingInfo.Frame) '.mat'];
    temp = load(nc_ref_name,'nc_ref_frame');
    samplingInfo.nc_ref_frame_raw = temp.nc_ref_frame;            

    spot_ref_name = [refPath 'spot_roi_frame_set' sprintf('%02d',samplingInfo.SetID) '_frame' sprintf('%03d',samplingInfo.Frame) '.mat'];
    temp = load(spot_ref_name,'spot_dist_frame');    
    spot_dist_frame = temp.spot_dist_frame;
    samplingInfo.spot_dist_frame = spot_dist_frame;

    % get indices of particles in current set/frame 
    samplingInfo.FrameSetFilter = RefStruct.setID_ref==samplingInfo.SetID&RefStruct.frame_ref==samplingInfo.Frame;
    samplingInfo.frame_set_indices = find(samplingInfo.FrameSetFilter);