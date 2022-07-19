function spot_struct = main01_compile_traces(projectName,varargin)

%
% nucleusStruct = main01_compile_traces(dataStatusTab,dropboxFolder,varargin)
%
% DESCRIPTION
% Function to compile relevant outputs from main mRNADynamics image
% analysis pipeline across multiple experiments
%
%
% ARGUMENTS
% projectName: master ID variable (should match a tab name in the
%              DataStatus.xlsx spreadsheet)
%
% OPTIONS
% firstNC: script defaults to taking only nc14 traces. If you want
%          earlier traces, pass 'firstNC', followed by desired nuclear
%          cycle number
% NC: allows specification of singleNC
%
% other: script allows any default variable to be set using format:
%       "VariableNameString", VariableValue
%
% dt: specifies the time interval in seconds to be used for interpolated traces in spot_struct.
%   This allows for downsampling of traces. Added by GM 3/9/21.
%
% OUTPUT
% spot_struct: compiled data set contain key nucleus and
%                 particle attributes, contains the following fields:
%
%     setID: ID number for this dataset (experiment), see set_key.mat to
%            match setID to experiment Prefix
%     sourcePath: full path to the experiment's DynamicsResults folder
%     APFlag: indicates whether or not there is AP position infor for this
%             experiment
%     threeDFlag: indicates whether or not there is 3D spot fitting data
%                 for this experiment
%     xDim: # of pixels in X
%     yDim: # of pixels in Y
%     zDim: # of steps in Z
%     zStep: size of steps in z, in um
%     pixelSize: size of xy pixels, in um
%     particleID:
%     xPosParticle:
%     yPosParticle:
%     zPosParticle:
%     APPosParticle:
%     spotFrames:
%     xPosParticle3D:
%     yPosParticle3D:
%     zPosParticle3D:
%     fluo3D:
%     fluo:
%     fluoOffset:
%     qcFlag:
%     N: ???
%     sparsity:
%     xPosNucleus: x position of center of nucleus
%     yPosNucleus: y position of center of nucleus
%     rawNCProtein:
%     frames:
%     nucleusID
%     ncID
%     ncStart
%     minDP
%     minTime
%     pctSparsity
%
%     time
%
%     sisterIndex
%     sisterParticleID

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Set defaults %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all force
addpath(genpath('utilities'));
% Defaults
firstNC = 14;   % first nuclear cycle to pull data from
lastNC = 14;
minDP = 14;     % particles with fewer than minDP points will be flagged
pctSparsity = 50;   %
twoSpotFlag = contains(projectName, '2spot');
minTime = 0*60; % take no fluorescence data prior to this point
tresInterp = 20;
tresInterpFloor = 20;
calculatePSF = false;
% projectName = projectName;
SpotChannelIndex = 1; % this does nothing at the moment, but can serve as a starting point if ever we wish to analyze two-spot-two-color data
ncRefVec = 9:14;

%% %%%%%%%%%%%%%%%%%%%%%%% Process input parameters %%%%%%%%%%%%%%%%%%%%%%%

for i = 1:numel(varargin)
    if ischar(varargin{i})
        if strcmpi(varargin{i}, 'firstNC')
            if  i < numel(varargin) && isnumeric(varargin{i+1}) && (9 <= varargin{i+1} && varargin{i+1} <= 14)
                firstNC = varargin{i+1};
            else
                error(['The ''firstNC'' option must be followed by an integer from 9 to 14, inclusive. Default firstNC is ' num2str(firstNC) '.'])
            end
        elseif strcmpi(varargin{i}, 'NC')
            if  i < numel(varargin) && isnumeric(varargin{i+1}) && (9 <= varargin{i+1} && varargin{i+1} <= 14)
                lastNC = varargin{i+1};
                firstNC = varargin{i+1};
            else
                error(['The ''NC'' option must be followed by an integer from 9 to 14, inclusive. Default NC is ' num2str(firstNC) '.'])
            end
        elseif strcmpi(varargin{i}, 'dt')
                dt = varargin{i+1};
   
        elseif strcmpi(varargin{i}, 'projectName')
            error('Pipeline currently does not support alternate names for projects (otehr than what is on the data status tab)')
        elseif i < numel(varargin)
            eval([varargin{i} '=varargin{i+1};']);
        end
    end
end

%% %%%%%%%%%%%%%%%%%%%%% Set Path Specs, ID Vars %%%%%%%%%%%%%%%%%%%%%%%%%%

[liveProject, numExperiments, dataName, hasAPInfo, has3DSpotInfo, hasProteinInfo] = headerFunction(projectName);

%% %%%%%%%%%%%%%%%%% Extract relevant processed data %%%%%%%%%%%%%%%%%%%%%%

% Generate master structure with info on all nuclei and traces in
% constituent sets
spot_struct = [];

% Add data from each experiment to the master structure
h = waitbar(0,'Compiling data ...');
nc_ap_flags = false(1,numExperiments);

for i = 1:numExperiments
    
    waitbar((i-1)/numExperiments,h, ['Compiling data for dataset ' num2str(i) ' of ' num2str(numExperiments)])
    
    setID = i;
    currExperiment = liveProject.includedExperiments{i};
    
    %%%%%%%% Read in processed data from main mRNADyanmics pipeline %%%%%%%
    
    schnitzcells = getSchnitzcells(currExperiment);
    processedData = getCompiledParticles(currExperiment);   % this structure contains all info compiled for this experiment
    compiledParticles = processedData.CompiledParticles;    % this is the cell array containing the CompiledParticles cell array
    if iscell(compiledParticles)
        compiledParticles = compiledParticles{SpotChannelIndex};               % this assumes there is only one channel with spots to analyze
    end
    
    % MT 2020-08-11: Right now we don't use any of the data in
    % CompiledNuclei (from mRNADynamics pipeline) because we generate all
    % of it ourselves in main02. This might change in the future, so I'm
    % including this section so it's easy to use should we need it
    %     compiledNuclei = processedData.CompiledNuclei;        % this is the CompiledNuclei structure we actually want to use
    %     if iscell(compiledNuclei)
    %         compiledNuclei = compiledNuclei{1};                     % unlike compiledParticles, I don't think this gets stored in a cell array, but just in case...
    %     end
    
    
    %%%%%%%%%%%%%%%%%%% Get basic data characteristics %%%%%%%%%%%%%%%%%%%%
    
    % Check if certain info is present for this experiment
    DetrendedZFlag = isfield(compiledParticles, 'zPosDetrended');   %Z pos that accounts for stack repositioning
    
    % Pull trace, time, and frame variables
    timeRaw = processedData.ElapsedTime*60; %time vector, [sec]
    framesRaw = 1:length(timeRaw); % Frame list
    tracesRaw = processedData.AllTracesVector;  %array with a column for each trace
    if iscell(tracesRaw)        %check to see if traces are stored in cell array
        tracesRaw = tracesRaw{SpotChannelIndex};
    end
    
    % Find start frames for nc's we wish to include. Pull this from
    % FrameInfo if not in movie database
    anaphaseFrames = currExperiment.anaphaseFrames;
    if all(isnan(anaphaseFrames))
        FrameInfo = getFrameInfo(currExperiment);
        ncVec = [FrameInfo.nc];
        ncIndex = unique(ncVec);
        for n = 1:length(ncIndex)
            anaphaseFrames(ncRefVec==ncIndex(n)) = find(ncVec==ncIndex(n),1);
        end
    end
    ncIndices = find(ncRefVec>=firstNC&ncRefVec<=lastNC&~isnan(anaphaseFrames'));
    ncStartFrameVec = anaphaseFrames(ncIndices)';
    ncStartFrameVec(end) = max([1 ncStartFrameVec(end)]);    
    if max(ncIndices) == 6
        ncStartFrameVec(end+1) = framesRaw(end)+1;
    else
        ncStartFrameVec(end+1) = anaphaseFrames(max(ncIndices)+1)+1;
    end
    if ncStartFrameVec(1) == 0 && ncStartFrameVec(2) > 1
        ncStartFrameVec(1) = 1;
    elseif ncStartFrameVec(1) == 0
        firstInd = find(ncStartFrameVec,1);
        if ncStartFrameVec(firstInd) == 1
            ncStartFrameVec = ncStartFrameVec(firstInd:end);
            ncIndices = ncIndices(firstInd:end);
        else
            ncStartFrameVec(firstInd-1) = 1;
            ncStartFrameVec = ncStartFrameVec(firstInd-1:end);
            ncIndices = ncIndices(firstInd-1:end);
        end
    end        
    firstTime = timeRaw(min(ncStartFrameVec));
    
    % iterate through nuclear cycles
    for ncInd = 1:length(ncStartFrameVec)-1
        
        % initialize structure to store nucleus and particle info
        compiledSchnitzCells = struct;
        
        nucleusCounter = 1;    %counter to ensure no empty spaces in the compiledSchnitzCells
        %struct if some schnitzcells have no frames in the
        %desired nuclear cycle
        
        firstNCFrame = ncStartFrameVec(ncInd);%max([1, processedData.(['nc' num2str(firstNC)])]);
        lastNCFrame = ncStartFrameVec(ncInd+1)-1;
        
        % filter reference vectors
        framesNC = framesRaw(firstNCFrame:lastNCFrame);
        tracesNC = tracesRaw(firstNCFrame:lastNCFrame,:);
        timeNC = timeRaw(firstNCFrame:lastNCFrame);
        timeNC = timeNC - firstTime; %normalize to start of first nc
        timeMovie = timeRaw(firstNCFrame:min((lastNCFrame+10), max(framesRaw)));
        %%%%%%%%%%%%%%%%%%%% Compile nucleus schnitz info %%%%%%%%%%%%%%%%%%%%%
        nc_ap_flags(i) = hasAPInfo && ~isfield(schnitzcells,'APpos');
        for s = 1:length(schnitzcells)
            schnitzFrames = schnitzcells(s).frames;
            
            % Filter for only the frames in this nc
            ncFilter = ismember(schnitzFrames,framesNC);
            rawNucleusFrames = schnitzFrames(ncFilter);
            
            if length(rawNucleusFrames) >= 1 & (schnitzcells(s).anaphaseFrame >= firstNCFrame-2) & (max(schnitzFrames) <= lastNCFrame+ 10)  %only grab data from the desired nuclear cycle(s)
               
                % Add info that's the same for all schnitzcells in this
                % experiment
                compiledSchnitzCells(nucleusCounter).setID = setID;
                compiledSchnitzCells(nucleusCounter).nc = ncRefVec(ncIndices(ncInd));
                
                % Initialize particle fields--will be set to particle values
                % for nuclei with matching particle
                compiledSchnitzCells = initializeParticleFields(...
                    compiledSchnitzCells,nucleusCounter,sum(ncFilter),has3DSpotInfo,hasAPInfo);
                
                % Add QC-related flags
                compiledSchnitzCells(nucleusCounter).TraceQCFlag = NaN;
                compiledSchnitzCells(nucleusCounter).FrameQCFlags = NaN(1,length(ncFilter));
                compiledSchnitzCells(nucleusCounter).N = NaN;
                compiledSchnitzCells(nucleusCounter).sparsity = NaN;
                
                % Add core nucleus info
                x = schnitzcells(s).cenx;
                y = schnitzcells(s).ceny;
                compiledSchnitzCells(nucleusCounter).xPosNucleus = x(ncFilter);
                compiledSchnitzCells(nucleusCounter).yPosNucleus = y(ncFilter);
                if hasAPInfo && isfield(schnitzcells,'APpos')
                    compiledSchnitzCells(nucleusCounter).apPosNucleus =  schnitzcells(s).APpos(ncFilter);
                end
                
                % Add protein and nucleus info
                if hasProteinInfo
                    compiledSchnitzCells(nucleusCounter).rawNCProtein = nanmax(schnitzcells(s).Fluo(ncFilter,:),[],2)';
                end
                compiledSchnitzCells(nucleusCounter).frames = rawNucleusFrames';
                compiledSchnitzCells(nucleusCounter).nucleusID = s;
                compiledSchnitzCells(nucleusCounter).ncID = eval([num2str(setID) '.' sprintf('%04d',s)]);
                compiledSchnitzCells(nucleusCounter).minDP = minDP;
                compiledSchnitzCells(nucleusCounter).minTime = minTime;
                compiledSchnitzCells(nucleusCounter).pctSparsity = pctSparsity;
                
                if ~isnan(schnitzcells(s).anaphaseFrame) & (schnitzcells(s).anaphaseFrame >= firstNCFrame-2) 
                    timeSchnitz = timeMovie - timeRaw(schnitzcells(s).anaphaseFrame);
                else
                    timeSchnitz = NaN(size(timeMovie));
                end
                % Add time and set info
                compiledSchnitzCells(nucleusCounter).time = timeNC(ismember(framesNC,rawNucleusFrames));
                if numel(unique(compiledSchnitzCells(nucleusCounter).time))~=numel(compiledSchnitzCells(nucleusCounter).time)
                    error(['Non-unique time values. Check FrameInfo for Prefix: ' currExperiment.Prefix])
                end
                
                compiledSchnitzCells(nucleusCounter).nuclear_time = timeSchnitz(ismember(framesNC,rawNucleusFrames));
                if numel(unique(compiledSchnitzCells(nucleusCounter).nuclear_time))~=numel(compiledSchnitzCells(nucleusCounter).nuclear_time)
                    error(['Non-unique nuclear time values. Check FrameInfo for Prefix: ' currExperiment.Prefix])
                end
                % initialize nucleus qc flag
                compiledSchnitzCells(nucleusCounter).missingNucleusFrames = 0;
                
                % Add schnitz cycle info -- weird errors arising without
                % this GM 7/18/22
                compiledSchnitzCells(nucleusCounter).schnitz_cycle = schnitzcells(s).cycle;
                % Add schnitz and particle approval--to be updated 
                compiledSchnitzCells(nucleusCounter).schnitz_approved = schnitzcells(s).Approved;
                compiledSchnitzCells(nucleusCounter).particle_approved = 1;
                % Add sister spot fields
                compiledSchnitzCells(nucleusCounter).sisterIndex = NaN;
                compiledSchnitzCells(nucleusCounter).sisterParticleID = NaN;
                
                % increment counter to ensure no empty spaces in the
                % compiledSchnitzCells struct if some schnitzcells have no
                % frames in the desired nuclear cycle
                nucleusCounter = nucleusCounter + 1;
            end
        end
        
        % Index vector to cross-ref w/ particles
        schnitzIndex = [compiledSchnitzCells.nucleusID];
        
        % Iterate through traces
        for j = 1:size(tracesNC,2)
            % Raw fluo trace
            rawTrace = tracesNC(:,j);
            
            % Get nucleus ID
            schnitz = compiledParticles(j).schnitz;
            
            % Find first and last expression frames, requiring that spots
            % cannot apear earlier than some fixed time (minTime)
            traceStart = find(timeNC>=minTime&~isnan(rawTrace'),1);
            traceStop = find(~isnan(rawTrace),1,'last');
            
            % skip particles not in specified nc range
            if isempty(traceStart)
                continue
            end
            
            %Create versions with all intervening frames present (missing frames
            %appear as NaNs)
            traceFull = rawTrace(traceStart:traceStop)';
            traceFramesFull = framesNC(traceStart:traceStop);
            
            % Perform qc tests
            nDP = sum(~isnan(traceFull));
            sparsity = prctile(diff(find(~isnan(traceFull))),pctSparsity);
            TraceQCFlag = nDP >= minDP && sparsity == 1;
            
            % trace-nucleus mapping may be many-to-1
            ncIndex = find(schnitzIndex==schnitz);
            if length(ncIndex) ~= 1
                continue
                warning(['Problem with particle-nucleus crossreference for Prefix: ' currExperiment.Prefix])
            end
            
            % Identifier variable
            particle = compiledParticles(j).OriginalParticle;
            particleID = eval([num2str(setID) '.' sprintf('%04d',particle)]);
            
            % check to see if a different particle has already been assigned
            % if so, create a new entry
            if ~isnan(compiledSchnitzCells(ncIndex).particleID) && twoSpotFlag
                compiledSchnitzCells(end+1) = compiledSchnitzCells(ncIndex);
                
                % update cross-reference variables
                compiledSchnitzCells(ncIndex).sisterParticleID = particleID;
                compiledSchnitzCells(end).sisterParticleID = compiledSchnitzCells(ncIndex).sisterParticleID;
                compiledSchnitzCells(ncIndex).sisterIndex = length(compiledSchnitzCells);
                compiledSchnitzCells(end).sisterIndex = ncIndex;
                
                % redefine index variables
                ncIndex = length(compiledSchnitzCells);
                
                % reset particle fields
                nEntries = length(compiledSchnitzCells(ncIndex).xPosParticle);
                compiledSchnitzCells(ncIndex).ParticleID = NaN;
                
                % Initialize particle fields--will be set to particle values
                % for nuclei with matching particle
                compiledSchnitzCells = initializeParticleFields(...
                    compiledSchnitzCells,ncIndex,nEntries);%,threeDFlag,APFlag,DetrendedZFlag);
            end
            
            % Record particle identifier
            compiledSchnitzCells(ncIndex).particleID = particleID;
            
            % Find overlap between nucleus and trace
            rawNucleusFrames = compiledSchnitzCells(ncIndex).frames;
            spotFilter = ismember(rawNucleusFrames,traceFramesFull);
            
            compiledSchnitzCells(ncIndex).spotFrames = spotFilter;
            if sum(spotFilter) < length(traceFramesFull)
                warning(['Inconsistent particle and nucleus frames for Prefix: ' currExperiment.Prefix ' .Flagging...'])
                compiledSchnitzCells(ncIndex).missingNucleusFrames = 1;
            end
            
            % record fluorescence info
            compiledSchnitzCells(ncIndex).fluo(spotFilter) = traceFull(ismember(traceFramesFull,rawNucleusFrames));
            
            % Find intersection btw full frame range and CP frames
            rawParticleFrames = compiledParticles(j).Frame;
            %             rawParticleFrames = rawParticleFrames(ismember(rawParticleFrames,traceFramesFull));
            
            % Add Manual Flagging Info for particle approval
            compiledSchnitzCells(ncIndex).particle_approved = compiledParticles(j).ManualApproved;
            % make filters
            ncSpotFilter1 = ismember(rawNucleusFrames,rawParticleFrames);
            ncSpotFilter2 = ismember(rawParticleFrames,rawNucleusFrames);
            
            % add offset info
            compiledSchnitzCells(ncIndex).fluoOffset(ncSpotFilter1) = compiledParticles(j).Off(ncSpotFilter2);
            % x, y, and z info
            compiledSchnitzCells(ncIndex).xPosParticle(ncSpotFilter1) = compiledParticles(j).xPos(ncSpotFilter2);
            compiledSchnitzCells(ncIndex).yPosParticle(ncSpotFilter1) = compiledParticles(j).yPos(ncSpotFilter2);
            compiledSchnitzCells(ncIndex).zPosParticle(ncSpotFilter1) = compiledParticles(j).zPos(ncSpotFilter2);
            if hasAPInfo
                compiledSchnitzCells(ncIndex).APPosParticle(ncSpotFilter1) = compiledParticles(j).APposParticle(ncSpotFilter2)*100;
            end
            % 3D info
            if has3DSpotInfo & isfield(compiledParticles, 'xPos3D')
                compiledSchnitzCells(ncIndex).xPosParticle3D(ncSpotFilter1) = compiledParticles(j).xPos3D(ncSpotFilter2);
                compiledSchnitzCells(ncIndex).yPosParticle3D(ncSpotFilter1) = compiledParticles(j).yPos3D(ncSpotFilter2);
                compiledSchnitzCells(ncIndex).zPosParticle3D(ncSpotFilter1) = compiledParticles(j).zPos3D(ncSpotFilter2);
                compiledSchnitzCells(ncIndex).fluo3D(ncSpotFilter1) = compiledParticles(j).Fluo3DRaw(ncSpotFilter2);
            end
            if DetrendedZFlag
                compiledSchnitzCells(ncIndex).zPosParticleDetrended(ncSpotFilter1) = compiledParticles(j).zPosDetrended(ncSpotFilter2);
            end
            % add qc info
            compiledSchnitzCells(ncIndex).N = nDP;
            compiledSchnitzCells(ncIndex).sparsity = sparsity;
            compiledSchnitzCells(ncIndex).TraceQCFlag = TraceQCFlag;
            compiledSchnitzCells(ncIndex).FrameQCFlags(ncSpotFilter1) = TraceQCFlag;
        end
        if length(compiledSchnitzCells) ~= length([compiledSchnitzCells.particleID])
            error('wtf')
        end
        spot_struct = [spot_struct  compiledSchnitzCells];
    end

    waitbar(i/numExperiments,h, ['Compiling data for dataset ' num2str(i) ' of ' num2str(numExperiments)])
end
close(h)

spot_struct = spot_struct(ismember([spot_struct(:).schnitz_cycle], firstNC:lastNC));
spot_struct = spot_struct([spot_struct(:).schnitz_approved] == 1 & [spot_struct(:).particle_approved] == 1);
% add nucleus AP info if appropriate
if any(nc_ap_flags)
    for i=1:length(spot_struct)    
        spot_struct(i).APPosNucleus = NaN(size(spot_struct(i).yPosNucleus));
    end
    set_vec = [spot_struct.setID];        
    for i = 1:numExperiments
        set_filter = set_vec==i;
        currExperiment = liveProject.includedExperiments{i};
        spot_struct(set_filter) = convertToFractionalEmbryoLength(currExperiment.Prefix,spot_struct(set_filter));
    end
end

% add additional fields
% NL: this does nothing at the moment, but leaving to preserve two-spot
% compatibility
for i = 1:numel(spot_struct)
    spot_struct(i).twoSpotFlag = twoSpotFlag;
    spot_struct(i).targetLocusFlag = NaN;
    spot_struct(i).controlLocusFlag = NaN;
end

disp('interpolating data...')

% generate interpolation fields
interpFields = {'fluo','xPosParticle','yPosParticle','zPosParticle',};
if has3DSpotInfo
    interpFields(end+1:end+4) = {'fluo3D','xPosParticle3D','yPosParticle3D','zPosParticle3D'};
end
if hasAPInfo
    interpFields(end+1:end+2) = {'APPosParticle','APPosNucleus'};
end
if hasProteinInfo
    interpFields(end+1) = {'rawNCProtein'};
end

% calculate interpolation time
%%
% sample true experimental res times
dt_vec_full = [];
multi_step_traces = find([spot_struct.N]>2);
for i = randsample(multi_step_traces,min([50,length(multi_step_traces)]),false)
    t_vec = spot_struct(i).time;
    dt_vec_full = [dt_vec_full diff(t_vec)];
end

% calculate rounded dt
dt_round = floor(median(dt_vec_full)/5)*5;
if exist('dt', 'var')
    if dt < dt_round
        warning('Specified dt implies a higher time resolution than is used in the original data.')
    end
    dt_round = dt;
end

tresInterp = max([tresInterpFloor,dt_round]);
interpGrid = 0:tresInterp:120*60;

% use 95th percentile of fluorescence values to set scale for flagging
% unusual jump events
fluo_scale = prctile([spot_struct.fluo],97.5);
big_blip_thresh = fluo_scale;

for i = 1:length(spot_struct)
    
    % First identify isolated "islands" at start or end of trace. These
    % tend to be artifactual
    fluoVec = spot_struct(i).fluo;
    frameIndex = 1:length(spot_struct(i).frames);
    %         spotFrames = spot_struct(i).spotFrames;
    
    % NL: these sizes generally work ok, but may need to change if time res
    % is >> or  << ~20 seconds or reporates with elongation times >> or <<
    % 2 min
    startIndex = [];
    stopIndex = [];
    
    vecNans = ~isnan(fluoVec);
    startIndexRaw = find(vecNans,1);
    stopIndexRaw = find(vecNans,1,'last');
    nFramesRaw = stopIndexRaw-startIndexRaw+1;
    
    % this is designed to get rid of blips at start and end of traces
    if length(fluoVec) > 1
        bSize = 8;
        kernelBig = ones(1,bSize);
        sSize = 3;
        kernelSmall = ones(1,sSize);
        
        vecConvBig = conv(kernelBig,vecNans);
        vecConvSmall = conv(kernelSmall,vecNans);
        
        % generate starting flags
        vecConvStartBig = vecConvBig(bSize:end);
        vecConvStartSmall = vecConvSmall(sSize:end);
        
        startIndex = find(vecConvStartSmall>=sSize | vecConvStartBig>=floor(bSize/2.1) & vecNans,1);
        
        % generate ending flags
        vecConvEndBig = vecConvBig(1:end-bSize+1);
        vecConvEndSmall = vecConvSmall(1:end-sSize+1);
        
        stopIndex = find(vecConvEndSmall>=sSize | vecConvEndBig>=floor(bSize/2.1) & vecNans,1,'last');
    end

    % update QC flags
    if ~isempty(startIndex) && ~isempty(stopIndex)
        spot_struct(i).FrameQCFlags(frameIndex<startIndex & vecNans) = false;
        spot_struct(i).FrameQCFlags(frameIndex>stopIndex & vecNans) = false;
    else
        spot_struct(i).FrameQCFlags(vecNans) = false;
    end
    
    if ~isnan(spot_struct(i).TraceQCFlag)
        spot_struct(i).TraceQCFlag =  spot_struct(i).TraceQCFlag && nansum(spot_struct(i).FrameQCFlags)>=minDP;
    end
    
    spot_struct(i).truncatedFlag = 0;
    if ~isempty(startIndex) && ~isempty(stopIndex)
        nFramesNew = stopIndex-startIndex+1;
        spot_struct(i).nFramesTruncated = nFramesRaw-nFramesNew;
        spot_struct(i).truncatedFlag = spot_struct(i).nFramesTruncated>0;
    elseif ~isempty(startIndexRaw) && ~isempty(stopIndexRaw)
        spot_struct(i).nFramesTruncated = nFramesRaw;
        spot_struct(i).truncatedFlag = 1;
    end
    
       % generate time reference vefctors for interpolation
    % generate time reference vefctors for interpolation
    timeVec = spot_struct(i).nuclear_time;
    timeVec = timeVec(startIndex:stopIndex);
    
    
    if length(timeVec)>1
        
        % add caps before and after
        first_ind = find(interpGrid>=timeVec(1),1);
        first_ind = max([1 first_ind-1]);
        last_ind = find(interpGrid<=timeVec(end),1,'last');
        last_ind = min([length(interpGrid) last_ind+1]);
        timeInterp = interpGrid(first_ind:last_ind);
                
        if isempty(timeInterp)
            timeInterp = NaN;
            for  j = 1:numel(interpFields)
                spot_struct(i).([interpFields{j} 'Interp']) = NaN;
            end
        else
            for  j = 1:length(interpFields)
                if ~isfield(spot_struct, interpFields{j})
                    continue 
                end
                vec = spot_struct(i).(interpFields{j})(startIndex:stopIndex);
                
                % perform basic QC on fluorescence fields
                if contains(interpFields{j},'fluo')
                    
                    %Look for clusters of 6 or more NaNs
                    kernel = [1,1,1,1,1];
                    vecNans = isnan(vec);
                    vecConv = conv(kernel,vecNans);
                    vecConv = vecConv(3:end-2);
                    zIDs = find(vecConv==numel(kernel));
                    zIDs = unique([zIDs-1 zIDs zIDs+1]); % get set of z_ids
                    vec(zIDs) = 0; % set clusters to zeros
                    vec(vec<0) = 0; % deal with negative values
                    
                    % find single dp "blips". These will be replaced via interpolation
                    % interpolate remaining NaNs
                    fluo_dd = [0 diff(vec,2) 0];
                    vec(abs(fluo_dd)>big_blip_thresh) = NaN;
                    
                    % find and remove suspsiciously large rises
                    fluo_d = [0 diff(vec,1) 0];
                    vec(abs(fluo_d)>big_blip_thresh*0.75) = NaN;
                    
                end
                
                % interpolate remaining NaNs
                referenceTime = timeVec(~isnan(vec));
                referenceVec = vec(~isnan(vec));
                
                if length(referenceTime)>1
                    % Interpolate to standardize spacing
                    vec_i = interp1(referenceTime,referenceVec,timeInterp,'linear','extrap');
                    vec_i(vec_i<0) = 0;
                    spot_struct(i).([interpFields{j} 'Interp']) = vec_i;
                elseif length(referenceTime)==1 && referenceTime >= timeInterp(1) && referenceTime <= timeInterp(end)
                    vecTo = NaN(size(timeInterp));
                    [~,mi] = min(abs(referenceTime-timeInterp));
                    vecTo(mi) = referenceVec;
                    spot_struct(i).([interpFields{j} 'Interp']) = vecTo;
                else
                    spot_struct(i).([interpFields{j} 'Interp']) = NaN(size(timeInterp));
                end
            end
        end
    else
        timeInterp = NaN;
        for  j = 1:numel(interpFields)
            spot_struct(i).([interpFields{j} 'Interp']) = NaN;
        end
    end  
    spot_struct(i).timeInterp = timeInterp;
    spot_struct(i).tresInterp = tresInterp;  
    
    % generate time reference vefctors for interpolation
    timeVec = spot_struct(i).time;
    timeVec = timeVec(startIndex:stopIndex);
    
    
    if length(timeVec)>1
        
        % add caps before and after
        first_ind = find(interpGrid>=timeVec(1),1);
        first_ind = max([1 first_ind-1]);
        last_ind = find(interpGrid<=timeVec(end),1,'last');
        last_ind = min([length(interpGrid) last_ind+1]);
        timeInterp = interpGrid(first_ind:last_ind);
                
        if isempty(timeInterp)
            timeInterp = NaN;
            for  j = 1:numel(interpFields)
                spot_struct(i).([interpFields{j} 'OrigInterp']) = NaN;
            end
        else
            for  j = 1:length(interpFields)
                if ~isfield(spot_struct, interpFields{j})
                    continue 
                end
                vec = spot_struct(i).(interpFields{j})(startIndex:stopIndex);
                
                % perform basic QC on fluorescence fields
                if contains(interpFields{j},'fluo')
                    
                    %Look for clusters of 6 or more NaNs
                    kernel = [1,1,1,1,1];
                    vecNans = isnan(vec);
                    vecConv = conv(kernel,vecNans);
                    vecConv = vecConv(3:end-2);
                    zIDs = find(vecConv==numel(kernel));
                    zIDs = unique([zIDs-1 zIDs zIDs+1]); % get set of z_ids
                    vec(zIDs) = 0; % set clusters to zeros
                    vec(vec<0) = 0; % deal with negative values
                    
                    % find single dp "blips". These will be replaced via interpolation
                    % interpolate remaining NaNs
                    fluo_dd = [0 diff(vec,2) 0];
                    vec(abs(fluo_dd)>big_blip_thresh) = NaN;
                    
                    % find and remove suspsiciously large rises
                    fluo_d = [0 diff(vec,1) 0];
                    vec(abs(fluo_d)>big_blip_thresh*0.75) = NaN;
                    
                end
                
                % interpolate remaining NaNs
                referenceTime = timeVec(~isnan(vec));
                referenceVec = vec(~isnan(vec));
                
                if length(referenceTime)>1
                    % Interpolate to standardize spacing
                    vec_i = interp1(referenceTime,referenceVec,timeInterp,'linear','extrap');
                    vec_i(vec_i<0) = 0;
                    spot_struct(i).([interpFields{j} 'OrigInterp']) = vec_i;
                elseif length(referenceTime)==1 && referenceTime >= timeInterp(1) && referenceTime <= timeInterp(end)
                    vecTo = NaN(size(timeInterp));
                    [~,mi] = min(abs(referenceTime-timeInterp));
                    vecTo(mi) = referenceVec;
                    spot_struct(i).([interpFields{j} 'OrigInterp']) = vecTo;
                else
                    spot_struct(i).([interpFields{j} 'OrigInterp']) = NaN(size(timeInterp));
                end
            end
        end
    else
        timeInterp = NaN;
        for  j = 1:numel(interpFields)
            spot_struct(i).([interpFields{j} 'OrigInterp']) = NaN;
        end
    end  
    spot_struct(i).origTimeInterp = timeInterp;
    spot_struct(i).tresInterp = tresInterp;    
    
    
 
end
%%


if lastNC ~= 14
    dataName = strrep(dataName, '\spot_struct.mat', [filesep, 'NC', num2str(lastNC),'\spot_struct.mat']);
    outpath = fileparts(dataName);
    if ~exist(outpath, 'dir')
        mkdir(outpath)
    end
end

if exist('dt', 'var')
    dataName = strrep(dataName, '\spot_struct.mat', ['\spot_struct_dt', num2str(dt), '.mat']);
end

% save
save(dataName ,'spot_struct')

if calculatePSF % NL: this needs to be fixed...
    disp('calculating psf dims...')
    % call function to calculate average psf difs
    calculate_average_psf(projectName,dropboxFolder);
end

disp('done.')
