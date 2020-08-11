function nucleus_struct = main01_compile_traces(dataStatusTab,dropboxFolder,varargin)

%
% nucleusStruct = main01_compile_traces(dataStatusTab,dropboxFolder,varargin)
%
% DESCRIPTION
% Function to compile relevant outputs from main mRNADynamics image 
% analysis pipeline across multiple experiments
%
%
% ARGUMENTS
% dataStatusTab: master ID variable (should match a tab name in the 
%                DataStatus.xlsx spreadsheet)
% dropboxFolder: full file path to folder containing compiled imaging
%                results
%
% OPTIONS
% firstNC: script defaults to taking only nc14 traces. If you want
%          earlier traces, pass 'firstNC', followed by desired nuclear 
%          cycle number
%
% OUTPUT
% nucleus_struct: compiled data set contain key nucleus and
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
%     rawNCPprotein: 
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

%% Set defaults
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

[rawResultsRoot, dataPath, ~] =   header_function(dropboxFolder, projectName);


%% %%%%%%%%%%%%%%%%%%%%% Set Path Specs, ID Vars %%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the DataStatus.xlsx file and grab only those datasets marked as
% approved by 'ReadyForEnrichmentAnalysis' flag
liveProject = LiveProject(projectName);
approvedFlag = 'ReadyForEnrichmentAnalysis';
approvedPrefixes = getCustomApprovedPrefixes(liveProject,approvedFlag);
prefixes = approvedPrefixes;
    
% Make the output filepath
mkdir(dataPath);
% Assign save names
nucleusName = [dataPath 'nucleus_struct.mat']; % names for compiled elipse struct


%% %%%%%%%%%%%%%%%%%%%%%%% Obtain Relevant Filepaths %%%%%%%%%%%%%%%%%%%%%%

% Get all the file paths for the relevant pipeline outputs
numExperiments = numel(prefixes);
compiledNucleiFilenames = cell(1,numExperiments); % nuclear protein
for i = 1:numExperiments
    currPrefix = prefixes{i};            
    compiledNucleiFilenames{i} = [rawResultsRoot, currPrefix, filesep, 'CompiledNuclei.mat'];      
end

% Generate set key data structure
setKey = array2table((1:numel(prefixes))','VariableNames',{'setID'});
% prefixesTable = cell2table(prefixes);
setKey.prefix = prefixes;
save([dataPath 'set_key.mat'],'setKey')

% Generate master structure with info on all nuclei and traces in
% constituent sets
h = waitbar(0,'Compiling data ...');
nucleus_struct = [];

% Loop through prefixes    
for i = 1:numExperiments
    waitbar(i/numExperiments,h, ['Compiling data for dataset ' num2str(i) ' of ' num2str(numExperiments)])
    
    setID = i;
    currExperiment = LiveExperiment(prefixes{i});
    
    % Read in processed data files from the main mRNADyanmics pipeline
    schnitzcells = getSchnitzcells(currExperiment);
    processedSpotData = getCompiledParticles(currExperiment);   % this structure contains all info compiled for this experiment
    compiledParticles = processedSpotData.CompiledParticles;    % this is the cell array containing the CompiledParticles cell array
    if iscell(compiledParticles)
        compiledParticles = compiledParticles{1};               % this is the CompiledParticles structure we actually want to use  
    end   
    
    % MT 2020-07-15 I don't think this is ever used ... maybe remove?
    processedNucleiData = load(compiledNucleiFilenames{i}); % processed nuclei 

    %%%%% Basic data characteristics %%%%%%
        
    % Check if certain info is present for this experiment    
    threeDFlag = isfield(compiledParticles,'Fluo3DRaw');    %3D spot fitting data        
    apFlag = isfield(compiledParticles, 'APpos');   %AP position info

    % Pull trace, time, and frame variables
    timeRaw = processedSpotData.ElapsedTime*60; %time vector, [sec] 
    framesRaw = 1:length(timeRaw); % Frame list
    tracesRaw = processedSpotData.AllTracesVector;  %array with a column for each trace
    if iscell(tracesRaw)        %check to see if traces are stored in cell array
        tracesRaw = tracesRaw{1};   
    end

    % Trim trace and time arrays to start at first frame of the first
    % nuclear cycle
    firstFrame = max([1, processedSpotData.(['nc' num2str(firstNC)])]);
    frames = framesRaw(firstFrame:end); 
    traces = tracesRaw(firstFrame:end,:);
    time = timeRaw(firstFrame:end);    
    time = time - min(time); %normalize to start of nc
    
    % temporary fix for problematic time vector
    if strcmp(prefixes{i},'2020-02-07-Dl_Ven_snaBAC_MCPmCh_F-F-F_Leica_Zoom2_7uW14uW_06') %wtf
        time = linspace(0,max(time),numel(time));
    end
       
    % Get basic frame info from the liveExperiment instance
    yDim = currExperiment.yDim;
    xDim = currExperiment.xDim;
    zDim = currExperiment.zDim;
    zStep = currExperiment.zStep_um;
    pixelSize = currExperiment.pixelSize_um; %um
    
    
    %%%%% Compile nucleus schnitz info %%%%%%
    compiledSchnitzCells = struct;
    
    sWithNC = 1;    %counter to ensure no empty spaces in the compiledSchnitzCells
                    %struct if some schnitzcells have no frames in the
                    %desired nuclear cycle
    for s = 1:length(schnitzcells)
        schnitzFrames = schnitzcells(s).frames;
        
        % Filter for only the frames in this nc
        ncFilter = ismember(schnitzFrames,frames);
        ncFrames = schnitzFrames(ncFilter); 
        
        if length(ncFrames) >= 1     %only grab data from the desired nuclear cycle(s)
         
            % Add info that's the same for all schnitzcells in this
            % experiment
            compiledSchnitzCells(sWithNC).setID = setID;
            compiledSchnitzCells(sWithNC).sourcePath = currExperiment.resultsFolder;
            compiledSchnitzCells(sWithNC).APFlag = 0;  %AP flag not relevant atm
            compiledSchnitzCells(sWithNC).threeDFlag = threeDFlag;
            compiledSchnitzCells(sWithNC).xDim = xDim;
            compiledSchnitzCells(sWithNC).yDim = yDim;
            compiledSchnitzCells(sWithNC).zDim = zDim;
            compiledSchnitzCells(sWithNC).zStep = zStep;
            compiledSchnitzCells(sWithNC).pixelSize = pixelSize;   %um
            
            % Initialize particle fields--will be set to particle values 
            % for nuclei with matching particle
            compiledSchnitzCells(sWithNC).particleID = NaN;
            compiledSchnitzCells(sWithNC).xPosParticle = NaN(1,sum(ncFilter));
            compiledSchnitzCells(sWithNC).yPosParticle = NaN(1,sum(ncFilter));
            compiledSchnitzCells(sWithNC).zPosParticle = NaN(1,sum(ncFilter)); 
            compiledSchnitzCells(sWithNC).APPosParticle = NaN(1,sum(ncFilter)); 
            compiledSchnitzCells(sWithNC).spotFrames = NaN(1,sum(ncFilter)); 
            if threeDFlag
                compiledSchnitzCells(sWithNC).xPosParticle3D = NaN(1,sum(ncFilter));
                compiledSchnitzCells(sWithNC).yPosParticle3D = NaN(1,sum(ncFilter));
                compiledSchnitzCells(sWithNC).zPosParticle3D = NaN(1,sum(ncFilter)); 
                compiledSchnitzCells(sWithNC).fluo3D = NaN(1,sum(ncFilter));
            end
            compiledSchnitzCells(sWithNC).fluo = NaN(1,sum(ncFilter)); 
            compiledSchnitzCells(sWithNC).fluoOffset = NaN(1,sum(ncFilter));  
            
            % Add QC-related flags
            compiledSchnitzCells(sWithNC).qcFlag = NaN;
            compiledSchnitzCells(sWithNC).N = NaN;
            compiledSchnitzCells(sWithNC).sparsity = NaN;
            
            % Add core nucleus info
            x = schnitzcells(s).cenx;            
            y = schnitzcells(s).ceny;
            compiledSchnitzCells(sWithNC).xPosNucleus = x(ncFilter);
            compiledSchnitzCells(sWithNC).yPosNucleus = y(ncFilter);
            
            % Add protein info
            compiledSchnitzCells(sWithNC).rawNCPprotein = nanmax(schnitzcells(s).Fluo(ncFilter,:),[],2);
            compiledSchnitzCells(sWithNC).frames = ncFrames';            
            compiledSchnitzCells(sWithNC).nucleusID = s;     
            compiledSchnitzCells(sWithNC).ncID = eval([num2str(setID) '.' sprintf('%04d',s)]);
            compiledSchnitzCells(sWithNC).ncStart = firstNC;
            compiledSchnitzCells(sWithNC).minDP = minDP;
            compiledSchnitzCells(sWithNC).minTime = minTime;
            compiledSchnitzCells(sWithNC).pctSparsity = pctSparsity;
            
            % Add time and set info
            compiledSchnitzCells(sWithNC).time = time(ismember(frames,ncFrames));
            if numel(unique(compiledSchnitzCells(sWithNC).time))~=numel(compiledSchnitzCells(sWithNC).time)
                error(['Non-unique time values. Check FrameInfo for Prefix: ' currExperiment.Prefix])
            end

            % Add sister spot fields
            compiledSchnitzCells(sWithNC).sisterIndex = NaN;
            compiledSchnitzCells(sWithNC).sisterParticleID = NaN; 
            
            % increment counter to ensure no empty spaces in the 
            % compiledSchnitzCells struct if some schnitzcells have no 
            % frames in the desired nuclear cycle
            sWithNC = sWithNC + 1;
            
        else
            %Do nothing with the schnitzcells not in the desired nuclear
            %cycle(s)
        end
    end
    
    % Index vector to cross-ref w/ particles            
    schnitzIndex = [compiledSchnitzCells.nucleusID];          
   
    % Iterate through traces
    for j = 1:size(traces,2)  
        % Raw fluo trace
        rawTrace = traces(:,j); 
        % Get nucleus ID
        schnitz = compiledParticles(j).schnitz;    
        % Find first and last expression frames, requiring that spots
        % cannot apear earlier than some fixed time (minTime)
        traceStart = find(time>=minTime&~isnan(rawTrace'),1);
        traceStop = find(~isnan(rawTrace),1,'last');        
        %Creat versions with all intervening frames present (missing frames
        %appear as NaNs)
        traceFull = rawTrace(traceStart:traceStop)';                                   
        traceFramesFull = frames(traceStart:traceStop); 
        % skip particles not in specified nc range
        if sum(~isnan(traceFull)) == 0
            continue
        end   
        
        % Find intersection btw full frame range and CP frames        
        rawProteinFrames = compiledParticles(j).Frame;
        rawProteinFrames = rawProteinFrames(ismember(rawProteinFrames,traceFramesFull));
        % Perform qc tests    
        nDP = sum(~isnan(traceFull));
        sparsity = prctile(diff(find(~isnan(traceFull))),pctSparsity);
        qcFlag = nDP >= minDP && sparsity == 1;     
        % trace-nucleus mapping may be many-to-1
        ncIndex = find(schnitzIndex==schnitz);  
        if length(ncIndex) ~= 1
            warning('Problem with Particle-Nucleus Crossref')
            continue
        end 
        
        % Identifier variable                        
        particle = compiledParticles(j).OriginalParticle;                        
        particleID = eval([num2str(setID) '.' sprintf('%04d',particle)]);
        
        % check to see if a different particle has already been assigned
        % if so, create a new entry
        if ~isnan(compiledSchnitzCells(ncIndex).particleID)
            twoSpotFlag = true;
            compiledSchnitzCells = [compiledSchnitzCells compiledSchnitzCells(ncIndex)];
            % update cross-reference variables
            compiledSchnitzCells(ncIndex).sisterParticleID = particleID;
            compiledSchnitzCells(end).sisterParticleID = compiledSchnitzCells(ncIndex).sisterParticleID;
            compiledSchnitzCells(ncIndex).sisterIndex = numel(compiledSchnitzCells);
            compiledSchnitzCells(end).sisterIndex = ncIndex;
            % redefine index variables
            ncIndex = numel(compiledSchnitzCells);
            % reset particle fields
            nEntries = numel(compiledSchnitzCells(ncIndex).xPosParticle);
            compiledSchnitzCells(ncIndex).ParticleID = NaN;
            compiledSchnitzCells(ncIndex).xPosParticle = NaN(1,nEntries);
            compiledSchnitzCells(ncIndex).yPosParticle = NaN(1,nEntries);
            compiledSchnitzCells(ncIndex).zPosParticle = NaN(1,nEntries);
            if threeDFlag
                compiledSchnitzCells(ncIndex).xPosParticle3D = NaN(1,sum(ncFilter));
                compiledSchnitzCells(ncIndex).yPosParticle3D = NaN(1,sum(ncFilter));
                compiledSchnitzCells(ncIndex).zPosParticle3D = NaN(1,sum(ncFilter));
                compiledSchnitzCells(ncIndex).fluo3D = NaN(1,nEntries);
            end
            compiledSchnitzCells(ncIndex).fluo = NaN(1,nEntries);    
            compiledSchnitzCells(ncIndex).fluoOffset = NaN(1,nEntries);  
        end
        % Record particle identifier
        compiledSchnitzCells(ncIndex).particleID = particleID;
        % find overlap between nucleus and trace
        ncFrames = compiledSchnitzCells(ncIndex).frames;         
        spotFilter = ismember(ncFrames,traceFramesFull);            
        compiledSchnitzCells(ncIndex).spotFrames = spotFilter;
        if sum(spotFilter) < numel(traceFramesFull)            
            error('Inconsistent particle and nucleus frames')
        end
        % record fluorescence info             
        compiledSchnitzCells(ncIndex).fluo(spotFilter) = traceFull; % note that
        % make filters
        ncSpotFilter1 = ismember(ncFrames,rawProteinFrames);
        ncSpotFilter2 = ismember(rawProteinFrames,ncFrames);
        % add offset info
        compiledSchnitzCells(ncIndex).fluoOffset(ncSpotFilter1) = compiledParticles(j).Off(ncSpotFilter2);
        % x, y, and z info                                
        compiledSchnitzCells(ncIndex).xPosParticle(ncSpotFilter1) = compiledParticles(j).xPos(ncSpotFilter2);
        compiledSchnitzCells(ncIndex).yPosParticle(ncSpotFilter1) = compiledParticles(j).yPos(ncSpotFilter2);   
        compiledSchnitzCells(ncIndex).zPosParticle(ncSpotFilter1) = compiledParticles(j).zPos(ncSpotFilter2);
        if apFlag
            compiledSchnitzCells(ncIndex).APPosParticle(ncSpotFilter1) = compiledParticles(j).APposParticle(ncSpotFilter2)*100;
        end
        % 3D info
        if threeDFlag
            compiledSchnitzCells(ncIndex).xPosParticle3D(ncSpotFilter1) = compiledParticles(j).xPosGauss3D(ncSpotFilter2);            
            compiledSchnitzCells(ncIndex).yPosParticle3D(ncSpotFilter1) = compiledParticles(j).yPosGauss3D(ncSpotFilter2);              
            compiledSchnitzCells(ncIndex).zPosParticle3D(ncSpotFilter1) = compiledParticles(j).zPosGauss3D(ncSpotFilter2);            
            compiledSchnitzCells(ncIndex).fluo3D(ncSpotFilter1) = compiledParticles(j).Fluo3DRaw(ncSpotFilter2);
        end
        % add qc info
        compiledSchnitzCells(ncIndex).N = nDP;
        compiledSchnitzCells(ncIndex).sparsity = sparsity;        
        compiledSchnitzCells(ncIndex).qcFlag = qcFlag;  
    end      
    nucleus_struct = [nucleus_struct  compiledSchnitzCells];        
end
close(h)

% add additional fields
for i = 1:numel(nucleus_struct)
    nucleus_struct(i).twoSpotFlag = twoSpotFlag;
    nucleus_struct(i).targetLocusFlag = NaN;
    nucleus_struct(i).controlLocusFlag = NaN;
end    

disp('interpolating data...')
% generate interpolation fields
if threeDFlag
    interpFields = {'fluo','fluo3D'};
else
    interpFields = {'fluo'};
end

interpGrid = 0:tresInterp:60*60;

for i = 1:numel(nucleus_struct)
    timeVec = nucleus_struct(i).time;
    fluoVec = nucleus_struct(i).fluo;
    startIndex = find(~isnan(fluoVec),1);
    stopIndex = find(~isnan(fluoVec),1,'last');    
    timeVec = timeVec(startIndex:stopIndex);       
    if numel(timeVec)>1%nucleus_struct(i).qc_flag == 1
        timeInterp = interpGrid(interpGrid>=timeVec(1)&interpGrid<=timeVec(end));
        for  j = 1:numel(interpFields)
            vec = nucleus_struct(i).(interpFields{j})(startIndex:stopIndex);
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
            queryPoints = timeVec(isnan(vec));
            interpT = timeVec(~isnan(vec));
            interpF = vec(~isnan(vec));
            if ~isempty(queryPoints)
                newF = interp1(interpT,interpF,queryPoints);  
                vec(ismember(timeVec,queryPoints)) = newF;   
            else
                vec = interpF;
            end                 
            % Interpolate to standardize spacing   
            try
                nucleus_struct(i).([interpFields{j} 'Interp']) = interp1(timeVec,vec,timeInterp);
            catch
                error('why?')
            end
        end
    else
        timeInterp = NaN;
        for  j = 1:numel(interpFields)
            nucleus_struct(i).([interpFields{j} 'Interp']) = NaN;
        end
    end
    nucleus_struct(i).timeInterp = timeInterp;
    nucleus_struct(i).tresInterp = tresInterp;
end
% call function to calculate average psf difs
% save
save(nucleusName ,'nucleus_struct') 

if calculatePSF
    disp('calculating psf dims...')
    calculate_average_psf(projectName,DropboxFolder);
end

disp('done.')
