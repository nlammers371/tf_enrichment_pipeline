function nucleus_struct = main01_compile_traces(projectName,varargin)

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
%
% other: script allows any default variable to be set using format:
%       "VariableNameString", VariableValue
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


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Set defaults %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all force
addpath(genpath('utilities'));
% Defaults
firstNC = 14;   % first nuclear cycle to pull data from
minDP = 15;     % what is this for?
pctSparsity = 50;   % 
twoSpotFlag = contains(projectName, '2spot');
minTime = 0*60; % take no fluorescence data prior to this point
tresInterp = 20; 
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
        elseif strcmpi(varargin{i}, 'projectName')
            error('Pipeline currently does not support alternate names for projects (otehr than what is on the data status tab)')
        elseif i < numel(varargin)                         
            eval([varargin{i} '=varargin{i+1};']);        
        end
    end
end

%% %%%%%%%%%%%%%%%%%%%%% Set Path Specs, ID Vars %%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% if isempty(livemRNAPath)
%   livemRNAPath = FindLivemRNA;
% end

% % copy ComputerFolders to local directory 
% copyfile([livemRNAPath '\ComputerFolders.csv'],'.');
    
% Find the DataStatus.xlsx file and grab only those datasets marked as
% approved by 'ReadyForEnrichmentAnalysis' flag
approvedFlag = 'ReadyForEnrichmentAnalysis';
[approvedPrefixes, dropboxFolder] = getProjectPrefixes(projectName,'customApproved',approvedFlag);
numExperiments = length(approvedPrefixes);
    
% Make the output filepath
dataPath = [dropboxFolder, projectName, '/'];
mkdir(dataPath);

% Assign save names
nucleusName = [dataPath 'nucleus_struct.mat']; % names for compiled elipse struct

% Generate set key data structure
setKey = array2table((1:numExperiments)','VariableNames',{'setID'});
% prefixesTable = cell2table(prefixes);
setKey.prefix = approvedPrefixes;
save([dataPath,'set_key.mat'],'setKey')

%% %%%%%%%%%%%%%%%%% Extract relevant processed data %%%%%%%%%%%%%%%%%%%%%%

% Generate master structure with info on all nuclei and traces in
% constituent sets
nucleus_struct = [];

% Add data from each experiment to the master structure
h = waitbar(0,'Compiling data ...');
for i = 1:numExperiments
    
    waitbar(i/numExperiments,h, ['Compiling data for dataset ' num2str(i) ' of ' num2str(numExperiments)])
    
    setID = i;
    currExperiment = LiveExperiment(approvedPrefixes{i});
    
    
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
    threeDFlag = isfield(compiledParticles,'Fluo3DRaw');    %3D spot fitting data        
    APFlag = isfield(compiledParticles, 'APpos');   %AP position info
    DetrendedZFlag = isfield(compiledParticles, 'zPosDetrended');   %Z pos that accounts for stack repositioning

    % Pull trace, time, and frame variables
    timeRaw = processedData.ElapsedTime*60; %time vector, [sec] 
    framesRaw = 1:length(timeRaw); % Frame list
    tracesRaw = processedData.AllTracesVector;  %array with a column for each trace
    if iscell(tracesRaw)        %check to see if traces are stored in cell array
        tracesRaw = tracesRaw{SpotChannelIndex};   
    end
    
    % Find start frames for nc's we wish to include
    ncIndices = find(ncRefVec>=firstNC);
    ncStartFrameVec = currExperiment.anaphaseFrames(ncIndices)';
    ncStartFrameVec(end) = max([1 ncStartFrameVec(end)]); 
    ncStartFrameVec(end+1) = framesRaw(end);
    firstTime = timeRaw(min(ncStartFrameVec));
    
    % Get basic frame info from the liveExperiment instance
    zStep = currExperiment.zStep_um;
    pixelSize = currExperiment.pixelSize_um; %um    
    
    % initialize structure to store nucleus and particle info
    compiledSchnitzCells = struct;
    
    sWithNC = 1;    %counter to ensure no empty spaces in the compiledSchnitzCells
                    %struct if some schnitzcells have no frames in the
                    %desired nuclear cycle
    
    % iterate through nuclear cycles
    for ncInd = 1:length(ncStartFrameVec)-1
        firstNCFrame = ncStartFrameVec(ncInd);%max([1, processedData.(['nc' num2str(firstNC)])]);      
        lastNCFrame = ncStartFrameVec(end);

        % filter reference vectors
        framesNC = framesRaw(firstNCFrame:lastNCFrame); 
        tracesNC = tracesRaw(firstNCFrame:lastNCFrame,:);
        timeNC = timeRaw(firstNCFrame:lastNCFrame);    
        timeNC = timeNC -firstTime; %normalize to start of first nc


        %%%%%%%%%%%%%%%%%%%% Compile nucleus schnitz info %%%%%%%%%%%%%%%%%%%%%

        for s = 1:length(schnitzcells)
            schnitzFrames = schnitzcells(s).frames;

            % Filter for only the frames in this nc
            ncFilter = ismember(schnitzFrames,framesNC);
            rawNucleusFrames = schnitzFrames(ncFilter); 

            if length(rawNucleusFrames) >= 1     %only grab data from the desired nuclear cycle(s)

                % Add info that's the same for all schnitzcells in this
                % experiment
                compiledSchnitzCells(sWithNC).setID = setID;
                compiledSchnitzCells(sWithNC).sourcePath = currExperiment.resultsFolder;
                compiledSchnitzCells(sWithNC).APFlag = 0;  %AP flag not relevant atm
                compiledSchnitzCells(sWithNC).threeDFlag = threeDFlag;
                compiledSchnitzCells(sWithNC).zStep = zStep;
                compiledSchnitzCells(sWithNC).pixelSize = pixelSize;   %um
                compiledSchnitzCells(sWithNC).nc = ncIndices(ncInd);

                % Initialize particle fields--will be set to particle values 
                % for nuclei with matching particle
                compiledSchnitzCells = initializeParticleFields(...
                  compiledSchnitzCells,sWithNC,sum(ncFilter));%,threeDFlag,APFlag,DetrendedZFlag);              
                
                % Add QC-related flags
                compiledSchnitzCells(sWithNC).qcFlag = NaN;
                compiledSchnitzCells(sWithNC).N = NaN;
                compiledSchnitzCells(sWithNC).sparsity = NaN;

                % Add core nucleus info
                x = schnitzcells(s).cenx;            
                y = schnitzcells(s).ceny;
                compiledSchnitzCells(sWithNC).xPosNucleus = x(ncFilter);
                compiledSchnitzCells(sWithNC).yPosNucleus = y(ncFilter);

                % Add protein and nucleus info
                compiledSchnitzCells(sWithNC).rawNCPprotein = nanmax(schnitzcells(s).Fluo(ncFilter,:),[],2);
                compiledSchnitzCells(sWithNC).frames = rawNucleusFrames';            
                compiledSchnitzCells(sWithNC).nucleusID = s;     
                compiledSchnitzCells(sWithNC).ncID = eval([num2str(setID) '.' sprintf('%04d',s)]);
  %               compiledSchnitzCells(sWithNC).ncStart = firstNC;
                compiledSchnitzCells(sWithNC).minDP = minDP;
                compiledSchnitzCells(sWithNC).minTime = minTime;
                compiledSchnitzCells(sWithNC).pctSparsity = pctSparsity;

                % Add time and set info
                compiledSchnitzCells(sWithNC).time = timeNC(ismember(framesNC,rawNucleusFrames));
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
            qcFlag = nDP >= minDP && sparsity == 1;     

            % trace-nucleus mapping may be many-to-1
            ncIndex = find(schnitzIndex==schnitz);  
            if length(ncIndex) ~= 1
                error(['Problem with particle-nucleus crossreference for Prefix: ' currExperiment.Prefix])          
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
                error(['Inconsistent particle and nucleus frames for Prefix: ' currExperiment.Prefix])  
            end
            
            % record fluorescence info             
            compiledSchnitzCells(ncIndex).fluo(spotFilter) = traceFull; 
                        
            % Find intersection btw full frame range and CP frames        
            rawParticleFrames = compiledParticles(j).Frame;
            rawParticleFrames = rawParticleFrames(ismember(rawParticleFrames,traceFramesFull));
            
            % make filters
            ncSpotFilter1 = ismember(rawNucleusFrames,rawParticleFrames);
            ncSpotFilter2 = ismember(rawParticleFrames,rawNucleusFrames);
            % add offset info
            compiledSchnitzCells(ncIndex).fluoOffset(ncSpotFilter1) = compiledParticles(j).Off(ncSpotFilter2);
            % x, y, and z info                                
            compiledSchnitzCells(ncIndex).xPosParticle(ncSpotFilter1) = compiledParticles(j).xPos(ncSpotFilter2);
            compiledSchnitzCells(ncIndex).yPosParticle(ncSpotFilter1) = compiledParticles(j).yPos(ncSpotFilter2);   
            compiledSchnitzCells(ncIndex).zPosParticle(ncSpotFilter1) = compiledParticles(j).zPos(ncSpotFilter2);
            if APFlag
                compiledSchnitzCells(ncIndex).APPosParticle(ncSpotFilter1) = compiledParticles(j).APposParticle(ncSpotFilter2)*100;
            end
            % 3D info
            if threeDFlag
                compiledSchnitzCells(ncIndex).xPosParticle3D(ncSpotFilter1) = compiledParticles(j).xPosGauss3D(ncSpotFilter2);            
                compiledSchnitzCells(ncIndex).yPosParticle3D(ncSpotFilter1) = compiledParticles(j).yPosGauss3D(ncSpotFilter2);              
                compiledSchnitzCells(ncIndex).zPosParticle3D(ncSpotFilter1) = compiledParticles(j).zPosGauss3D(ncSpotFilter2);            
                compiledSchnitzCells(ncIndex).fluo3D(ncSpotFilter1) = compiledParticles(j).Fluo3DRaw(ncSpotFilter2);
            end
            if DetrendedZFlag
              compiledSchnitzCells(ncIndex).zPosParticleDetrended(ncSpotFilter1) = compiledParticles(j).zPosDetrended(ncSpotFilter2);            
            end            
            % add qc info
            compiledSchnitzCells(ncIndex).N = nDP;
            compiledSchnitzCells(ncIndex).sparsity = sparsity;        
            compiledSchnitzCells(ncIndex).qcFlag = qcFlag;  
        end      
        nucleus_struct = [nucleus_struct  compiledSchnitzCells];
    end
end
close(h)

% add additional fields
% NL: this does nothing at the moment, but leaving to preserve two-spot
% compatibility
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

% use 95th percentile of fluorescence values to set scale for flagging 
% unusual jump events
fluo_scale = prctile([nucleus_struct.fluo],97.5);
big_blip_thresh = fluo_scale;

for i = 1:numel(nucleus_struct)
    timeVec = nucleus_struct(i).time;
    fluoVec = nucleus_struct(i).fluo;
    startIndex = find(~isnan(fluoVec),1);
    stopIndex = find(~isnan(fluoVec),1,'last');    
    timeVec = timeVec(startIndex:stopIndex);       
    if length(timeVec)>1
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
%             vec(vec<0) = 0; % deal with negative values    
            
            % find single dp "blips". These will be replaced via interpolation
            % interpolate remaining NaNs  
            fluo_dd = [0 diff(vec,2) 0];
            vec(abs(fluo_dd)>big_blip_thresh) = NaN;
            
            % find and remove suspsiciously large rises
            fluo_d = [0 diff(vec,1) 0];
            vec(abs(fluo_d)>big_blip_thresh*0.75) = NaN;

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
            nucleus_struct(i).([interpFields{j} 'Interp']) = interp1(timeVec,vec,timeInterp);
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

% save
save(nucleusName ,'nucleus_struct') 

if calculatePSF % NL: this needs to be fixed...
    disp('calculating psf dims...')
    % call function to calculate average psf difs
    calculate_average_psf(projectName,dropboxFolder);
end

% revmove temproary copy of ComputerFolders
% if ~strcmpi(pwd,livemRNAPath)
%   delete ComputerFolders.csv
% end
disp('done.')
