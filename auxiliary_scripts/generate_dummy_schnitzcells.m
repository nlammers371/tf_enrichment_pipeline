function schnitzcells = generate_dummy_schnitzcells(projectName,varargin)

% %%%%%%%%%%%%%%%%%%%%% Set Path Specs, ID Vars %%%%%%%%%%%%%%%%%%%%%%%%%%
[liveProject, numExperiments, ~, hasAPInfo, ~, ~] = headerFunction(projectName);
%hasProteinInfo = false;

% %%%%%%%%%%%%%%%%% Extract relevant processed data %%%%%%%%%%%%%%%%%%%%%%

% Generate master structure with info on all nuclei and traces in
% constituent sets
% spot_struct = [];

% Add data from each experiment to the master structure
h = waitbar(0,'Generating schnitcells sets...');
% nc_ap_flags = false(1,numExperiments);

for i = 1:numExperiments
    
    waitbar((i)/numExperiments,h, ['Generating dummy set ' num2str(i) ' of ' num2str(numExperiments)])
        
    currExperiment = liveProject.includedExperiments{i};
    
    %%%%%%%% Read in processed data from main mRNADyanmics pipeline %%%%%%%
    
    % check to see if schnitzcells extists
    schnitzcells = getSchnitzcells(currExperiment);
    FluoDim = size(schnitzcells(1),2);
%     if ~isempty(schnitzcells)
%         save([currExperiment.resultsFolder, currExperiment.Prefix, '_lin_orig.mat'],'schnitzcells')
%     end
    processedData = getCompiledParticles(currExperiment);   % this structure contains all info compiled for this experiment
    compiledParticles = processedData.CompiledParticles;    % this is the cell array containing the CompiledParticles cell array
%     save([currExperiment.resultsFolder, 'CompiledParticles_orig.mat'],'CompiledParticles')
    if iscell(compiledParticles)
        compiledParticles = compiledParticles{1};               % this assumes there is only one channel with spots to analyze
    end        
    
    % generate new schnitzcells
    schnitzcells = struct;
    for c = 1:length(compiledParticles)
      
        schnitz = compiledParticles(c).schnitz;
        
        schnitzcells(schnitz).frames = compiledParticles(c).Frame(1):compiledParticles(c).Frame(end);
        schnitzcells(schnitz).APpos = NaN(size(schnitzcells(schnitz).frames));
        schnitzcells(schnitz).Fluo = NaN(length(schnitzcells(schnitz).frames),FluoDim);
        schnitzcells(schnitz).cenx = NaN(size(schnitzcells(schnitz).frames));
        schnitzcells(schnitz).ceny = NaN(size(schnitzcells(schnitz).frames));
        frame_filter = ismember(schnitzcells(schnitz).frames,compiledParticles(c).Frame);
        if hasAPInfo
            schnitzcells(schnitz).APpos(frame_filter) = compiledParticles(c).APPos;
        end
        schnitzcells(schnitz).cenx(frame_filter) = compiledParticles(c).xPos;
        schnitzcells(schnitz).ceny(frame_filter) = compiledParticles(c).yPos;
        
%         % update CompiledParticles 
%         compiledParticles(c)..schnitz
    end
    
    
    save([currExperiment.resultsFolder, currExperiment.Prefix, '_lin_dummy.mat'],'schnitzcells')
end    
delete(h)