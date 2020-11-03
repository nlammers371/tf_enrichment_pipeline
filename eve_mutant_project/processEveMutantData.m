% Script to generate custom stripe-based grouping variables for eveBAC
% mutant datasets
% Run this after running main01
clear
close all

startTime = 25*60;

% list of project names
projectList = {'eveWT','EveGtSL','EveS1Null','EveGtSL-S1Null'};
ectopicGroups = {[],[1.5 -1],[1.5 -1],[1.5 -1]};
endogenousGroups = {1:3,1:3,2:3,2:3};
customDataFolders = {'EveWT_CompiledParticles_segmented',...
  'eveGt_Analysis2.0_Segmented\eveGtSL-CompiledParticles2.0',...
  'eveS1Null_CompiledParticles_analysis',...
  'CompiledParticles-eveGtSL-eS1_2.1'};




% iterate through list of projects and generate required variables
for p = 2%1:length(projectList)
  
  liveProject = LiveProject(projectList{p});
  
  % load compiled traces data
  load([liveProject.dataPath filesep 'spot_struct.mat']);
  
  % initialize grouping variables
  for i = 1:length(spot_struct)
    spot_struct(i).APPosParticleNorm = NaN(size(spot_struct(i).xPosParticle));
    spot_struct(i).StripeVec = NaN(size(spot_struct(i).xPosParticle));
    spot_struct(i).ectopicFlag = NaN;
    spot_struct(i).Stripe = NaN;
  end
  
  % generate reference vector of setIDS
  setIDVec = [spot_struct.setID];
  particleIDVec = [spot_struct.particleID];
  
  if liveProject.hasTracesCompiled
      % pull normalized AP positions from custom-segmented data sets
      numExperiments = length(liveProject.includedExperiments);
      resultsDir = liveProject.includedExperiments{1}.userResultsFolder;
      customDataPath = [resultsDir filesep customDataFolders{p} filesep];

      % iterate through each individual experiment
      for e = 1:numExperiments
        currExperiment = liveProject.includedExperiments{e};

        % dearch for corresponding custom data file
        customFile = dir([customDataPath filesep '*' currExperiment.Prefix '*']);
        if length(customFile) == 1
          % load custom file
          load([customFile.folder filesep customFile.name]);
          if iscell(CompiledParticles)
            CompiledParticles = CompiledParticles{1};
          end
          origParticleList = [CompiledParticles.OriginalParticle];
          % get list of indices in combined data set that correspond to this
          % embryo
          expIndices = find(setIDVec==e & ~isnan(particleIDVec));
          particleIDs = particleIDVec(expIndices);
          % now loop through each particle
          for ind = expIndices
            % cross reference with custom set
            ParticleID = spot_struct(ind).particleID;

            originalParticle = round(ParticleID * 1e4 - e*1e4);          

            % cross reference frames and particle position
            frameFilterTo = ismember(spot_struct(ind).frames,CompiledParticles(origParticleList==originalParticle).Frame);
            frameFilterFrom = ismember(CompiledParticles(origParticleList==originalParticle).Frame,spot_struct(ind).frames);

            % quick consistency check
            if all(spot_struct(ind).xPosParticle(frameFilterTo) == CompiledParticles(origParticleList==originalParticle).xPos(frameFilterFrom))
              % assign            
              if strcmp(projectList{p},'EveGtSL-S1Null')
                spot_struct(ind).APPosParticleNorm(frameFilterTo) = CompiledParticles(origParticleList==originalParticle).APPos_Normalized(frameFilterFrom);
              end                            
              spot_struct(ind).StripeVec(frameFilterTo) = CompiledParticles(origParticleList==originalParticle).Stripe(frameFilterFrom);
            else
              error('Problem with cross-referencing')
            end
          end
        else
          error(['Could not find custom file for Prefix: ' currExperiment.Prefix])
        end    
      end

      % perform project-specific post-processing
      endogenousRegions = endogenousGroups{p};
      ectopicRegions = ectopicGroups{p};

      for i = 1:length(spot_struct)
          if ~isnan(spot_struct(i).particleID)
              timeFilter = spot_struct(i).time >= startTime;

              spot_struct(i).Stripe = mode(spot_struct(i).StripeVec(timeFilter));
              
              if strcmp(projectList{p},'EveGtSL-S1Null')
                meanAP = mean(spot_struct(i).APPosParticleNorm(timeFilter));
                % readjust boundary between 1.5 and 2
                if meanAP <= 0.36 && spot_struct(i).Stripe == 2
                  spot_struct(i).Stripe = 1.5;
                elseif meanAP > 0.36 && spot_struct(i).Stripe == 1.5
                  spot_struct(i).Stripe = 2;
                end
              end

              % assign ectopic/endogenous status
              if ismember(spot_struct(i).Stripe,endogenousRegions)
                spot_struct(i).ectopicFlag = false;
              elseif ismember(spot_struct(i).Stripe,ectopicRegions) 
                spot_struct(i).ectopicFlag = true;
              end
          end
      end
      
      % record data options
      customDataOptions.endogenousRegions = endogenousRegions;
      customDataOptions.ectopicRegions = ectopicRegions;
      customDataOptions.startTime = startTime;

      % save
      save([liveProject.dataPath filesep 'spot_struct.mat'],'spot_struct');
      save([liveProject.dataPath filesep 'customDataOptions.mat'],'customDataOptions');
  end
  
  
end