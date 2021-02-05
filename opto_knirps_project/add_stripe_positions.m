% Script to incorporate curvature-adjusted AP positions 
% Run this after running main02
clear
close all


projectList = {'2xDl-Ven_snaBAC-mCh'}; % Cell array containing all projects you wish to process

% iterate through list of projects and generate required variables
for p = 1:length(projectList)
  
  % get liveProject object
  liveProject = LiveEnrichmentProject(projectList{p});
  
  % load compiled traces dataot
  disp('Loading spot structure...')
  load([liveProject.dataPath filesep 'spot_struct_protein.mat']);
  
  % load AP map file
  disp('load AP map file here...')
  
  % initialize grouping variables and add new position
  for i = 1:length(spot_struct_protein)
    % rename fields
    spot_struct_protein(i).APPosParticleInterpOrig = spot_struct_protein(i).APPosParticleInterp; % rename native AP variable
    spot_struct_protein(i).APPosParticleInterp = NaN(size(spot_struct_protein(i).xPosParticle));
    % extract position info
    xPosVec = spot_struct_protein(i).xPosParticleInterp;
    yPosVec = spot_struct_protein(i).yPosParticleInterp;
    zPosVec = spot_struct_protein(i).zPosParticleInterp;
    for j = 1:length(xPosVec)      
        %NL: do lookup here
        %spot_struct_protein(i).APPosParticleInterp(j) = something...
    end
  end
  
  % save
  save([liveProject.dataPath filesep 'spot_struct_protein.mat'],'spot_struct_protein');      
  
end