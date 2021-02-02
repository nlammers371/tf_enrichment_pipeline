% Script to incorporate curvature-adjusted AP positions 
% Run this after running main02
clear
close all


projectList = {'opto_knirps_WT'}; % Cell array containing all projects you wish to process

% iterate through list of projects and generate required variables
for p = 1:length(projectList)
  
  % get liveProject object
  liveProject = LiveEnrichmentProject(projectList{p});
  
  % load compiled traces data
  load([liveProject.dataPath filesep 'spot_struct_protein.mat']);
  
  % load AP map file
  disp('load AP map file here...')
  
  % initialize grouping variables
  for i = 1:length(spot_struct_protein)
    spot_struct_protein(i).APPosParticleInterpOrig = spot_struct_protein(i).APPosParticleInterp; % rename native AP variable
    spot_struct_protein(i).APPosParticleInterp = NaN(size(spot_struct_protein(i).xPosParticle));
  end
  
  % generate reference vector of setIDS
  setIDVec = [spot_struct_protein.setID];
  particleIDVec = [spot_struct_protein.particleID];
    
  % save
  save([liveProject.dataPath filesep 'spot_struct_protein.mat'],'spot_struct_protein');
  save([liveProject.dataPath filesep 'customDataOptions.mat'],'customDataOptions');
  
  
  
end