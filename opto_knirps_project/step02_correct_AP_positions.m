% Script to incorporate curvature-adjusted AP positions 
% Run this after running main02
clear
close all

projectList = {'optokni_eve4+6_ON'}; % Cell array containing all projects you wish to process

% iterate through list of projects and generate required variables
for p = 1:length(projectList)
  
  % get liveProject object
  liveProject = LiveEnrichmentProject(projectList{p});
  
  % load compiled traces dataot
  disp('loading spot structure...')
  load([liveProject.dataPath filesep 'spot_struct.mat']);
  load([liveProject.dataPath filesep 'spot_struct_protein.mat']);
  
  % load AP map file
  disp('load AP map file here...')
  num_exp = length(liveProject.includedExperimentNames);
  for i = 1:num_exp
      ap_map(i) = load(['ap_map' filesep] + liveProject.includedExperimentNames(i) + ['_APmap.mat']);
  end
  
  disp('performing AP correction...')
 
  % initialize grouping variables and add new position
  for i = 1:length(spot_struct_protein)
    % rename fields
    spot_struct_protein(i).APPosParticleInterpOrig = spot_struct_protein(i).APPosParticleInterp; % rename native AP variable
    spot_struct_protein(i).APPosParticleInterp = NaN(size(spot_struct_protein(i).xPosParticleInterp));
    % extract position info
    xPosVec = spot_struct_protein(i).xPosParticleInterp;
    yPosVec = spot_struct_protein(i).yPosParticleInterp;
    zPosVec = spot_struct_protein(i).zPosParticleInterp;
    
    expID = spot_struct_protein(i).setID;
    
    for j = 1:length(xPosVec)    
        %NL: do lookup here
        if isnan(xPosVec)
            spot_struct_protein(i).APPosParticleInterp(j) = NaN;
        else
            spot_struct_protein(i).APPosParticleInterp(j) = ap_map(expID).APmap(ceil(yPosVec(j)),ceil(xPosVec(j)));
        end
    end
  end

  
  for i = 1:length(spot_struct)
    % rename fields
    spot_struct(i).APPosParticleInterpOrig = spot_struct(i).APPosParticleInterp; % rename native AP variable
    spot_struct(i).APPosParticleInterp = NaN(size(spot_struct(i).xPosParticleInterp));
    
    spot_struct(i).APPosNucleuspOrig = spot_struct(i).APPosNucleus; % rename native AP variable
    spot_struct(i).APPosNucleus = NaN(size(spot_struct(i).xPosNucleus));
    
    % extract position info
    xPosVec = spot_struct(i).xPosParticleInterp;
    yPosVec = spot_struct(i).yPosParticleInterp;
    zPosVec = spot_struct(i).zPosParticleInterp;
    
    xPosNucVec = spot_struct(i).xPosNucleus;
    yPosNucVec = spot_struct(i).yPosNucleus;
    
    expID = spot_struct(i).setID;
    
    for j = 1:length(xPosVec)    
        %NL: do lookup here
        if isnan(xPosVec)
            spot_struct(i).APPosParticleInterp(j) = NaN;
        else
            spot_struct(i).APPosParticleInterp(j) = ap_map(expID).APmap(ceil(yPosVec(j)),ceil(xPosVec(j)));
        end
    end
    
    for j = 1:length(xPosNucVec)    
        %NL: do lookup here
        if isnan(xPosNucVec)
            spot_struct(i).APPosNucleus(j) = NaN;
        else
            spot_struct(i).APPosNucleus(j) = ap_map(expID).APmap(ceil(yPosNucVec(j)),ceil(xPosNucVec(j)));
        end
    end
    
  end
  
  % save
  save([liveProject.dataPath filesep 'spot_struct_protein_corrected.mat'],'spot_struct_protein');
  save([liveProject.dataPath filesep 'spot_struct_corrected.mat'],'spot_struct');
  
end