function [liveProject, numExperiments, nucleusName, hasAPInfo, has3DSpotInfo] = headerFunction(projectName)

% check to see if there is a full mRNADynamics repo on the working path. If
% so remove it to prevent function cross-talk
cleanUpmRNADynamics;
       
% Determine key project characteristics
liveProject = LiveProject(projectName);
numExperiments = length(liveProject.includedExperimentNames);

% Make the output filepath
mkdir(liveProject.dataPath);

% Assign save names
nucleusName = [liveProject.dataPath 'nucleus_struct.mat']; % names for compiled elipse struct

% Check for optiuonal parameter felds
hasAPInfo = all(liveProject.hasAPInfo);
has3DSpotInfo = all(liveProject.has3DSpotInfo);