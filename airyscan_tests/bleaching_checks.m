clear
close all 

projectName = 'Bcd-GFP_hbMS2-mCh_AiryscanTest';

% get paths 
liveProject = LiveEnrichmentProject(projectName);
resultsRoot = [liveProject.dataPath filesep];

% load data
load([resultsRoot filesep 'spot_struct.mat'])

%% Check Frame over Frame trends, first for a single set
