% script to call functions that build set of reference data fro parameter
% sweeps
clear
close all

% build reactivation set
projectNameRA = 'optokni_eve4+6_ON'; 
disp('Building reactivation reference set...')
tic
io_ref_ra = build_ra_ref_set(projectNameRA);
toc
disp('Done.')

% build WT set
projectNameWT = 'optokni_eve4+6_WT_FUN';  
disp('Building wild-type reference set...')
tic
io_ref_wt = build_wt_ref_set(projectNameWT);
toc
disp('Done.')

% save combined dataset to master directory
liveProject = LiveEnrichmentProject(projectNameWT);
resultsPath = [liveProject.dataPath];
slashes = strfind(resultsPath,'\');
ResultsRoot = resultsPath(1:slashes(end-1));
sweepDir = [ResultsRoot 'parameterSweeps' filesep];
mkdir(sweepDir);
disp('Saving...')
save([sweepDir 'io_ref_wt.mat'],'io_ref_wt')
save([sweepDir 'io_ref_ra.mat'],'io_ref_ra')
disp('Done.')
save([sweepDir 'io_ref_wt.mat'],'io_ref_wt')

