% script to call functions that build set of reference data fro parameter
% sweeps
clear
close all
addpath('functions');

dataRoot = 'S:\Jake\Dropbox\ProcessedEnrichmentData\';

% build reactivation set
projectNameRA = 'optokni_eve4+6_ON'; 
disp('Building reactivation reference set...')
tic
io_ref_ra = build_ra_ref_set_v2(projectNameRA,dataRoot);
toc
disp('Done.')

% build combined io dataset
disp('Building combined io set...')
tic
io_ref_cm = build_combined_io_set(dataRoot);
toc
disp('Done.')


disp('Saving...')
dataRootOut = 'S:\Nick\Dropbox\ProcessedEnrichmentData\';
sweepDir = [dataRootOut 'parameterSweeps_v2' filesep];
mkdir(sweepDir);
save([sweepDir 'io_ref_cm.mat'],'io_ref_cm')
save([sweepDir 'io_ref_ra.mat'],'io_ref_ra')
disp('Done.')

