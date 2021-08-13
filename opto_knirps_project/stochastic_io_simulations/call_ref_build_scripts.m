% script to call functions that build set of reference data fro parameter
% sweeps
clear
close all

% build reactivation set
projectNameRA = 'optokni_eve4+6_ON'; 
disp('Building reactivation reference set...')
io_ref_ra = build_ra_ref_set(projectNameRA);
disp('Done.')

% build WT set
projectNameWT = 'optokni_eve4+6_WT_FUN';  
disp('Building wild-type reference set...')
io_ref_wt = build_wt_ref_set(projectNameWT);
disp('Done.')

% save combined dataset to master directory


