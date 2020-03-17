% Script to estimate elongation time using average auto-correlation of MS2
% traces
clear
close all
addpath('../utilities')
% define core ID variables
project = 'Dl-Ven_snaBAC-mCh_v3';
DropboxFolder =  'S:\Nick\Dropbox\';
[~, DataPath, FigureRoot] =   header_function(DropboxFolder, project);

FigPath = [FigureRoot '\absolute_calibration\' project '\'];
mkdir(FigPath)

load([DataPath 'nucleus_struct.mat'])
load("S:\Nick\Dropbox\ProcessedEnrichmentData\absolute_calibration\calibration_info.mat")
% assumptions
max_rnap = 70;

% calculate approximate calibration
fluo_vec = [nucleus_struct.fluo];
prctile99 = prctile(fluo_vec(~isnan(fluo_vec)),99);

au_per_rnap = prctile99/max_rnap

% transcription calibration (E left)
sna_burst_axis = [.4 .8 1.2] / au_per_rnap * 60


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dorsal calibration
VoxelSize = (nucleus_struct(1).PixelSize)^2*nucleus_struct(1).zStep;
integration_disk_size = 9;
% dorsal trace(C)
trace_axis = [-1:1] / calibration_info.venus_au_per_molecule * integration_disk_size

% surge (E left)
surge_axis = round([-20:10:20] * VoxelSize / calibration_info.venus_au_per_molecule * integration_disk_size,2)

% surge amplitude
surge_amp_axis = round([-.5:.5:1.5]/ calibration_info.venus_au_per_molecule * integration_disk_size / 6,2)

sna_burst_axis2 = [2 3 4] / au_per_rnap / 3 
