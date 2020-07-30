clear
close all

addpath('../utilities/')
% define basic path variables
DropboxFolder =  'S:\Nick\Dropbox (Personal)\ProcessedEnrichmentData\';
snail_project = '2xDl-Ven_snaBAC-mCh_v4';
np_project = 'Dl-Ven_hbP2P-mCh_v2';

% load snail data first
load([DropboxFolder snail_project '\nucleus_struct.mat'])
