clear 
close all
addpath('utilities')

% set ID variables and paths

project = 'Rbp1-GFP_snaBAC-mCh';
DropboxFolder = 'S:\Nick\Dropbox\';

[~, DataPath, FigureRoot] =   header_function(DropboxFolder, project); 
FigPath = [FigureRoot '\' project '\exploratory_analyses\'];
mkdir(FigPath)

% load data
load([DataPath 'nucleus_struct_protein.mat'])

%% plot MS2 fluorescence vs Pol II
ms2_fluo_vec = [nucleus_struct_protein.fluo];
polII_fluo_vec = [nucleus_struct_protein.spot_protein_vec];

scatter_fig = figure;
scatter(ms2_fluo_vec,polII_fluo_vec)