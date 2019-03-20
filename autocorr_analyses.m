% Script to probe compare loca protein dynamics of control and target loci
clear 
close all
% define ID variables
project_cell = {'Dl_Venus_snaBAC_mCherry_Leica_hp','Dl_Venus_hbP2P_mCherry_LeicaZoom2x_HighPower'};
dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
primary_project = project_cell{1};
% dataPath = [dropboxFolder '\ProcessedEnrichmentData\' project '\'];
figPath = [dropboxFolder '\LocalEnrichmentFigures\' primary_project '\'];
mkdir(figPath)
n_lags = 15;
% loop through each project and extract local protein vectors
autocorr_cell = cell(1,numel(project_cell));
for i = 1:numel(project_cell)
    load([dropboxFolder 'ProcessedEnrichmentData/' project_cell{i} '/nucleus_struct_protein.mat'])
    auto_mat = NaN(n_lags+1,numel(nucleus_struct_protein));
    auto_weights = NaN(n_lags+1,numel(nucleus_struct_protein));
    for j = 1:numel(nucleus_struct_protein)
        pt_vec = nucleus_struct_proteinj(j).spot_protein_vec;
        
        
end