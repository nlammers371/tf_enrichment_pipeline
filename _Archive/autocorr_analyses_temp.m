% Script to probe compare loca protein dynamics of control and target loci
clear 
close all
% define ID variables
project = 'Dl_Venus_snaBAC_MCPmCherry_Leica_Zoom2_7uW14uW';
dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
% dataPath = [dropboxFolder '\ProcessedEnrichmentData\' project '\'];
figPath = [dropboxFolder '\LocalEnrichmentFigures\' project '\'];
mkdir(figPath)
n_lags = 15;
% load data
load([dropboxFolder 'ProcessedEnrichmentData/' project '/nucleus_struct_protein.mat'])


%%
auto_mat_spot = NaN(n_lags+1,numel(nucleus_struct_protein));
auto_mat_null = NaN(n_lags+1,numel(nucleus_struct_protein));
auto_weights = NaN(n_lags+1,numel(nucleus_struct_protein));
for j = 1:numel(nucleus_struct_protein)
    time = nucleus_struct_protein(j).time;
    pt_vec_spot = nucleus_struct_protein(j).spot_protein_vec;
    pt_vec_null = nucleus_struct_protein(j).serial_null_protein_vec;     
    start_ind = find(~isnan(pt_vec_null)&~isnan(pt_vec_spot),1);
    stop_ind = find(~isnan(pt_vec_null)&~isnan(pt_vec_spot),1,'last');
    % perform qc checks
    pt_vec_null = pt_vec_null(start_ind:stop_ind);
    if numel(pt_vec_null) > n_lags+1 && sum(~isnan(pt_vec_null)) / numel(pt_vec_null) > .8
        time = time(start_ind:stop_ind);
        pt_vec_spot = pt_vec_spot(start_ind:stop_ind);
        nan_spot = isnan(pt_vec_spot);
        % interpolate
        pt_spot_interp = pt_vec_spot;
        pt_spot_interp(nan_spot) = interp1(time(~nan_spot),pt_vec_spot(~nan_spot),time(nan_spot));
        nan_null = isnan(pt_vec_null);
        pt_null_interp = pt_vec_null;
        pt_null_interp(nan_null) = interp1(time(~nan_null),pt_vec_null(~nan_null),time(nan_null));
        %
        ac_spot = xcov(pt_spot_interp,'coeff',n_lags);
        auto_mat_spot(:,j) = ac_spot(n_lags+1:end);
        ac_null = xcov(pt_null_interp,'coeff',n_lags);
        auto_mat_null(:,j) = ac_null(n_lags+1:end);
        auto_weights(:,j) = fliplr(numel(pt_spot_interp)-n_lags:numel(pt_spot_interp));
    end
end
        
