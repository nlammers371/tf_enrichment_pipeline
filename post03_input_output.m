% Script to probe relationship between input protein concentration and
% output transcriptional response
clear 
close all
% define ID variables
K = 3;
w = 7;
project = 'Dl-Ven x hbP2P';
nBoots = 100;
% project = 'Dl_Venus_hbP2P_MCPmCherry_Zoom2_7uW14uW';
dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
% dropboxFolder = 'C:\Users\nlamm\Dropbox (Garcia Lab)\';
dataPath = [dropboxFolder '\ProcessedEnrichmentData\' project '\'];
figPath = [dropboxFolder '\LocalEnrichmentFigures\' project '\hmm_input_output_K' num2str(K) '_w' num2str(w) '\'];
mkdir(figPath)
load([dataPath 'input_output_snips.mat'])
gene_name = 'hbP2P';
protein_name = 'Dorsal';
%%
%%% Make time-dependent cross-covariance plots
% define some colors
yw = [234 194 100]/256; % yellow
bl = [115 143 193]/256; % blue
rd = [213 108 85]/256; % red
gr = [191 213 151]/256; % green
br = [207 178 147]/256; % brown
n_lags = floor(numel(input_output_snips(3).hmm_ctrl_xcov)/2);
feature_cell = {'fluo_rise','fluo_fall','protein_peak','protein_trough'};
flag_names = {'fluo_change_flags','fluo_change_flags','pt_peak_flag','pt_trough_flag'};
val_vec = [1 -1 1 1];


results_struct = struct;
for i = 1:numel(feature_cell)
    f_string = feature_cell{i};
    % initialize arrays
    protein_spot_diff_mat = NaN(nBoots,2*n_lags+1);
    protein_spot_mat = NaN(nBoots,2*n_lags+1);
    protein_mf_mat = NaN(nBoots,2*n_lags+1);
    protein_swap_spot_mat = NaN(nBoots,2*n_lags+1);
    protein_serial_mat = NaN(nBoots,2*n_lags+1);
    protein_swap_serial_mat = NaN(nBoots,2*n_lags+1);
    fluo_spot_mat = NaN(nBoots,2*n_lags+1);
    fluo_swap_spot_mat = NaN(nBoots,2*n_lags+1); 
        
    feature_indices = find([input_output_snips.(flag_names{i})]==val_vec(i));
    for n = 1:nBoots    
        % sample indices  
        boot_indices = randsample(feature_indices,numel(feature_indices),true);
        % extract
        p_vec_spot = nanmean(vertcat(input_output_snips(boot_indices).spot_protein_vec));
        p_vec_mf = nanmean(vertcat(input_output_snips(boot_indices).mf_protein_vec));
        p_vec_serial = nanmean(vertcat(input_output_snips(boot_indices).serial_protein_vec));
        
        
        p_vec_spot_swap = nanmean(vertcat(input_output_snips(boot_indices).swap_spot_protein_vec));                
        p_vec_mf_swap = nanmean(vertcat(input_output_snips(boot_indices).swap_mf_protein_vec));
        p_vec_serial_swap = nanmean(vertcat(input_output_snips(boot_indices).swap_serial_protein_vec));
        
        f_vec_spot = nanmean(vertcat(input_output_snips(boot_indices).fluo_vec));
        f_vec_swap = nanmean(vertcat(input_output_snips(boot_indices).swap_fluo_vec));
        % record        
        protein_mf_mat(n,:) = p_vec_mf;
        protein_spot_mat(n,:) = p_vec_spot;
        protein_serial_mat(n,:) = p_vec_serial;
        protein_swap_spot_mat(n,:) = p_vec_spot_swap;        
        fluo_spot_mat(n,:) = f_vec_spot;
        fluo_swap_spot_mat(n,:) = f_vec_swap;
    end
    % save to structure
    results_struct(i).feature = f_string;
    
    results_struct(i).protein_mf_mean = nanmean(protein_mf_mat);
    results_struct(i).protein_mf_ste = nanstd(protein_mf_mat);
    
    results_struct(i).protein_spot_mean = nanmean(protein_spot_mat);
    results_struct(i).protein_spot_ste = nanstd(protein_spot_mat);
    
    results_struct(i).protein_serial_mean = nanmean(protein_serial_mat);
    results_struct(i).protein_serial_ste = nanstd(protein_serial_mat);
    
    results_struct(i).protein_diff_mean = nanmean(protein_spot_mat - protein_serial_mat);
    results_struct(i).protein_diff_ste = nanstd(protein_spot_mat - protein_serial_mat);
    
    results_struct(i).protein_swap_spot_mean = nanmean(protein_swap_spot_mat);
    results_struct(i).protein_swap_spot_ste = nanstd(protein_swap_spot_mat);
    
    results_struct(i).fluo_spot_mean = nanmean(fluo_spot_mat);
    results_struct(i).fluo_spot_ste = nanstd(fluo_spot_mat);
    
    results_struct(i).fluo_swap_spot_mean = nanmean(fluo_swap_spot_mat);
    results_struct(i).fluo_swap_spot_ste = nanstd(fluo_swap_spot_mat);
    
    results_struct(i).fluo_spot_diff_mean = nanmean(fluo_spot_mat-fluo_swap_spot_mat);
    results_struct(i).fluo_spot_diff_ste = nanstd(fluo_spot_mat-fluo_swap_spot_mat);
end    

save([dataPath 'input_output_results.mat'],'results_struct')
