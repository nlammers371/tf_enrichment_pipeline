% Script to generate figures establishing presence of enrichment bursts at
% start of transcription bursts
clear 
% close all
addpath('utilities')
% set ID variables
targetProject = 'Dl-Ven_snaBAC-mCh';
DropboxFolder = 'E:\Nick\LivemRNA\Dropbox (Personal)\';
[~, DataPath, FigureRoot] =   header_function(DropboxFolder, targetProject); 

% load data
load([DataPath 'hmm_input_output_results.mat'])
FigPath = [FigureRoot 'data_controls\'];
mkdir(FigPath)

% analysis parameters
Tres = 20; % seconds
roi_window = 6; 
window_size = 15;
start = window_size + 2;
nBoots = 100; % number of bootstrap samples to use
min_pause_len = 6; % minimum length of preceding OFF period (in time steps)
min_burst_len = 2;
sat_sets = [1 2 4 6];
% extract relevant arrays from target project 
lag_dur_vec_target = results_struct.lag_dur_vec;
lead_dur_vec_target = results_struct.lead_dur_vec;
hmm_array = results_struct.hmm_array;
spot_array_dt = results_struct.spot_array_dt;
surge_size_vec = nansum(results_struct.spot_array_dt(:,start:start + roi_window),2);
mf_protein_vec = results_struct.mf_protein_vec;
feature_sign_vec_target = results_struct.feature_sign_vec;
set_vec = floor(results_struct.particle_id_vec);
% generate basic filter for target locus and computational controls
burst_ft_primary = feature_sign_vec_target == 1&lead_dur_vec_target>=min_pause_len&lag_dur_vec_target>min_burst_len; % filter for rise events

%% run basic regressions to determine predictive power of each covariate
% apply feature filter
surge_size_vec_ft = surge_size_vec(burst_ft_primary);
mf_protein_vec_ft = mf_protein_vec(burst_ft_primary)';
burst_dur_vec_ft = lag_dur_vec_target(burst_ft_primary)';

lm_mf = fitlm(mf_protein_vec_ft,surge_size_vec_ft)



lm_both = fitlm([mf_protein_vec_ft burst_dur_vec_ft],surge_size_vec_ft)
