% Script to build data set for systematic input/output analyses
clear 
close all
% define ID variables
K = 3;
w = 6;
project = 'Dl_Venus_snaBAC_MCPmCherry_Leica_Zoom2_7uW14uW';
dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
% dropboxFolder = 'C:\Users\nlamm\Dropbox (Garcia Lab)\';
dataPath = [dropboxFolder '\ProcessedEnrichmentData\' project '\'];
% load data set
load([dataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'])
% set size of time series "reads"
window_size_time = 5*60;
Tres = hmm_input_output(1).Tres;
window_size_steps = round(window_size_time/Tres);
pt_sm_kernel = 1; % size of kernel used to smooth protein channel
%%
min_peak_prom_mcp = .05

test = hmm_input_output(15).r_vec;
