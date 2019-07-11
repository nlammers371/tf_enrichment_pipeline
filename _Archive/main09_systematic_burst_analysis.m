clear 
close all

% define core ID variables
project = 'Dl-Ven_snaBAC-mCh';
% project = 'Dl-Ven_hbP2P-mCh';
dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
dataPath = [dropboxFolder 'ProcessedEnrichmentData\' project '/'];
% load input-output data set
load([dataPath 'hmm_burst_table.mat'])

% examine location of burst "center or mass" as a function of burst
% duration
burst_class_vec = results_table.burst_class;
burst_dur_vec = results_table.tr_dur;
burst_amp_vec = results_table.tr_amp;
burst_dur_prev_vec = results_table.tr_dur_prev;
pt_cm_spot_vec = results_table.pt_cm_spot;
pt_net2_spot_vec = results_table.pt_net2_spot;
pt_net_spot_vec = results_table.pt_net_spot;
pt_cm_swap_vec = results_table.pt_cm_swap;
mf_vec = results_table.pt_mf;
fluo_vec = results_table.fluo;
ft_vec = burst_class_vec==1& burst_dur_prev_vec >=6 & burst_dur_vec >=3 & ~isnan(pt_cm_spot_vec);
