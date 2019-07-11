clear 
close all

% define core ID variables
project = 'Dl-Ven_snaBAC-mCh';
% project = 'Dl-Ven_hbP2P-mCh';
dropboxFolder =  'E:\Nick\Dropbox\';
dataPath = [dropboxFolder 'ProcessedEnrichmentData\' project '/'];
% load input-output data set
K = 3;
w = 7;
load([dataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'],'hmm_input_output')

%% Examin variation in burst duration and frequency as a function of [Dl]

z_dur_lag_vec_master = [hmm_input_output.z_dur_lag_vec];
z_diff_vec_master = [hmm_input_output.z_diff_vec];
mf_vec_master = [hmm_input_output.mf_protein];

% apply filter
burst_dur_vec = z_dur_lag_vec_master(z_diff_vec_master==1);
mf_burst_dur_vec = mf_vec_master(z_diff_vec_master==1);

burst_sep_vec = z_dur_lag_vec_master(z_diff_vec_master==-1);
mf_burst_sep_vec = mf_vec_master(z_diff_vec_master==-1);

mf_grid = 0:5:400;
burst_sep_mean = NaN(size(mf_grid));
burst_dur_mean = NaN(size(mf_grid));
mf_sigma = 5;
for m = 1:numel(mf_grid)
    d_mf_dur = exp(-(mf_grid(m)-mf_burst_dur_vec).^2 / 2 /mf_sigma^2);
    d_mf_sep = exp(-(mf_grid(m)-mf_burst_sep_vec).^2 / 2 /mf_sigma^2);
    
    burst_sep_mean(m) = nansum(burst_sep_vec.*d_mf_sep) ./ nansum(d_mf_sep);
    burst_dur_mean(m) = nansum(burst_dur_vec.*d_mf_dur) ./ nansum(d_mf_dur);
end



