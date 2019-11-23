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
surge_size_vec_ft = surge_size_vec(burst_ft_primary)';
mf_protein_vec_ft = mf_protein_vec(burst_ft_primary);
burst_dur_vec_ft = lag_dur_vec_target(burst_ft_primary);
set_vec_ft = set_vec(burst_ft_primary);

% calculate feature sizes and mf for two groups
mf_sat_vec = NaN(1,nBoots);
sz_sat_vec = NaN(1,nBoots);
mf_nsat_vec = NaN(1,nBoots);
sz_nsat_vec = NaN(1,nBoots);

sat_ids = find(ismember(set_vec_ft,sat_sets));
nsat_ids = find(~ismember(set_vec_ft,sat_sets));

for i = 1:nBoots
    s_ids = randsample(sat_ids,numel(sat_ids),true);
    ns_ids = randsample(nsat_ids,numel(nsat_ids),true);

    mf_sat_vec(i) = nanmean(mf_protein_vec_ft(s_ids));
    sz_sat_vec(i) = nanmean(surge_size_vec_ft(s_ids));

    mf_nsat_vec(i) = nanmean(mf_protein_vec_ft(ns_ids));
    sz_nsat_vec(i) = nanmean(surge_size_vec_ft(ns_ids));
end
    
mf_sat_mean = nanmean(mf_sat_vec);
mf_sat_se = nanstd(mf_sat_vec);
sz_sat_mean = nanmean(sz_sat_vec);
sz_sat_se = nanstd(sz_sat_vec);

mf_nsat_mean = nanmean(mf_nsat_vec);
mf_nsat_se = nanstd(mf_nsat_vec);
sz_nsat_mean = nanmean(sz_nsat_vec);
sz_nsat_se = nanstd(sz_nsat_vec);

% plot mean trend to get a feel for the relationship
nPoints = 15;
mf_index = linspace(prctile(mf_protein_vec_ft,1),prctile(mf_protein_vec_ft,95),nPoints+1);
surge_sz_vec = NaN(1,numel(mf_index)-1);

for m = 1:numel(mf_index)-1
    surge_sz_vec(m) = nanmean(surge_size_vec_ft(mf_protein_vec_ft<mf_index(m+1)&mf_protein_vec_ft>=mf_index(m)));
end
    

p = polyfit(mf_protein_vec_ft,surge_size_vec_ft,1);
surge_pd1 = polyval(p,mf_index);

close all
fig = figure;
hold on
plot(mf_index,surge_pd1);
s1 = scatter(mf_index(1:end-1),surge_sz_vec);
errorbar(mf_sat_mean,sz_sat_mean,-sz_sat_se,sz_sat_se,-mf_sat_se,mf_sat_se,'o','Color','black')
errorbar(mf_nsat_mean,sz_nsat_mean,-sz_nsat_se,sz_nsat_se,-mf_nsat_se,mf_nsat_se,'o','Color','black')
s2 = scatter(mf_sat,sz_sat,'filled');
s3 = scatter(mf_nsat,sz_nsat,'filled');
xlabel('average Dl concentration')
ylabel('surge size')
legend([s2 s3],'"saturating sets"','"non-saturating sets"','Location','northwest')
grid on
box on
xlim([90 275])
set(gca,'Fontsize',14)
% save
saveas(fig,[FigPath 'MCP_vs_burst_size_control.tif'])
saveas(fig,[FigPath 'MCP_vs_burst_size_control.pdf'])