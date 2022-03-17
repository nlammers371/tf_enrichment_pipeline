% Script to generate figures establishing presence of enrichment bursts at
% start of transcription bursts
clear 
% close all
addpath('../utilities')
% set ID variables
DropboxFolder = 'S:\Nick\Dropbox\';
project = 'Dl-Ven_snaBAC-mCh_v4';
[~, DataPath, FigureRoot] =   header_function(DropboxFolder, project); 

fluo_dim = 2;
protein_dim = 2;

% Params
K = 3;
w = 7;
% load data
% final results
load([DataPath 'hmm_input_output_results_w' num2str(w) '_K' num2str(K) '_f' num2str(fluo_dim) 'D_p' num2str(protein_dim) 'D.mat'])
% intermediate input/output set
load([DataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '_f' num2str(fluo_dim) 'D_p' num2str(protein_dim) 'D_dt.mat'],'hmm_input_output')
% raw compiled data
load([DataPath 'nucleus_struct.mat'])

% make figure directory
FigPath = [FigureRoot '\input_output_driver\' project '\'];
mkdir(FigPath)


%% create analysis filters
analysis_struct = struct;
Tres = 20; % seconds
min_pause_len = 5; % minimum length of preceding OFF period (in time steps)
max_pause_len = 1000;
min_burst_len = 2;
max_burst_len = 1000;

lag_dur_vec = results_struct.lag_dur_vec;
lead_dur_vec = results_struct.lead_dur_vec;
feature_sign_vec = results_struct.feature_sign_vec;
% generate basic filter for target locus and computational controls
burst_ft = results_struct.feature_sign_vec == 1&results_struct.lead_dur_vec>=min_pause_len&results_struct.lead_dur_vec<=max_pause_len...
    &results_struct.lag_dur_vec>=min_burst_len&results_struct.lag_dur_vec<=max_burst_len;%    

% extract raw vectors
hmm_array = results_struct.hmm_array(burst_ft,:);
fluo_array = results_struct.fluo_array(burst_ft,:);
mf_array = results_struct.mf_array(burst_ft,:);     
time_vec = results_struct.center_time_vec(burst_ft)/60;
spot_array_dt = results_struct.spot_array_dm(burst_ft,:);

% create surge size vector, along with potential explanatory vectors
window_size = floor(size(hmm_array,2)/2);
window_vec = -window_size:window_size;
time_axis = 20*(window_vec);
b_ft = ismember(window_vec,-5:0);
t_ft = ismember(window_vec,1:6);

pt_surge_size_vec = nanmean(spot_array_dt(:,t_ft),2) - nanmean(spot_array_dt(:,b_ft),2);
mf_protein_vec = nanmean(mf_array(:,b_ft | t_ft),2);
mf_fluo_vec = nanmean(mf_array(:,b_ft | t_ft),2);
fluo_surge_size_vec = nanmean(fluo_array(:,t_ft),2) - nanmean(fluo_array(:,b_ft),2);

%% Look at distribution of surge event sizes
bins = linspace(-2,3);
surgeModel = fitgmdist(pt_surge_size_vec,2);

mu1 = surgeModel.mu(1);
sig1 = surgeModel.Sigma(1);
cp1 = surgeModel.ComponentProportion(1);

mu2 = surgeModel.mu(2);
sig2 = surgeModel.Sigma(2);
cp2 = surgeModel.ComponentProportion(2);

g1 = exp(-.5*((mu1-bins)./sig1).^2);
g1 = g1/sum(g1)*cp1;
g2 = exp(-.5*((mu2-bins)./sig2).^2);
g2 = g2/sum(g2)*cp2;

surge_hist = figure;
hold on
histogram(pt_surge_size_vec,bins,'Normalization','Probability')
xlabel('surge size (au)')
ylabel('share')
set(gca,'Fontsize',14)
saveas(surge_hist,[FigPath 'surge_size_hist.png'])
% plot(bins,g1)
% plot(bins,g2)