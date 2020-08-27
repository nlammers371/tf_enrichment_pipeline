% Script to attempt a systematic dissection of various factors driving
% proteinxtranscription burst coincidence
clear
close all
addpath('../utilities')
% define core ID variables
project = 'Dl-Ven_snaBAC-mCh';
% project = 'Dl-Ven_hbP2P-mCh';
DropboxFolder =  'E:\Nick\LivemRNA\Dropbox (Personal)\';
[~, DataPath, FigRoot] =   header_function(DropboxFolder, project);
FigPath = [FigRoot '\' project '\burst_snip_movie\'];
mkdir(FigPath)
% load data
load([DataPath 'hmm_input_output_results.mat'])
load([DataPath 'nucleus_struct_protein.mat'])
w = 7;
K = 3;
load([DataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'],'hmm_input_output');

% define size of window of interest
roi_window = 6; 
window_size = 15;
start = window_size + 2;
% extract roi_vectors from wapo and locus arrays
locus_protein_vec = nansum(results_struct.spot_array_dt(:,start:start + roi_window),2);% - ...
%     nanmean(results_struct.spot_array_dt(:,start-2-roi_window:start-2),2);
% swap_protein_vec = nanmean(results_struct.swap_array_dt(:,start:start + roi_window),2);
mf_protein_vec = nanmean(results_struct.mf_array(:,start:start + roi_window),2);
% pull other trend vectors
feature_sign_vec = results_struct.feature_sign_vec';
lag_size_vec = results_struct.lag_size_vec';
lead_size_vec = results_struct.lead_size_vec';
lag_dur_vec = results_struct.lag_dur_vec';
lead_dur_vec = results_struct.lead_dur_vec';
tr_burst_size_vec = lag_dur_vec.*lag_size_vec;
%%
% set basic analyisis parameters
min_pause_len = 6; % minimum length of preceding OFF period (in time steps)
min_burst_short = 2;
min_burst_long = 5;
long_flag = true;
suffix = 'long';
if ~long_flag
    suffix = 'short';
end

if long_flag
    % generate basic filter for target locus and computational controls
    burst_ft = feature_sign_vec == 1&lead_dur_vec>=min_pause_len&lag_dur_vec>=min_burst_long;%&mf_protein_vec>175; % filter for rise events
else
    burst_ft = feature_sign_vec == 1&lead_dur_vec>=min_pause_len&lag_dur_vec<min_burst_long & lag_dur_vec>=2;
end

% get raw protein snips for particles and time points of interest
nc_pt_id_index = [nucleus_struct_protein.ParticleID];
burst_rise_pt_ids = results_struct.particle_id_vec(burst_ft);
burst_rise_times = results_struct.center_time_vec(burst_ft);

snip_sz = size(nucleus_struct_protein(1).spot_protein_snips,1);
% initialize array to store snips 
snip_array = NaN(snip_sz, snip_sz, numel(burst_rise_pt_ids), 2*window_size + 1);

% iterate through events and pull snips
for i = 1:numel(burst_rise_pt_ids)
    ind = find(nc_pt_id_index==burst_rise_pt_ids(i));
    c_time = burst_rise_times(i);
    time_vec = nucleus_struct_protein(ind).time;
    % find apprixmate center time
    [~, c_ind] = min(abs(time_vec-c_time));
    full_range = c_ind-window_size:c_ind+window_size;
    trunc_range = full_range(full_range>0 & full_range <= numel(time_vec));
    % extract samples
    snip_array(:,:,i,ismember(full_range,trunc_range)) = nucleus_struct_protein(ind).spot_protein_snips(:,:,trunc_range);
end
    
% generate de-trended versions 
close all
window_vec = -window_size:window_size;
mean_snip_array = nanmean(snip_array,3);
trend_vec = NaN(size(window_vec));
dist_array = zeros(snip_sz,snip_sz);
dist_array(ceil(snip_sz/2),ceil(snip_sz/2)) = 1;
dist_array = bwdist(dist_array);
for i = 1:numel(trend_vec)
    snip = mean_snip_array(:,:,1,i);
    trend_vec(i) = nanmean(snip(dist_array <= 2));
end

mdl = fitlm(window_vec',trend_vec');
offset_vec = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2)*window_vec;
mean_snip_array_dt = mean_snip_array;
for i = 1:numel(window_vec)
    mean_snip_array_dt(:,:,1,i) = mean_snip_array_dt(:,:,1,i) - offset_vec(i);
end

time_vec = window_vec * 20;
%%% adjust contrast settings
hm_cm = flipud(brewermap([],'RdYlBu'));
hm_mid = hm_cm(129,:);
hm_max = hm_cm(256,:);
hm_trunc = interp1([0,1]',vertcat(hm_mid,hm_max),linspace(0,1,56)');

hm_new = vertcat(hm_cm(1:128,:),hm_trunc,repmat(hm_trunc(end,:),128-56,1));

close all
for w = 9:window_size+10
    snip_fig = figure;
    colormap(hm_new)
    imagesc(mean_snip_array_dt(:,:,1,w));
    caxis([-.5 .5]);
    h = colorbar;
    ylabel(h,'Dorsal Concentration')
    title([num2str(time_vec(w)) 's']);
    saveas(snip_fig,[FigPath 'snip_frame_' suffix '_' sprintf('%02d', w) '.tif'])    
end


    
