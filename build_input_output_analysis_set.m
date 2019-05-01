% Script to build data set for systematic input/output analyses
clear 
close all
% define ID variables
K = 3;
w = 6;
project = 'Dl_Venus_snaBAC_MCPmCherry_Leica_Zoom2_7uW14uW';
% dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
dropboxFolder = 'C:\Users\nlamm\Dropbox (Garcia Lab)\';
dataPath = [dropboxFolder '\ProcessedEnrichmentData\' project '\'];
% load data set
load([dataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'])
% set size of time series "reads"
window_size_time = 5*60;
Tres = 20;%hmm_input_output(1).Tres;
window_size = round(window_size_time/Tres);
pt_sm_kernel = 1; % size of kernel used to smooth protein channel
%%%
min_change_size_mcp = prctile(vertcat(hmm_input_output.r_vec),35); % eyeballed
min_peak_prom_mcp = prctile(vertcat(hmm_input_output.r_vec),10);
pt_vec = [hmm_input_output.spot_protein];
min_peak_prom_pt = prctile(pt_vec,10) - prctile(pt_vec,5);
%%
base_in_out_struct = [];
for i = 1:numel(hmm_input_output)
    r_vec = hmm_input_output(i).r_vec;
    if numel(r_vec) < 2*window_size+1
        continue
    end
    temp = hmm_input_output(i); % inistialize temporary structure
    % find and document features in transcriptional output
    mcp_peak_flags = false(size(r_vec));
    mcp_trough_flags = false(size(r_vec));
    mcp_feature_width = NaN(size(r_vec));
    mcp_feature_prom = NaN(size(r_vec));
    mcp_feature_sep = NaN(size(r_vec));
    % find peaks and troughs
    [~,peak_loc,peak_width,peak_prom] = findpeaks(r_vec,'MinPeakProminence',min_peak_prom_mcp);
    peak_sep = nanmean(vertcat([NaN diff(peak_loc')],[diff(peak_loc') NaN]));
    [~,trough_loc,trough_width,trough_prom] = findpeaks(-r_vec,'MinPeakProminence',min_peak_prom_mcp);
    trough_sep = nanmean(vertcat([NaN diff(trough_loc')],[diff(trough_loc') NaN]));
    % record    
	mcp_peak_flags(peak_loc) = true;
    mcp_trough_flags(trough_loc) = true;
    indices = [peak_loc' trough_loc'];
    mcp_feature_width(indices) = [peak_width' trough_width'];
    mcp_feature_prom(indices) = [peak_prom' trough_prom'];
    mcp_feature_sep(indices) = [peak_sep trough_sep]; 
    % find changepoints
    mcp_change_flags = false(size(r_vec));    
    mcp_change_size = NaN(size(r_vec));
    mcp_change_sep = NaN(size(r_vec));
    mcp_change_points = findchangepts(r_vec,'MinThreshold',min_change_size_mcp);
    mcp_change_sep(mcp_change_points) = nanmean(vertcat([NaN diff(mcp_change_points')],[diff(mcp_change_points') NaN]));
    change_sign = sign(vertcat(0,diff(r_vec)));
    change_size = abs(vertcat(0,diff(r_vec)));
    % record
    mcp_change_flags(mcp_change_points) = change_sign(mcp_change_points);
    mcp_change_size(mcp_change_points) = change_size(mcp_change_points);
    %%% store in temp
    temp.mcp_feature_sep = mcp_feature_sep;
    temp.mcp_peak_flags = mcp_peak_flags;
    temp.mcp_trough_flags = mcp_trough_flags;
    temp.mcp_feature_width = mcp_feature_width;
    temp.mcp_feature_prom = mcp_feature_prom;
    temp.mcp_change_sep = mcp_change_sep;
    temp.mcp_change_flags = mcp_change_flags;
    temp.mcp_change_size = mcp_change_size;
    
    %%% Now protein features
    p_vec = hmm_input_output(i).spot_protein_all;
    mf_vec = hmm_input_output(i).mf_protein_all;
    pt_vec_sm = imgaussfilt(p_vec-mf_vec,pt_sm_kernel);
    % find peaks and troughs
    pt_peak_flags = false(size(r_vec));
    pt_trough_flags = false(size(r_vec));
    pt_feature_width = NaN(size(r_vec));
    pt_feature_prom = NaN(size(r_vec));
    pt_feature_sep = NaN(size(r_vec));
    [~,peak_loc,peak_width,peak_prom] = findpeaks(pt_vec_sm,'MinPeakProminence',min_peak_prom_pt);    
    [~,trough_loc,trough_width,trough_prom] = findpeaks(-pt_vec_sm,'MinPeakProminence',min_peak_prom_pt);
    % record    
    pt_peak_flags(peak_loc) = true;
    pt_trough_flags(trough_loc) = true;
    indices = [peak_loc trough_loc];
    pt_feature_width(indices) = [peak_width trough_width];
    pt_feature_prom(indices) = [peak_prom trough_prom];
    % store in temp
    temp.pt_peak_flags = pt_peak_flags;
    temp.pt_trough_flags = pt_trough_flags;
    temp.pt_feature_width = pt_feature_width;
    temp.pt_feature_prom = pt_feature_prom;
    temp.delta_protein_smooth = pt_vec_sm;
    temp.n_reads = numel(r_vec) - 2*window_size;
    % add
    base_in_out_struct = [base_in_out_struct temp];
end
%% generate data set of reads
% input_output_snips = struct;
iter = sum([base_in_out_struct.n_reads]);
% base_var_cell = {'ParticleID','mf_protein','spot_protein'
tic
for i = 1:numel(base_in_out_struct)  
    % extract set of vectors
    t_vec = base_in_out_struct(i).time;
    r_vec = base_in_out_struct(i).r_vec;
    p_vec = base_in_out_struct(i).spot_protein;
    mf_vec = base_in_out_struct(i).mf_protein;
    serial_vec = base_in_out_struct(i).serial_protein;
    fluo_vec = base_in_out_struct(i).fluo;
    delta_protein_smooth = base_in_out_struct(i).delta_protein_smooth;        
    % feature vectors
    pt_peak_flags = base_in_out_struct(i).pt_peak_flags;
    pt_trough_flags = base_in_out_struct(i).pt_trough_flags;
    pt_feature_width = base_in_out_struct(i).pt_feature_width;
    pt_feature_prom = base_in_out_struct(i).pt_feature_prom;
    mcp_feature_sep = base_in_out_struct(i).mcp_feature_sep;
    mcp_peak_flags = base_in_out_struct(i).mcp_peak_flags;
    mcp_trough_flags = base_in_out_struct(i).mcp_trough_flags;
    mcp_feature_width = base_in_out_struct(i).mcp_feature_width;
    mcp_feature_prom = base_in_out_struct(i).mcp_feature_prom;
    mcp_change_sep = base_in_out_struct(i).mcp_change_sep;
    mcp_change_flags = base_in_out_struct(i).mcp_change_flags;
    mcp_change_size = base_in_out_struct(i).mcp_change_size;
    % iterate through single reads
    n_reads = base_in_out_struct(i).n_reads;
    for j = 1:n_reads        
        input_output_snips(iter).ParticleID = base_in_out_struct(i).ParticleID;
        % record data vectors
        input_output_snips(iter).time_vec = t_vec(j:j+2*window_size);
        input_output_snips(iter).hmm_vec = r_vec(j:j+2*window_size)';
        input_output_snips(iter).fluo_vec = fluo_vec(j:j+2*window_size);
        input_output_snips(iter).serial_protein_vec = serial_vec(j:j+2*window_size);
        input_output_snips(iter).mf_protein_vec = mf_vec(j:j+2*window_size);
        input_output_snips(iter).spot_protein_vec = p_vec(j:j+2*window_size);
        input_output_snips(iter).delta_protein_vec = delta_protein_smooth(j:j+2*window_size);
        % pull feature info
        input_output_snips(iter).t_center = t_vec(j+window_size);
        input_output_snips(iter).pt_peak_flag = pt_peak_flags(j+window_size);
        input_output_snips(iter).pt_trough_flag = pt_trough_flags(j+window_size);
        input_output_snips(iter).pt_feature_width = pt_feature_width(j+window_size);
        input_output_snips(iter).pt_feature_prom = pt_feature_prom(j+window_size);
        input_output_snips(iter).mcp_feature_sep = mcp_feature_sep(j+window_size);
        input_output_snips(iter).mcp_peak_flag = mcp_peak_flags(j+window_size);
        input_output_snips(iter).mcp_trough_flag = mcp_trough_flags(j+window_size);
        input_output_snips(iter).mcp_feature_width = mcp_feature_width(j+window_size);
        input_output_snips(iter).mcp_feature_prom = mcp_feature_prom(j+window_size);
        input_output_snips(iter).mcp_change_sep = mcp_change_sep(j+window_size);
        input_output_snips(iter).mcp_change_flags = mcp_change_flags(j+window_size);
        input_output_snips(iter).mcp_change_size = mcp_change_size(j+window_size);
        % calculate cross correlation
        spot_delta = input_output_snips(iter).spot_protein_vec - input_output_snips(iter).mf_protein_vec;
        ctrl_delta = input_output_snips(iter).serial_protein_vec - input_output_snips(iter).mf_protein_vec;
        input_output_snips(iter).spot_xcov = xcov(input_output_snips(iter).hmm_vec,...
            spot_delta,window_size,'Coeff');
        input_output_snips(iter).ctrl_xcov = xcov(input_output_snips(iter).hmm_vec,...
            ctrl_delta,window_size,'Coeff');
        iter = iter - 1;
    end    
end
toc