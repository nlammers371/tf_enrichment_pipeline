% Script to build data set for systematic input/output analyses
% DESCRIPTION
% Script to conduct HMM inference
%
% ARGUMENTS
% project: master ID variable 
%
% wInf: memory used for inference
%
% KInf: number of states used for inference
%
% OPTIONS
% dropboxFolder: Path to data folder where you wish to save
%                pipeline-generated data sets and figures. If this
%                var is not specified, output will be saved one level
%                above git repo in the folder structure
%
% controlProject: specifies a project to use as an external control
%
% OUTPUT: hmm_input_output, structure containing vectors of protein and MS2
% intensities, along with corresponding HMM-decoded activity trajectories

function input_output_snips = main07_build_input_output_analysis_set(project,varargin)

close all
dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
dataPath = [dropboxFolder 'ProcessedEnrichmentData\' project '/'];
K = 3;
w = 7;

%%%%%%%%%%%%%%
for i = 1:numel(varargin)    
    if strcmpi(varargin{i},'dropboxFolder')
        dataRoot = [varargin{i+1} 'ProcessedEnrichmentData\'];
    end
    if ischar(varargin{i}) && i ~= numel(varargin)
        if ismember(varargin{i},{'dpBootstrap','controlProject'})
            eval([varargin{i} '=varargin{i+1};']);
        end
    end
end
disp('building "final" input-output analysis set...')
tic
% load data set
load([dataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'])
% set size of time series "reads"
window_size_time = 5*60;
Tres = 20;%hmm_input_output(1).Tres;
window_size = round(window_size_time/Tres);
pt_sm_kernel = 1; % size of kernel used to smooth protein channel
%%%
% set scales for feature identification
fluo_scale = prctile([hmm_input_output.fluo],40);
hmm_scale = prctile([hmm_input_output.r_vec],40);

pt_vec = [hmm_input_output.spot_protein];
min_peak_prom_pt = prctile(pt_vec,10) - prctile(pt_vec,5);

n_vec = NaN(size(hmm_input_output));
for i = 1:numel(hmm_input_output)
    n_vec(i) = numel(hmm_input_output(i).time);
end
n_ids = find(n_vec >=2*window_size+1);
iter = numel(n_ids);
% base_in_out_struct = [];
for i = n_ids
    r_vec = hmm_input_output(i).r_vec;
    temp = hmm_input_output(i); % inistialize temporary structure
    % find and document features in hmm-decoded transcriptional output
    hmm_peak_flags = false(size(r_vec));
    hmm_trough_flags = false(size(r_vec));
    hmm_feature_width = NaN(size(r_vec));
    hmm_feature_prom = NaN(size(r_vec));
    hmm_feature_sep = NaN(size(r_vec));
    % find peaks and troughs                
    [~,peak_loc,peak_width,peak_prom] = findpeaks(r_vec,'MinPeakProminence',hmm_scale);
    [~,trough_loc,trough_width,trough_prom] = findpeaks(-r_vec,'MinPeakProminence',hmm_scale);    
    % record    
	hmm_peak_flags(peak_loc) = true;
    hmm_trough_flags(trough_loc) = true;
    indices = [peak_loc trough_loc];
    hmm_feature_width(indices) = [peak_width trough_width];
    hmm_feature_prom(indices) = [peak_prom trough_prom];    
    % find changepoints
    hmm_change_flags = zeros(size(r_vec));    
    hmm_change_size = NaN(size(r_vec));
    hmm_change_sep = NaN(size(r_vec));
    hmm_change_points = findchangepts(r_vec,'MinThreshold',hmm_scale);
    if numel(hmm_change_points) > 1        
        hmm_change_sep(hmm_change_points) = nanmean(vertcat([NaN diff(hmm_change_points)],[diff(hmm_change_points) NaN]));
    end
    change_sign = sign([0,diff(r_vec)]);
    change_size = abs([0,diff(r_vec)]);
    % record
    hmm_change_flags(hmm_change_points) = change_sign(hmm_change_points);
    hmm_change_size(hmm_change_points) = change_size(hmm_change_points);
    %%% store in temp
    temp.hmm_feature_sep = hmm_feature_sep;
    temp.hmm_peak_flags = hmm_peak_flags;
    temp.hmm_trough_flags = hmm_trough_flags;
    temp.hmm_feature_width = hmm_feature_width;
    temp.hmm_feature_prom = hmm_feature_prom;
    temp.hmm_change_sep = hmm_change_sep;
    temp.hmm_change_flags = hmm_change_flags;
    temp.hmm_change_size = hmm_change_size;
    
    % find and document features in raw transcriptional output
    fluo_vec = hmm_input_output(i).fluo;
    fluo_peak_flags = false(size(r_vec));
    fluo_trough_flags = false(size(r_vec));
    fluo_feature_width = NaN(size(r_vec));
    fluo_feature_prom = NaN(size(r_vec));
    fluo_feature_sep = NaN(size(r_vec));
    % find peaks and troughs
    [~,peak_loc,peak_width,peak_prom] = findpeaks(fluo_vec,'MinPeakProminence',fluo_scale);
    [~,trough_loc,trough_width,trough_prom] = findpeaks(-fluo_vec,'MinPeakProminence',fluo_scale);    
    % record    
	fluo_peak_flags(peak_loc) = true;
    fluo_trough_flags(trough_loc) = true;
    indices = [peak_loc trough_loc];
    fluo_feature_width(indices) = [peak_width trough_width];
    fluo_feature_prom(indices) = [peak_prom trough_prom];
    % find changepoints
    fluo_change_flags = zeros(size(r_vec));    
    fluo_change_size = NaN(size(r_vec));
    fluo_change_sep = NaN(size(r_vec));
    fluo_change_points = findchangepts(fluo_vec,'MinThreshold',fluo_scale);
    if numel(fluo_change_points) > 1
        fluo_change_sep(fluo_change_points) = nanmean(vertcat([NaN diff(fluo_change_points)],[diff(fluo_change_points) NaN]));
    end
    change_sign = sign([0,diff(fluo_vec)]);
    change_size = abs([0,diff(fluo_vec)]);
    % record
    fluo_change_flags(fluo_change_points) = change_sign(fluo_change_points);
    fluo_change_size(fluo_change_points) = change_size(fluo_change_points);
    %%% store in temp
    temp.fluo_feature_sep = fluo_feature_sep;
    temp.fluo_peak_flags = fluo_peak_flags;
    temp.fluo_trough_flags = fluo_trough_flags;
    temp.fluo_feature_width = fluo_feature_width;
    temp.fluo_feature_prom = fluo_feature_prom;
    temp.fluo_change_sep = fluo_change_sep;
    temp.fluo_change_flags = fluo_change_flags;
    temp.fluo_change_size = fluo_change_size;
    %%% Now protein features
    p_vec = hmm_input_output(i).spot_protein;
    mf_vec = hmm_input_output(i).mf_protein;
    pt_vec_sm = imgaussfilt(p_vec-mf_vec,pt_sm_kernel);
    % find peaks and troughs
    pt_peak_flags = false(size(r_vec));
    pt_trough_flags = false(size(r_vec));
    pt_feature_width = NaN(size(r_vec));
    pt_feature_prom = NaN(size(r_vec));    
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
    base_in_out_struct(iter) = temp;
    iter = iter - 1;
end

%%% generate data set of reads
% input_output_snips = struct;
iter = sum([base_in_out_struct.n_reads]);
% base_var_cell = {'ParticleID','mf_protein','spot_protein'

for i = 1:numel(base_in_out_struct)  
    % extract set of vectors
    t_vec = base_in_out_struct(i).time;
    r_vec = base_in_out_struct(i).r_vec;
    p_vec = base_in_out_struct(i).spot_protein;
    mf_vec = base_in_out_struct(i).mf_protein;
    serial_vec = base_in_out_struct(i).serial_protein;
    swap_spot_protein_vec = base_in_out_struct(i).swap_spot_protein;
    swap_serial_protein_vec = base_in_out_struct(i).swap_serial_protein;
    swap_mf_vec = base_in_out_struct(i).swap_mf_protein;
    fluo_vec = base_in_out_struct(i).fluo;
    swap_fluo_vec = base_in_out_struct(i).swap_fluo;
    swap_hmm_vec = base_in_out_struct(i).swap_hmm;
    delta_protein_smooth = base_in_out_struct(i).delta_protein_smooth;        
    % feature vectors
    pt_peak_flags = base_in_out_struct(i).pt_peak_flags;
    pt_trough_flags = base_in_out_struct(i).pt_trough_flags;
    pt_feature_width = base_in_out_struct(i).pt_feature_width;
    pt_feature_prom = base_in_out_struct(i).pt_feature_prom;
    
    hmm_feature_sep = base_in_out_struct(i).hmm_feature_sep;
    hmm_peak_flags = base_in_out_struct(i).hmm_peak_flags;
    hmm_trough_flags = base_in_out_struct(i).hmm_trough_flags;
    hmm_feature_width = base_in_out_struct(i).hmm_feature_width;
    hmm_feature_prom = base_in_out_struct(i).hmm_feature_prom;
    hmm_change_sep = base_in_out_struct(i).hmm_change_sep;
    hmm_change_flags = base_in_out_struct(i).hmm_change_flags;
    hmm_change_size = base_in_out_struct(i).hmm_change_size;
    
    fluo_feature_sep = base_in_out_struct(i).fluo_feature_sep;
    fluo_peak_flags = base_in_out_struct(i).fluo_peak_flags;
    fluo_trough_flags = base_in_out_struct(i).fluo_trough_flags;
    fluo_feature_width = base_in_out_struct(i).fluo_feature_width;
    fluo_feature_prom = base_in_out_struct(i).fluo_feature_prom;
    fluo_change_sep = base_in_out_struct(i).fluo_change_sep;
    fluo_change_flags = base_in_out_struct(i).fluo_change_flags;
    fluo_change_size = base_in_out_struct(i).fluo_change_size;
    % iterate through single reads
    n_reads = base_in_out_struct(i).n_reads;
    for j = 1:n_reads        
        input_output_snips(iter).ParticleID = base_in_out_struct(i).ParticleID;
        % record data vectors
        input_output_snips(iter).time_vec = t_vec(j:j+2*window_size);
        input_output_snips(iter).hmm_vec = r_vec(j:j+2*window_size);
        input_output_snips(iter).fluo_vec = fluo_vec(j:j+2*window_size);
        input_output_snips(iter).swap_hmm_vec = swap_hmm_vec(j:j+2*window_size);
        input_output_snips(iter).swap_fluo_vec = swap_fluo_vec(j:j+2*window_size);
        input_output_snips(iter).serial_protein_vec = serial_vec(j:j+2*window_size);
        input_output_snips(iter).mf_protein_vec = mf_vec(j:j+2*window_size);
        input_output_snips(iter).swap_mf_protein_vec = swap_mf_vec(j:j+2*window_size);
        input_output_snips(iter).spot_protein_vec = p_vec(j:j+2*window_size);
        input_output_snips(iter).swap_spot_protein_vec = swap_spot_protein_vec(j:j+2*window_size);
        input_output_snips(iter).swap_serial_protein_vec = swap_serial_protein_vec(j:j+2*window_size);
        input_output_snips(iter).delta_protein_vec = delta_protein_smooth(j:j+2*window_size);
        % pull feature info
        input_output_snips(iter).t_center = t_vec(j+window_size);
        input_output_snips(iter).pt_peak_flag = pt_peak_flags(j+window_size);
        input_output_snips(iter).pt_trough_flag = pt_trough_flags(j+window_size);
        input_output_snips(iter).pt_feature_width = pt_feature_width(j+window_size);
        input_output_snips(iter).pt_feature_prom = pt_feature_prom(j+window_size);
        % hmm
        input_output_snips(iter).hmm_feature_sep = hmm_feature_sep(j+window_size);
        input_output_snips(iter).hmm_peak_flag = hmm_peak_flags(j+window_size);
        input_output_snips(iter).hmm_trough_flag = hmm_trough_flags(j+window_size);
        input_output_snips(iter).hmm_feature_width = hmm_feature_width(j+window_size);
        input_output_snips(iter).hmm_feature_prom = hmm_feature_prom(j+window_size);
        input_output_snips(iter).hmm_change_sep = hmm_change_sep(j+window_size);
        input_output_snips(iter).hmm_change_flags = hmm_change_flags(j+window_size);
        input_output_snips(iter).hmm_change_size = hmm_change_size(j+window_size);
        % fluo
        input_output_snips(iter).fluo_feature_sep = fluo_feature_sep(j+window_size);
        input_output_snips(iter).fluo_peak_flag = fluo_peak_flags(j+window_size);
        input_output_snips(iter).fluo_trough_flag = fluo_trough_flags(j+window_size);
        input_output_snips(iter).fluo_feature_width = fluo_feature_width(j+window_size);
        input_output_snips(iter).fluo_feature_prom = fluo_feature_prom(j+window_size);
        input_output_snips(iter).fluo_change_sep = fluo_change_sep(j+window_size);
        input_output_snips(iter).fluo_change_flags = fluo_change_flags(j+window_size);
        input_output_snips(iter).fluo_change_size = fluo_change_size(j+window_size);
        % calculate cross correlation
        spot_delta = input_output_snips(iter).spot_protein_vec - input_output_snips(iter).mf_protein_vec;
        ctrl_delta = input_output_snips(iter).serial_protein_vec - input_output_snips(iter).mf_protein_vec;
        swap_delta = input_output_snips(iter).swap_spot_protein_vec - input_output_snips(iter).swap_mf_protein_vec;
        % fluo cross-corr
        input_output_snips(iter).fluo_spot_xcov = xcov(input_output_snips(iter).fluo_vec,...
            spot_delta,window_size,'Coeff');
        input_output_snips(iter).fluo_ctrl_xcov = xcov(input_output_snips(iter).fluo_vec,...
            ctrl_delta,window_size,'Coeff');
        input_output_snips(iter).fluo_swap_xcov = xcov(input_output_snips(iter).fluo_vec,...
            swap_delta,window_size,'Coeff');
        % hmm cross-corr
        input_output_snips(iter).hmm_spot_xcov = xcov(input_output_snips(iter).hmm_vec,...
            spot_delta,window_size,'Coeff');
        input_output_snips(iter).hmm_ctrl_xcov = xcov(input_output_snips(iter).hmm_vec,...
            ctrl_delta,window_size,'Coeff');
        input_output_snips(iter).hmm_swap_xcov = xcov(input_output_snips(iter).hmm_vec,...
            swap_delta,window_size,'Coeff');
        % increment
        iter = iter - 1;
    end    
end
disp('done.')
toc
% save
save([dataPath 'input_output_snips.mat'],'input_output_snips');