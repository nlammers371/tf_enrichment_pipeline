% script to check validity of enrichment vs. burst duration trend using
% stochastic simulations
clear
close all
% define core ID variables
project = 'Dl-Ven_snaBAC-mCh';
dropboxFolder =  'E:\Nick\LivemRNA\Dropbox\';
dataPath = [dropboxFolder 'ProcessedEnrichmentData\' project '/'];
figPath = [dropboxFolder 'LocalEnrichmentFigures\' project '\simulation_checks\'];
mkdir(figPath)
% HMM params
w = 7;
K = 3;  
alphaFrac = 1302 / 6000;
alpha = alphaFrac*w;
Tres = 20;
seq_length = 60 / Tres * 60;
n_traces = 250;
% protein sim params
protein_burst_dur = 60*2; % in seconds
window_size = 15;
% Set write path (inference results are now written to external directory)
hmm_suffix =  ['hmm_inference/w' num2str(w) '_K' num2str(K) '/']; 
file_list = dir([dataPath hmm_suffix 'hmm_results*.mat']);
if numel(file_list) > 1
    warning('multiple inference files detected. Ignoring all but first')
end

inference_results = load([dataPath hmm_suffix file_list(1).name]);
inference_results = inference_results.output;

% extract parameters
disp('extracting HMM parameters...')
[r_vec ,si] = sort(inference_results.r); 
A_mat = inference_results.A_mat(si,si);
v_vec = r_vec*Tres;
noise = sqrt(inference_results.noise);
pi0_vec = inference_results.pi0(si); 
eps = 1e-4;
R_raw = logm(A_mat) / Tres;
R = R_raw;
if ~isreal(R) || sum(R(:)<0) > K
    out = prob_to_rate_fit_sym(A_mat, Tres, 'gen', .005, 1);            
    R = out.R_out;     
end
% simulate transcription trajectories
disp('simualting transcription traces...')
master_struct = struct;
for i = 1:n_traces
    gillespie = synthetic_rate_gillespie(seq_length, alpha, ...
                                K, w, R, Tres, r_vec, noise, pi0_vec);
    fnames = fieldnames(gillespie);
    for f = 1:numel(fnames)
        master_struct(i).(fnames{f}) = gillespie.(fnames{f});
    end
end          
%%%
disp('simulating protein traces...')
t_vec = Tres:Tres:Tres*seq_length;
for i = 1:n_traces
    state_vec = master_struct(i).naive_states > 1;
    jump_vec = master_struct(i).transition_times;
    protein_vec_raw = zeros(size(jump_vec));    
    
    % find burst starts
    state_dd = [0 diff(state_vec)];
    burst_start_times = jump_vec(state_dd==1);
    
    % generate grid of protein levels
    protein_vec = zeros(size(t_vec));
    for t = 1:numel(t_vec)
        t_curr = t_vec(t);
        % find closest preceding rise
        last_id = find(burst_start_times < t_curr,1,'last');
        if ~isempty(last_id)
            protein_vec(t) = ((t_curr-burst_start_times(last_id))<=protein_burst_dur)*1.0;
        end
    end
    master_struct(i).protein_grid = protein_vec;
%     master_struct(i).state_grid = state_vec_interp;
end

% performing trace fits
disp('conducting single trace fits...')
A_log = log(inference_results.A_mat);
v = inference_results.r*Tres;
sigma = sqrt(inference_results.noise);
pi0_log = log(inference_results.pi0); 

fluo_values = cell(n_traces,1);
for i = 1:numel(master_struct)
    fluo = master_struct(i).fluo_MS2;
    fluo_values{i} = fluo;
end    
tic 
local_em_outputs = local_em_MS2_reduced_memory (fluo_values, ...
              v', sigma, pi0_log, A_log, K, w, alpha, 1, eps);
toc
soft_fit_struct = local_em_outputs.soft_struct;

%%% Convert soft fits to binary activity traces   
soft_state_cell = soft_fit_struct.p_z_log_soft;
for i = 1:n_traces
    % convert to binary trace
    ss_fit = exp(soft_state_cell{i});
    [~, hard_fit] = max(ss_fit);
    z_vec = hard_fit > 1;
    % record inferred trajectories
    master_struct(i).soft_vec_fit = ss_fit;
    master_struct(i).hard_vec_fit = hard_fit;
    master_struct(i).bin_vec_fit = z_vec;
    
    % identify burst characteristics       
    z_prob_vec = sum(ss_fit(2:3));
    zd_full = [0 diff(z_vec)];
    change_points = find(zd_full~=0);
    rise_points = find(zd_full>0);
    % calculate duration of bursts
    burst_dist_vec_lead = diff([NaN rise_points]);
    dur_vec_lag = diff([change_points NaN]);
    dur_vec_lead = diff([NaN change_points]);      
    % generate full-length vectors    
    z_dur_lag_vec_full = NaN(size(z_vec));
    z_dur_lead_vec_full = NaN(size(z_vec));
    z_dur_lag_vec_full(change_points) = dur_vec_lag;  
    z_dur_lead_vec_full(change_points) = dur_vec_lead;         

    % record    
    master_struct(i).z_dur_lag_vec = z_dur_lag_vec_full;  
    master_struct(i).z_dur_lead_vec = z_dur_lead_vec_full;    
    master_struct(i).change_points = change_points;
    master_struct(i).z_diff_vec = zd_full;
    master_struct(i).z_prob_vec = z_prob_vec';
    master_struct(i).time = t_vec;
end

%%% Now compile snips for average burst dynamics analyses
% generate master set of vectors for feature classification
time_vec = round([master_struct.time]/60);
z_dur_lag_list = [master_struct.z_dur_lag_vec];
z_dur_lead_list = [master_struct.z_dur_lead_vec];
z_diff_list = [master_struct.z_diff_vec];
% initialize results structure
results_struct = struct;   
% calculate expected number of features for pre-allocation
n_entries = sum(z_diff_list~=0);

% initialize arrays
hmm_array = NaN(n_entries,2*window_size+1);
spot_array = NaN(n_entries,2*window_size+1);
lag_dur_vec = NaN(1,n_entries);
lead_dur_vec = NaN(1,n_entries);
feature_sign_vec = NaN(1,n_entries);
iter = 1;
for j = 1:numel(master_struct)       
    time = master_struct(j).time;
    % feature classification vectors        
    z_dur_lag_vec_full = master_struct(j).z_dur_lag_vec;            
    z_dur_lead_vec_full = master_struct(j).z_dur_lead_vec;      
    z_diff_vec = master_struct(j).z_diff_vec;      
    hmm_vec = master_struct(i).soft_vec_fit;
    % protein 
    spot_protein = master_struct(j).protein_grid;        
    % find features
    id_list = find(z_diff_vec~=0);
    for id = id_list
        full_range = id - window_size:id+window_size;
        true_range = full_range(full_range>0&full_range<=numel(spot_protein));
        % record
        ft1 = ismember(full_range,true_range);
        % extract rde-trended fragments 
        spot_fragment = spot_protein(true_range);   
        % output and mf fragments        
        hmm_fragment = hmm_vec(true_range);        
        % protein snip               
        spot_array(iter,ft1) = spot_fragment;        
        % save hmm snip
        hmm_array(iter,ft1) = hmm_fragment; 
        % save other info
        lag_dur_vec(iter) = z_dur_lag_vec_full(id);
        lead_dur_vec(iter) = z_dur_lead_vec_full(id);   
        feature_sign_vec(iter) = sign(z_diff_vec(id));
        % increment
        iter = iter + 1; 
    end       
end     

%%% perform analysis
% (2) enrichment vs burst duration
nBoots = 100;
roi_window = 6;
start = window_size + 2;
min_prev_lag = 5;
% generate protein vec
spot_protein_vec = nanmean(spot_array(:,start:start + roi_window-1),2) - ...
    nanmean(spot_array(:,start-roi_window-1:start-2),2);
% filter
analysis_ft = feature_sign_vec==1&lead_dur_vec>min_prev_lag;
lag_dur_analysis = lag_dur_vec(analysis_ft);
spot_analysis = spot_protein_vec(analysis_ft);
% calculate bins
dur_bins = linspace(0,12,25);
dur_sigma = .5*median(diff(dur_bins));
% obtain average protein enrichment 
spot_dur_array = NaN(nBoots,numel(dur_bins));
sample_index = 1:sum(analysis_ft);
for n = 1:nBoots
    s_ids = randsample(sample_index,sum(analysis_ft),true);    
    spot_boot = spot_analysis(s_ids);
    dur_boot = lag_dur_analysis(s_ids);
    for d = 1:numel(dur_bins)
        dur_val = dur_bins(d);
        dur_weights = exp(-.5*((dur_boot-dur_val)/dur_sigma).^2);
        spot_dur_array(n,d) = nansum(spot_boot.*dur_weights') / nansum(dur_weights);
    end
end

% calculate bootstrap mean and standard error
spot_dur_mean = nanmean(spot_dur_array);
spot_dur_ste = nanstd(spot_dur_array);

window_vec = -window_size:window_size;
% make figure showing average profile
profile_fig = figure;
hold on
plot(window_vec*20/60,nanmean(spot_array(analysis_ft,:))-...
    nanmean(nanmean((spot_array(analysis_ft,:)))),'LineWidth',2)
plot(window_vec*20/60,nanmean(spot_array(feature_sign_vec==1&~analysis_ft,:))-...
    nanmean(nanmean(spot_array(feature_sign_vec==1&~analysis_ft,:))),'LineWidth',2)
grid on
xlabel('distance from start of burst (min)')
ylabel('protein level (au)')
set(gca,'Fontsize',12);
legend('lag-filtered','all bursts')
saveas(profile_fig, [figPath 'sim_protein_profiles.tif'])
saveas(profile_fig, [figPath 'sim_protein_profiles.pdf'])

% duration correlation
dur_fig = figure;
hold on
e = errorbar(dur_bins*20/60,spot_dur_mean,spot_dur_ste,'Color','black');
scatter(dur_bins*20/60,spot_dur_mean,'filled','MarkerEdgeColor','black')
e.CapSize = 0;
grid on
xlabel('transcription burst duration (min)')
ylabel('protein level (au)')
set(gca,'Fontsize',12);
saveas(profile_fig, [figPath 'sim_dur_trend.tif'])
saveas(profile_fig, [figPath 'sim_dur_trend.pdf'])