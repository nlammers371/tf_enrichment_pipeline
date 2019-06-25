clear
close all

% define core ID variables
project = 'Dl-Ven_snaBAC-mCh';
% project = 'Dl-Ven_hbP2P-mCh';
dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
dataPath = [dropboxFolder 'ProcessedEnrichmentData\' project '/'];
w = 7;
K = 3;
% window analysis params
window_size = 15; 
% arrays for fitting linear offsets
window_vec = -window_size:window_size;
fit_array = [ones(numel(window_vec),1) window_vec'];
% load input-output data set
load([dataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'],'hmm_input_output')

% perform time averaging to get protein trends over time
gap_vec_full = [hmm_input_output.dt_filter_gap];
pt_vec_full = [hmm_input_output.spot_protein];
pt_vec_full = pt_vec_full(~gap_vec_full);
virtual_vec_full = [hmm_input_output.serial_protein];
virtual_vec_full = virtual_vec_full(~gap_vec_full);
time_vec_full = [hmm_input_output.time];
time_vec_full = time_vec_full(~gap_vec_full);

time_index = unique(time_vec_full);
spot_protein_trend = NaN(size(time_index));
virtual_protein_trend = NaN(size(time_index));
t_sigma = 120;
for t = 1:numel(time_index)
    dt_vec = exp(-(time_vec_full-time_index(t)).^2 / 2 / t_sigma^2);
    % locus
    p_mean_spot = nansum(pt_vec_full.*dt_vec) ./ nansum(dt_vec);
    spot_protein_trend(t) = p_mean_spot;
    % virtual spot
    p_mean_virtual = nansum(virtual_vec_full.*dt_vec) ./ nansum(dt_vec);
    virtual_protein_trend(t) = p_mean_virtual;
end
options = optimoptions('lsqnonlin','Display','off');
n_features = 0;
% iterate
for i = 1:numel(hmm_input_output)    
    z_vec = hmm_input_output(i).z_vec' > 1;
    r_vec = hmm_input_output(i).r_vec';
    hmm_input_output(i).z_vec = z_vec;
    z_prob_vec = sum(hmm_input_output(i).z_mat(:,2:3),2);
    zd_full = [0 diff(z_vec)];
    change_points = find(zd_full~=0);
    rise_points = find(zd_full>0);
    % calculate duration of bursts
    burst_dist_vec_lead = diff([NaN rise_points]);
    dur_vec_lag = diff([change_points NaN]);
    dur_vec_lead = diff([NaN change_points]);
    % calculate intensity of bursts
    sz_vec_lag = NaN(size(dur_vec_lag));
    sz_vec_lead = NaN(size(dur_vec_lag));
    for j = 1:numel(change_points)-1
        r_mean = nanmean(r_vec(change_points(j):change_points(j+1)-1));
        sz_vec_lag(j) = r_mean;
        sz_vec_lead(j+1) = r_mean;
    end         
    % generate full-length vectors
    z_burst_dist_vec_full = NaN(size(z_vec));
    z_dur_lag_vec_full = NaN(size(z_vec));
    z_dur_lead_vec_full = NaN(size(z_vec));
    z_dur_lag_vec_full(change_points) = dur_vec_lag;  
    z_dur_lead_vec_full(change_points) = dur_vec_lead;  
    z_burst_dist_vec_full(rise_points) = burst_dist_vec_lead;
    sz_lag_vec_full = NaN(size(z_vec));
    sz_lead_vec_full = NaN(size(z_vec));
    sz_lag_vec_full(change_points) = sz_vec_lag;  
    sz_lead_vec_full(change_points) = sz_vec_lead; 
    
    % record
    hmm_input_output(i).z_burst_dist_vec = z_burst_dist_vec_full;  
    hmm_input_output(i).z_dur_lag_vec = z_dur_lag_vec_full;  
    hmm_input_output(i).z_dur_lead_vec = z_dur_lead_vec_full;
    hmm_input_output(i).sz_lag_vec = sz_lag_vec_full;  
    hmm_input_output(i).sz_lead_vec = sz_lead_vec_full;
    hmm_input_output(i).change_points = change_points;
    hmm_input_output(i).z_diff_vec = zd_full;
    hmm_input_output(i).z_prob_vec = z_prob_vec';
    % increment
    if ~isempty(change_points)
        n_features = n_features + numel(change_points)-1;
    end
    % generate detrended protein vectors
    dt_filter_gap = hmm_input_output(i).dt_filter_gap;
    % locus protein
    spot_pt_vec = hmm_input_output(i).spot_protein;
    nan_ft = ~isnan(spot_pt_vec)&~dt_filter_gap;
    time_vec = hmm_input_output(i).time;
    % get corresponding pt trend points
    spot_trend_vec = spot_protein_trend(ismember(time_index,time_vec));
    virtual_trend_vec = virtual_protein_trend(ismember(time_index,time_vec));
    % perform simple fit to acount for differences in offset and dynamic
    % range
    ob_fun = @(params) spot_pt_vec(nan_ft) - params(1)*spot_trend_vec(nan_ft) - params(2);
    fit = lsqnonlin(ob_fun,[1 0],[.2 5],[-300 300],options);    
    spot_pt_trend = fit(1)*spot_trend_vec(nan_ft) + fit(2);
    spot_protein_dt = NaN(size(dt_filter_gap));
    spot_protein_dt(nan_ft) = spot_pt_vec(nan_ft) - spot_pt_trend;
    hmm_input_output(i).spot_protein_dt = spot_protein_dt;
    
    % virtual spot protein
    virtual_pt_vec = hmm_input_output(i).serial_protein;    
    nan_ft = ~isnan(virtual_pt_vec)&~dt_filter_gap;
    ob_fun = @(params) virtual_pt_vec(nan_ft) - params(1)*virtual_trend_vec(nan_ft) - params(2);
    fit = lsqnonlin(ob_fun,[1 0],[.2 5],[-300 300],options);   
    virtual_pt_trend = fit(1)*virtual_trend_vec(nan_ft) + fit(2);
    serial_protein_dt = NaN(size(dt_filter_gap));
    serial_protein_dt(nan_ft) = virtual_pt_vec(nan_ft) - virtual_pt_trend;
    hmm_input_output(i).serial_protein_dt = serial_protein_dt;
    
    % swap spot protein
    swap_pt_vec = hmm_input_output(i).swap_spot_protein;  
    nan_ft = ~isnan(swap_pt_vec)&~dt_filter_gap;
    ob_fun = @(params) swap_pt_vec(nan_ft) - params(1)*spot_trend_vec(nan_ft) - params(2);
    fit = lsqnonlin(ob_fun,[1 0],[.2 5],[-300 300],options);   
    swap_pt_trend = fit(1)*spot_trend_vec(nan_ft) + fit(2);
    swap_spot_protein_dt = NaN(size(dt_filter_gap));
    swap_spot_protein_dt(nan_ft) = swap_pt_vec(nan_ft) - swap_pt_trend;
    hmm_input_output(i).swap_spot_protein_dt = swap_spot_protein_dt;
end

%%% initialize lists to store burst characteristics
% protein-weighted center of mass
pt_cm_spot_vec = NaN(n_features,1); 
pt_cm_swap_vec = NaN(n_features,1);
pt_cm_virtual_vec = NaN(n_features,1);
% net amount of protein under burst
pt_net_spot_vec = NaN(n_features,1); 
pt_net_swap_vec = NaN(n_features,1);
pt_net_virtual_vec = NaN(n_features,1);
% net amount of protein in first 2 minute of burst
pt_net2_spot_vec = NaN(n_features,1); 
pt_net2_swap_vec = NaN(n_features,1);
pt_net2_virtual_vec = NaN(n_features,1);
% average nucleus protein
pt_mf_vec = NaN(n_features,1);
% transcription burst intensity
tr_amp_vec = NaN(n_features,1);
tr_amp_prev_vec = NaN(n_features,1);
% tr burst duration
tr_dur_vec = NaN(n_features,1);
tr_dur_prev_vec = NaN(n_features,1);
burst_dist_prev_vec = NaN(n_features,1);
% tr burst class (feature or gap?)
tr_burst_class_vec = NaN(n_features,1);
%%% initialize ID vectors
feature_time_vec = NaN(n_features,1);
particle_id_vec = NaN(n_features,1);
fluo_vec = NaN(n_features,1);
% iterate through structure
iter = 1;                    
for j = 1:numel(hmm_input_output)
    % core ID vectors
    ParticleID = hmm_input_output(j).ParticleID;
    time = hmm_input_output(j).time;
    % feature classification vectors
    change_points = hmm_input_output(j).change_points;  
    gap_filter = hmm_input_output(j).dt_filter_gap;            
    z_dur_lag_vec_full = hmm_input_output(j).z_dur_lag_vec(change_points);            
    z_dur_lead_vec_full = hmm_input_output(j).z_dur_lead_vec(change_points);  
    sz_lag_vec_full = hmm_input_output(j).sz_lag_vec(change_points);            
    sz_lead_vec_full = hmm_input_output(j).sz_lead_vec(change_points);
    z_diff_vec = hmm_input_output(j).z_diff_vec(change_points);
    z_burst_dist_vec = hmm_input_output(j).z_burst_dist_vec(change_points);
    % activity
    fluo = hmm_input_output(j).fluo;
    r_vec = hmm_input_output(j).r_vec';
    % protein fields                        
    spot_protein = hmm_input_output(j).spot_protein_dt;
    mf_protein = hmm_input_output(j).mf_protein;
    swap_spot_protein = hmm_input_output(j).swap_spot_protein_dt;
    virtual_protein = hmm_input_output(j).serial_protein_dt;        
    % apply filter             
    spot_protein(gap_filter) = NaN;
    swap_spot_protein(gap_filter) = NaN;
    virtual_protein(gap_filter) = NaN;
    
    for i = 1:numel(change_points)-1        
        full_range = change_points(i):change_points(i+1)-1;       
        start_time = time(change_points(i));
        % record ID features
        particle_id_vec(iter) = ParticleID;
        fluo_vec(iter) = nanmean(fluo(full_range));
        feature_time_vec(iter) = nanmean(time(full_range));
        pt_mf_vec(iter) = nanmean(mf_protein(full_range));
        % burst class
        tr_burst_class_vec(iter) = z_diff_vec(i);
        % extract protein snips
        time_fragment = time(full_range);
        spot_fragment = spot_protein(full_range);
        swap_fragment = swap_spot_protein(full_range);
        virtual_fragment = virtual_protein(full_range);
        % calculate protein-weighted center of mass
        pt_cm_spot_vec(iter) = nansum(time_fragment.*spot_fragment) ./ nansum(spot_fragment) - start_time; 
        pt_cm_swap_vec(iter) = nansum(time_fragment.*swap_fragment) ./ nansum(swap_fragment) - start_time; 
        pt_cm_virtual_vec(iter) = nansum(time_fragment.*virtual_fragment) ./ nansum(virtual_fragment) - start_time; 
        % calculate net protein
        pt_net_spot_vec(iter) = nansum(spot_fragment);
        pt_net_swap_vec(iter) = nansum(swap_fragment);
        pt_net_virtual_vec(iter) = nansum(virtual_fragment);
        % calculate net protein in first 2 minutes
        early_filter = time_fragment - start_time <=60;
        pt_net2_spot_vec(iter) = nansum(spot_fragment(early_filter));
        pt_net2_swap_vec(iter) = nansum(swap_fragment(early_filter));
        pt_net2_virtual_vec(iter) = nansum(virtual_fragment(early_filter));
        % record transcription burst features
        tr_amp_vec(iter) = sz_lag_vec_full(i);
        tr_amp_prev_vec(iter) = sz_lead_vec_full(i);
        tr_dur_vec(iter) = z_dur_lag_vec_full(i);
        tr_dur_prev_vec(iter) = z_dur_lead_vec_full(i);
        burst_dist_prev_vec(iter) = z_burst_dist_vec(i);
        % increment
        iter = iter + 1;  
        if mod(iter,1000) == 0
            disp(iter)
        end
    end    
end        
% make results table
results_table = array2table([particle_id_vec, feature_time_vec, tr_burst_class_vec, fluo_vec, tr_amp_vec, tr_amp_prev_vec, tr_dur_vec, tr_dur_prev_vec,...
    pt_mf_vec, pt_cm_spot_vec, pt_cm_swap_vec, pt_cm_virtual_vec,...
    pt_net_spot_vec, pt_net_swap_vec, pt_net_virtual_vec, pt_net2_spot_vec, pt_net2_swap_vec, pt_net2_virtual_vec,burst_dist_prev_vec],...
    'VariableNames',{'ParticleID', 'time', 'burst_class', 'fluo', 'tr_amp', 'tr_amp_prev', 'tr_dur', 'tr_dur_prev', 'pt_mf', 'pt_cm_spot', 'pt_cm_swap', 'pt_cm_virtual',...
    'pt_net_spot', 'pt_net_swap', 'pt_net_virtual','pt_net2_spot', 'pt_net2_swap', 'pt_net2_virtual','prev_burst_dist'});

% save
save([dataPath 'hmm_burst_table.mat'],'results_table')

 
%%% Now compile snips for average burst dynamics analyses
% generate master set of vectors for feature classification
gap_filter_vec = [hmm_input_output.dt_filter_gap];
time_vec = round([hmm_input_output.time]/60);
z_dur_lag_list = [hmm_input_output.z_dur_lag_vec];
z_dur_lead_list = [hmm_input_output.z_dur_lead_vec];
sz_lag_list = [hmm_input_output.sz_lag_vec];
sz_lead_list = [hmm_input_output.sz_lead_vec];
z_diff_list = [hmm_input_output.z_diff_vec];
c = [hmm_input_output.z_diff_vec];
% initialize results structure
results_struct = struct;

   
% calculate expected number of features for pre-allocation
n_entries = sum(z_diff_list~=0&~gap_filter_vec);
% initialize arrays to store time series snip
fluo_array = NaN(n_entries,2*window_size+1);
hmm_array = NaN(n_entries,2*window_size+1);
spot_array = NaN(n_entries,2*window_size+1);
swap_array = NaN(n_entries,2*window_size+1);
virtual_array = NaN(n_entries,2*window_size+1);        
% initialize arrays to store id variables
particle_id_vec = NaN(1,n_entries);
center_time_vec = NaN(1,n_entries);
lag_dur_vec = NaN(1,n_entries);
lead_dur_vec = NaN(1,n_entries);
lag_size_vec = NaN(1,n_entries);
lead_size_vec = NaN(1,n_entries);
feature_sign_vec = NaN(1,n_entries);
prev_burst_dist_vec = NaN(1,n_entries);
% iterate through structure
iter = 1;                    
for j = 1:numel(hmm_input_output)
    % core ID vectors
    ParticleID = hmm_input_output(j).ParticleID;
    time = hmm_input_output(j).time;
    % feature classification vectors
    gap_filter = hmm_input_output(j).dt_filter_gap;      
    z_dur_lag_vec_full = hmm_input_output(j).z_dur_lag_vec;            
    z_dur_lead_vec_full = hmm_input_output(j).z_dur_lead_vec;  
    sz_lag_vec_full = hmm_input_output(j).sz_lag_vec;            
    sz_lead_vec_full = hmm_input_output(j).sz_lead_vec;
    z_diff_vec = hmm_input_output(j).z_diff_vec;  
    z_burst_dist_vec = hmm_input_output(j).z_burst_dist_vec;
    % activity
    fluo = hmm_input_output(j).fluo;
    r_vec = hmm_input_output(j).r_vec';
    % protein fields                        
    spot_protein = hmm_input_output(j).spot_protein_dt;
    swap_spot_protein = hmm_input_output(j).swap_spot_protein_dt;
    virtual_protein = hmm_input_output(j).serial_protein_dt;        
    % apply filter             
    spot_protein(gap_filter) = NaN;
    swap_spot_protein(gap_filter) = NaN;
    virtual_protein(gap_filter) = NaN;
    % find features
    id_list = find(z_diff_vec~=0&~gap_filter);
    for id = id_list
        full_range = id - window_size:id+window_size;
        true_range = full_range(full_range>0&full_range<=numel(virtual_protein));
        % record
        ft1 = ismember(full_range,true_range);
        % qc check
        if sum(~isnan(spot_protein(true_range))) >= window_size && sum(~isnan(swap_spot_protein(true_range)))...
                >= window_size && sum(~isnan(virtual_protein(true_range))) >= window_size
            % extract raw fragments
            spot_fragment = spot_protein(true_range);
            swap_fragment = swap_spot_protein(true_range);
            virtual_fragment = virtual_protein(true_range);
            fluo_fragment = fluo(true_range);
            hmm_fragment = r_vec(true_range);
            % fit linear offset          
            fit_sub_array = fit_array(ft1,:);
            fluo_fit = fit_sub_array(~isnan(fluo_fragment),:) \ fluo_fragment(~isnan(fluo_fragment))';
%           % save time snips               
            spot_array(iter,ft1) = spot_fragment;
            swap_array(iter,ft1) = swap_fragment;
            virtual_array(iter,ft1) = virtual_fragment;
            fluo_array(iter,ft1) = fluo_fragment - fluo_fit(1) - fluo_fit(2)*window_vec(ft1);
            hmm_array(iter,ft1) = hmm_fragment;
            % save other info
            particle_id_vec(iter) = ParticleID;
            center_time_vec(iter) = time(id);
            lag_dur_vec(iter) = z_dur_lag_vec_full(id);
            lead_dur_vec(iter) = z_dur_lead_vec_full(id);
            lag_size_vec(iter) = sz_lag_vec_full(id);
            lead_size_vec(iter) = sz_lead_vec_full(id);
            prev_burst_dist_vec(iter) = z_burst_dist_vec(id);
            feature_sign_vec(iter) = sign(z_diff_vec(id));
        end
        % increment
        iter = iter + 1;  
        if mod(iter,100) == 0
            disp(iter)
        end
    end    
end        
% record data
results_struct.spot_array = spot_array;
results_struct.swap_array = swap_array;
results_struct.virtual_array = virtual_array;
results_struct.fluo_array = fluo_array;
results_struct.hmm_array = hmm_array;      
particle_id_vec(iter) = ParticleID;
results_struct.center_time_vec = center_time_vec;
results_struct.lag_dur_vec = lag_dur_vec;
results_struct.lead_dur_vec = lead_dur_vec;
results_struct.lag_size_vec = lag_size_vec;
results_struct.lead_size_vec = lead_size_vec;
results_struct.feature_sign_vec = feature_sign_vec;
results_struct.prev_burst_dist = prev_burst_dist_vec;
% save
save([dataPath 'hmm_input_output_results.mat'],'results_struct')