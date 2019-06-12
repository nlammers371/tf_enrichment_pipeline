clear
close all

% define core ID variables
project = 'Dl-Ven_snaBAC-mCh';
% project = 'Dl-Ven_hbP2P-mCh';
dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
dataPath = [dropboxFolder 'ProcessedEnrichmentData\' project '/'];
w = 7;
K = 3;

% load input-output data set
load([dataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'],'hmm_input_output')
n_features = 0;
% iterate
for i = 1:numel(hmm_input_output)    
    z_vec = hmm_input_output(i).z_vec' > 1;
    r_vec = hmm_input_output(i).r_vec';
    hmm_input_output(i).z_vec = z_vec;
    z_prob_vec = sum(hmm_input_output(i).z_mat(:,2:3),2);
    zd = [0 diff(z_vec)];
    change_points = find(zd~=0);
    % calculate duration of bursts
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
    z_dur_lag_vec = NaN(size(z_vec));
    z_dur_lead_vec = NaN(size(z_vec));
    z_dur_lag_vec(change_points) = dur_vec_lag;  
    z_dur_lead_vec(change_points) = dur_vec_lead;  
    sz_lag_vec = NaN(size(z_vec));
    sz_lead_vec = NaN(size(z_vec));
    sz_lag_vec(change_points) = sz_vec_lag;  
    sz_lead_vec(change_points) = sz_vec_lead;  
    % record
    hmm_input_output(i).z_dur_lag_vec = dur_vec_lag;  
    hmm_input_output(i).z_dur_lead_vec = dur_vec_lead;
    hmm_input_output(i).sz_lag_vec = sz_vec_lag;  
    hmm_input_output(i).sz_lead_vec = sz_vec_lag;
    hmm_input_output(i).change_points = change_points;
    hmm_input_output(i).z_diff_vec = zd(change_points);
    hmm_input_output(i).z_prob_vec = z_prob_vec';
    % increment
    if ~isempty(change_points)
        n_features = n_features + numel(change_points)-1;
    end
    % generate detrended protein vectors
%     dt_filter_gap = hmm_input_output(i).dt_filter_gap;
    % locus protein
    spot_pt_vec = hmm_input_output(i).spot_protein;%(~dt_filter_gap);
    time_vec = hmm_input_output(i).time;%(~dt_filter_gap);
    p = polyfit(time_vec,spot_pt_vec,2);
    spot_pt_trend = polyval(p,time_vec);
    hmm_input_output(i).spot_protein_dt = spot_pt_vec - spot_pt_trend; 
    % virtual spot protein
    virtual_pt_vec = hmm_input_output(i).serial_protein;%(~dt_filter_gap);    
    p = polyfit(time_vec,virtual_pt_vec,2);
    virtual_pt_trend = polyval(p,time_vec);
    hmm_input_output(i).serial_protein_dt = virtual_pt_vec - virtual_pt_trend;
    % swap spot protein
    swap_pt_vec = hmm_input_output(i).swap_spot_protein;%(~dt_filter_gap);    
    p = polyfit(time_vec,swap_pt_vec,2);
    swap_pt_trend = polyval(p,time_vec);
    hmm_input_output(i).swap_spot_protein_dt = swap_pt_vec - swap_pt_trend;
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
    gap_filter = hmm_input_output(j).dt_filter_gap;      
    change_points = hmm_input_output(j).change_points;    
    z_dur_lag_vec = hmm_input_output(j).z_dur_lag_vec;            
    z_dur_lead_vec = hmm_input_output(j).z_dur_lead_vec;  
    sz_lag_vec = hmm_input_output(j).sz_lag_vec;            
    sz_lead_vec = hmm_input_output(j).sz_lead_vec;
    z_diff_vec = hmm_input_output(j).z_diff_vec;  
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
        spot_fragment_norm = spot_fragment-nanmin(spot_fragment);
        swap_fragment_norm = swap_fragment-nanmin(swap_fragment);
        virtual_fragment_norm = virtual_fragment-nanmin(virtual_fragment);
        % calculate protein-weighted center of mass
        pt_cm_spot_vec(iter) = nansum(time_fragment.*spot_fragment_norm) ./ nansum(spot_fragment_norm) - start_time; 
        pt_cm_swap_vec(iter) = nansum(time_fragment.*swap_fragment_norm) ./ nansum(swap_fragment_norm) - start_time; 
        pt_cm_virtual_vec(iter) = nansum(time_fragment.*virtual_fragment_norm) ./ nansum(virtual_fragment_norm) - start_time; 
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
        tr_amp_vec(iter) = sz_lag_vec(i);
        tr_amp_prev_vec(iter) = sz_lead_vec(i);
        tr_dur_vec(iter) = z_dur_lag_vec(i);
        tr_dur_prev_vec(iter) = z_dur_lead_vec(i);        
        % increment
        iter = iter + 1;  
        if mod(iter,100) == 0
            disp(iter)
        end
    end    
end        
% make results table
results_table = array2table([particle_id_vec, feature_time_vec, tr_burst_class_vec, fluo_vec, tr_amp_vec, tr_amp_prev_vec, tr_dur_vec, tr_dur_prev_vec,...
    pt_mf_vec, pt_cm_spot_vec, pt_cm_swap_vec, pt_cm_virtual_vec,...
    pt_net_spot_vec, pt_net_swap_vec, pt_net_virtual_vec, pt_net2_spot_vec, pt_net2_swap_vec, pt_net2_virtual_vec],...
    'VariableNames',{'ParticleID', 'time', 'burst_class', 'fluo', 'tr_amp', 'tr_amp_prev', 'tr_dur', 'tr_dur_prev', 'pt_mf', 'pt_cm_spot', 'pt_cm_swap', 'pt_cm_virtual',...
    'pt_net_spot', 'pt_net_swap', 'pt_net_virtual','pt_net2_spot', 'pt_net2_swap', 'pt_net2_virtual'});

% save
save([dataPath 'hmm_burst_table.mat'],'results_table')