clear
close all

% define core ID variables
project = 'Dl-Ven_snaBAC-mCh';
% project = 'Dl-Ven_hbP2P-mCh';
DropboxFolder =  'E:\Nick\LivemRNA\Dropbox (Personal)\';
[~, DataPathTarget, ~] =   header_function(DropboxFolder, project); 
DataPath = [DropboxFolder 'ProcessedEnrichmentData\' project '/'];
w = 7;
K = 3;
window_size = 15; % number of lags and leads over which to track protein/fluo dynamics
nBoots = 100; % number of bootstraps to use for calculating se
% load input-output data set
load([DataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '_dt.mat'])
% set scale for feature identification
fluo_scale = prctile([hmm_input_output.fluo],40);
% set window size
ref_vec = -window_size:window_size;
% define features to track
feature_cell = {'high','low','rise','fall'};
feature_titles = {'peaks', 'troughs', 'rises', 'falls'};

% initialize results structure
iter = 1;
feature_struct = struct;
tic
for i = 1:numel(hmm_input_output)
    %%% core data vectors
    time_vec = hmm_input_output(i).time/60;   
    frame_vec = 1:numel(time_vec);
    fluo_vec = hmm_input_output(i).fluo;
    pt_spot = hmm_input_output(i).spot_protein_dt;
    pt_serial = hmm_input_output(i).serial_protein_dt;
    pt_swap = hmm_input_output(i).swap_spot_protein_dt;
    qc_filter = hmm_input_output(i).dt_filter_gap;
    qc_filter_swap = hmm_input_output(i).dist_swap_dt_filter_gap;
    
    %%% find features
    % find fluo peaks and troughs
    [~,fluo_high_ids] = findpeaks(fluo_vec,frame_vec,'MinPeakProminence',fluo_scale); % NL: eye-balled atm
    [~,fluo_low_ids] = findpeaks(-fluo_vec,frame_vec,'MinPeakProminence',fluo_scale); % NL: eye-balled atm

    % identify changepoints     
    fluo_d_vec = sign([0 diff(fluo_vec)]);    
    ipt_fluo = findchangepts(fluo_vec,'MinThreshold',fluo_scale);

    fluo_rise_ids = ipt_fluo(fluo_d_vec(ipt_fluo)==1);
    fluo_fall_ids = ipt_fluo(fluo_d_vec(ipt_fluo)==-1);
    
    % apply filters 
    pt_spot(qc_filter) = NaN;
    pt_serial(qc_filter) = NaN;
    pt_swap(qc_filter_swap) = NaN;
    fluo_vec(qc_filter) = NaN;
    
    for j = 1:numel(feature_cell)
        feature_str = feature_cell{j};                    
        % get variables
        eval(['ids =' 'fluo_' feature_str '_ids;'])        
        % initialize arrays
        response_temp = NaN(numel(ids),2*window_size+1);            
        spot_protein_temp = NaN(numel(ids),2*window_size+1);
        serial_protein_temp = NaN(numel(ids),2*window_size+1);    
        swap_protein_temp = NaN(numel(ids),2*window_size+1);    
        time_temp = NaN(1,numel(ids));
        for m = 1:numel(ids)                
            raw_ind = ref_vec+ids(m); 
            ft_vec1 = raw_ind > 0 & raw_ind <= numel(fluo_vec);
            ft_vec2 = raw_ind(ft_vec1);
            % record
            response_temp(m,ft_vec1) = fluo_vec(ft_vec2);          
            spot_protein_temp(m,ft_vec1) = pt_spot(ft_vec2);
            serial_protein_temp(m,ft_vec1) = pt_serial(ft_vec2);  
            swap_protein_temp(m,ft_vec1) = pt_swap(ft_vec2);
            time_temp(m) = time_vec(ids(m));               
        end
        % add to main cell structures
        feature_struct(iter).(['fluo_' feature_str '_response']) = response_temp;
        feature_struct(iter).(['fluo_' feature_str '_spot_protein']) = spot_protein_temp;
        feature_struct(iter).(['fluo_' feature_str '_serial_protein']) = serial_protein_temp;
        feature_struct(iter).(['fluo_' feature_str '_swap_protein']) = swap_protein_temp;
        feature_struct(iter).(['fluo_' feature_str '_time']) = time_temp;                    
    end
    iter = iter + 1;
end
toc
tic
%%% Calculate bootstrap estimates of Mean and SE
results_struct = struct;
iter = 1;
for j = 1:numel(feature_cell)
    feature_str = feature_cell{j};        
    % get variables
    response_mat = vertcat(feature_struct.(['fluo_' feature_str '_response']));
    spot_protein_mat = vertcat(feature_struct.(['fluo_' feature_str '_spot_protein']));
    serial_protein_mat = vertcat(feature_struct.(['fluo_' feature_str '_serial_protein']));
    swap_protein_mat = vertcat(feature_struct.(['fluo_' feature_str '_swap_protein']));
    % bootstrap variables
    boot_index_vec = 1:size(response_mat,1);
    response_boot_mat = NaN(nBoots,size(response_mat,2));
    spot_protein_boot_mat = NaN(nBoots,size(response_mat,2));
    serial_protein_boot_mat = NaN(nBoots,size(response_mat,2));
    swap_protein_boot_mat = NaN(nBoots,size(response_mat,2));
    for n = 1:nBoots
        boot_ids = randsample(boot_index_vec,numel(boot_index_vec),true);
        response_boot_mat(n,:) = nanmean(response_mat(boot_ids,:));
        spot_protein_boot_mat(n,:) = nanmean(spot_protein_mat(boot_ids,:));
        serial_protein_boot_mat(n,:) = nanmean(serial_protein_mat(boot_ids,:));
        swap_protein_boot_mat(n,:) = nanmean(swap_protein_mat(boot_ids,:));
    end
    results_struct(iter).ID = ['fluo_' feature_titles{j}];
    results_struct(iter).fn = ['fluo_' feature_titles{j}];
    results_struct(iter).data_name = 'fluo';
    % response
    results_struct(iter).response_mean = nanmean(response_boot_mat);
    results_struct(iter).response_ste = nanstd(response_boot_mat);
    % spot protein
    results_struct(iter).spot_protein_mean = nanmean(spot_protein_boot_mat);
    results_struct(iter).spot_protein_ste = nanstd(spot_protein_boot_mat);
    % swap protein
    results_struct(iter).swap_protein_mean = nanmean(swap_protein_boot_mat);
    results_struct(iter).swap_protein_ste = nanstd(swap_protein_boot_mat);
    % serial protein
    results_struct(iter).serial_protein_mean = nanmean(serial_protein_boot_mat);
    results_struct(iter).serial_protein_ste = nanstd(serial_protein_boot_mat);
    % difference
    results_struct(iter).diff_protein_mean = nanmean(spot_protein_boot_mat - serial_protein_boot_mat);
    results_struct(iter).diff_protein_ste = nanstd(spot_protein_boot_mat - serial_protein_boot_mat);
    % save key data info
    results_struct(iter).fluo_scale = fluo_scale;
    results_struct(iter).ref_vec = ref_vec;
    iter = iter + 1  
end
toc
    
% save
save([DataPath 'fluo_input_output.mat'],'results_struct')