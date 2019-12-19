clear 
close all

% define core ID variables
project = 'Dl-Ven_snaBAC-mCh';
% project = 'Dl-Ven_hbP2P-mCh';
dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
dataPath = [dropboxFolder 'ProcessedEnrichmentData\' project '/'];
% load input-output data set
K = 3;
w = 7;
load([dataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'],'hmm_input_output')

% Attempt to track protein "burst" characteristics over time and [MF]
% at locus
spot_pt_burst_sep_vec = NaN(numel([hmm_input_output.spot_protein]),1);
spot_pt_burst_dur_vec = NaN(numel([hmm_input_output.spot_protein]),1);
spot_pt_burst_amp_vec = NaN(numel([hmm_input_output.spot_protein]),1);
spot_pt_mf_vec = NaN(numel([hmm_input_output.spot_protein]),1);
spot_pt_time_vec = NaN(numel([hmm_input_output.spot_protein]),1);
% at control spot
virt_pt_burst_sep_vec = NaN(numel([hmm_input_output.spot_protein]),1);
virt_pt_burst_dur_vec = NaN(numel([hmm_input_output.spot_protein]),1);
virt_pt_burst_amp_vec = NaN(numel([hmm_input_output.spot_protein]),1);
virt_pt_mf_vec = NaN(numel([hmm_input_output.spot_protein]),1);
virt_pt_time_vec = NaN(numel([hmm_input_output.spot_protein]),1);
% iterate through spot preotein first
iter = 1;
for i = 1:numel(hmm_input_output)
    % extract vectors
    spot_pt_vec = imgaussfilt(hmm_input_output(i).spot_protein_dt,1);
%     spot_pt_vec = imgaussfilt(spot_pt_vec,g1)-imgaussfilt(spot_pt_vec,g2);    
    time_vec = hmm_input_output(i).time;
    mf_vec = hmm_input_output(i).mf_protein;
    gap_vec = hmm_input_output(i).dt_filter_gap;
    spot_pt_vec(gap_vec) = NaN;
    % make binary vec
    spot_pt_bin = spot_pt_vec >= 0;
    dp_vec = [0 diff(spot_pt_bin)];
    ids = find(dp_vec~=1);
    rise_ids = find(dp_vec==1);
    fall_ids = find(dp_vec==-1);
    for j = rise_ids
        prev_fall = fall_ids(find(fall_ids < j,1,'last'));
        prev_nan = sum(isnan(spot_pt_vec(prev_fall+1:j)));
        if ~isempty(prev_fall) && prev_nan == 0
            spot_pt_burst_sep_vec(iter) = j-prev_fall;
        end
        next_fall = fall_ids(find(fall_ids > j,1));
        next_nan = sum(isnan(spot_pt_vec(j:next_fall-1)));
        if ~isempty(next_fall) && ~next_nan
            spot_pt_burst_dur_vec(iter) = next_fall-j;
            spot_pt_burst_amp_vec(iter) = sum(spot_pt_vec(j:next_fall-1)) / spot_pt_burst_dur_vec(iter);
        end
        spot_pt_mf_vec(iter) = mf_vec(j);
        spot_pt_time_vec(iter) = time_vec(j);
        % increment
        iter = iter + 1;
    end
end
last_id = find(~isnan(spot_pt_time_vec),1,'last');
spot_pt_mf_vec = spot_pt_mf_vec(1:last_id);
spot_pt_time_vec = spot_pt_time_vec(1:last_id);
spot_pt_burst_sep_vec = spot_pt_burst_sep_vec(1:last_id);
spot_pt_burst_amp_vec = spot_pt_burst_amp_vec(1:last_id);
spot_pt_burst_dur_vec = spot_pt_burst_dur_vec(1:last_id);
spot_pt_size_vec = spot_pt_burst_dur_vec.*spot_pt_burst_amp_vec;

% now control spot
iter = 1;
for i = 1:numel(hmm_input_output)
    % extract vectors
    spot_pt_vec = imgaussfilt(hmm_input_output(i).serial_protein_dt,1);
%     spot_pt_vec = imgaussfilt(spot_pt_vec,g1)-imgaussfilt(spot_pt_vec,g2);    
    time_vec = hmm_input_output(i).time;
    mf_vec = hmm_input_output(i).mf_protein;
    gap_vec = hmm_input_output(i).dt_filter_gap;
    spot_pt_vec(gap_vec) = NaN;
    % make binary vec
    spot_pt_bin = spot_pt_vec >= 0;
    dp_vec = [0 diff(spot_pt_bin)];
    ids = find(dp_vec~=1);
    rise_ids = find(dp_vec==1);
    fall_ids = find(dp_vec==-1);
    for j = rise_ids
        prev_fall = fall_ids(find(fall_ids < j,1,'last'));
        prev_nan = sum(isnan(spot_pt_vec(prev_fall+1:j)));
        if ~isempty(prev_fall) && prev_nan == 0
            virt_pt_burst_sep_vec(iter) = j-prev_fall;
        end
        next_fall = fall_ids(find(fall_ids > j,1));
        next_nan = sum(isnan(spot_pt_vec(j:next_fall-1)));
        if ~isempty(next_fall) && ~next_nan
            virt_pt_burst_dur_vec(iter) = next_fall-j;
            virt_pt_burst_amp_vec(iter) = sum(spot_pt_vec(j:next_fall-1)) / virt_pt_burst_dur_vec(iter);
        end
        virt_pt_mf_vec(iter) = mf_vec(j);
        virt_pt_time_vec(iter) = time_vec(j);
        % increment
        iter = iter + 1;
    end
end
last_id = find(~isnan(spot_pt_time_vec),1,'last');
virt_pt_mf_vec = virt_pt_mf_vec(1:last_id);
virt_pt_time_vec = virt_pt_time_vec(1:last_id);
virt_pt_burst_sep_vec = virt_pt_burst_sep_vec(1:last_id);
virt_pt_burst_amp_vec = virt_pt_burst_amp_vec(1:last_id);
virt_pt_burst_dur_vec = virt_pt_burst_dur_vec(1:last_id);
virt_pt_size_vec = virt_pt_burst_dur_vec.*virt_pt_burst_amp_vec;

mf_vec_full = [hmm_input_output.mf_protein];
time_vec_full = [hmm_input_output.time];

%%
close all
size_bins = 0:40:440;
mf_bins = 0:30:360;
total_hist_spot = histc(mf_vec_full,mf_bins);

size_hist_spot = hist3([spot_pt_size_vec spot_pt_mf_vec],'Edges',{size_bins,mf_bins});
size_hist_spot = size_hist_spot./total_hist_spot;

size_hist_virt = hist3([virt_pt_size_vec virt_pt_mf_vec],'Edges',{size_bins,mf_bins});
size_hist_virt = size_hist_virt./total_hist_spot;



virt_fig = figure;
hm_cm = flipud(brewermap([],'RdYlBu'));
colormap(hm_cm);
imagesc(size_hist_virt)
colorbar

spot_fig = figure;
colormap(hm_cm);
imagesc(size_hist_spot)
colorbar