clear
close all

% define core ID variables
project = 'Dl-Ven_New_snaBAC-mCh';
% project = 'Dl-Ven_hbP2P-mCh';
DropboxFolder =  'E:\Nick\LivemRNA\Dropbox (Personal)\';
dataPath = [DropboxFolder 'ProcessedEnrichmentData\' project '/'];
figPath = [DropboxFolder 'LocalEnrichmentFigures\' project '/'];
mkdir(figPath);
w = 7;
K = 3;
% establish time classes to investigate
time_class_cell = {1:60};%1:15,16:25,26:35,36:45,46:60};
fluo_feature_quantile_vec = [11];%[33,11,11,11,11,11];
% window analysis params
window_size = 15; 
out_quantiles = 11; % number of quantiles to use for output arrays
% load input-output data set
load([dataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'],'hmm_input_output')
% generate average fluo change vector
kernel_size = 1;
for i = 1:numel(hmm_input_output)    
    fluo_interp = hmm_input_output(i).fluo;
    fluo_diff = [0 diff(fluo_interp)];
    hmm_input_output(i).fluo_diff_sm = imgaussfilt(fluo_diff,kernel_size);    
end
% get fluo change deciles
fluo_sm_vec = [hmm_input_output.fluo_diff_sm];
gap_filter_vec = [hmm_input_output.dt_filter_gap];
time_vec = round([hmm_input_output.time]/60);
% initialize results structure
results_struct = struct;
for ti = 1:numel(time_class_cell)
    fluo_feature_quantiles = fluo_feature_quantile_vec(ti); % number of quantiles to use for feature classification
    time_range = time_class_cell{ti};
    center_time = round(mean(time_range));
    time_filter = ismember(time_vec,time_range); 
    iter_filter = ~gap_filter_vec&time_filter;
    fluo_diff_quantile_vec = [nanmin(fluo_sm_vec(iter_filter)) quantile(fluo_sm_vec(iter_filter),fluo_feature_quantiles-2) nanmax(fluo_sm_vec(iter_filter))];
    % calculate median and quartiles of protein trends for each decile range
    % fluorescence
    fluo_activity_array = NaN(numel(fluo_diff_quantile_vec)-1,2*window_size+1);
    % target locus
    spot_protein_array = NaN(numel(fluo_diff_quantile_vec)-1,2*window_size+1);
    % swap locus
    swap_protein_array = NaN(numel(fluo_diff_quantile_vec)-1,2*window_size+1);
    % virtual spot
    virtual_protein_array = NaN(numel(fluo_diff_quantile_vec)-1,2*window_size+1);
    window_vec = -window_size:window_size;
    fit_array = [ones(numel(window_vec),1) window_vec'];
    for i = 1:numel(fluo_diff_quantile_vec)-1
        lb = fluo_diff_quantile_vec(i);
        ub = fluo_diff_quantile_vec(i+1);
        n_entries = sum(fluo_sm_vec(~gap_filter_vec)>=lb&fluo_sm_vec(~gap_filter_vec)<ub);
        % initialize arrays
        fluo_array = NaN(n_entries,2*window_size+1);
        spot_array = NaN(n_entries,2*window_size+1);
        swap_array = NaN(n_entries,2*window_size+1);
        virtual_array = NaN(n_entries,2*window_size+1);
        % iterate through structure
        iter = 1;
        for j = 1:numel(hmm_input_output)
            gap_filter = hmm_input_output(j).dt_filter_gap;
            fluo_diff_sm = hmm_input_output(j).fluo_diff_sm;
            fluo = hmm_input_output(j).fluo;
            % protein fields                        
            spot_protein = hmm_input_output(j).spot_protein;
            swap_spot_protein = hmm_input_output(j).dist_swap_spot_protein;
            virtual_protein = hmm_input_output(j).serial_protein;        
            % apply filter
            fluo_diff_sm(gap_filter) = NaN;
            spot_protein(gap_filter) = NaN;
            swap_spot_protein(gap_filter) = NaN;
            virtual_protein(gap_filter) = NaN;
            % find features
            id_list = find(fluo_diff_sm>=lb&fluo_diff_sm<ub);
            for id = id_list
                full_range = id - window_size:id+window_size;
                true_range = full_range(full_range>0&full_range<=numel(virtual_protein));
                % record
                ft1 = ismember(full_range,true_range);
                if sum(~isnan(spot_protein)) >= window_size && sum(~isnan(swap_spot_protein)) >= window_size && sum(~isnan(virtual_protein)) >= window_size
                    spot_fragment = spot_protein(true_range);
                    swap_fragment = swap_spot_protein(true_range);
                    virtual_fragment = virtual_protein(true_range);
                    fluo_fragment = fluo(true_range);
                    % fit linear offsets            
                    fit_sub_array = fit_array(ft1,:);
                    spot_fit = fit_sub_array(~isnan(spot_fragment),:) \ spot_fragment(~isnan(spot_fragment))';
                    swap_fit = fit_sub_array(~isnan(swap_fragment),:) \ swap_fragment(~isnan(swap_fragment))';
                    virtual_fit = fit_sub_array(~isnan(virtual_fragment),:) \ virtual_fragment(~isnan(virtual_fragment))';                
                    fluo_fit = fit_sub_array(~isnan(fluo_fragment),:) \ fluo_fragment(~isnan(fluo_fragment))';
                    % save                
                    spot_array(iter,ft1) = spot_fragment - spot_fit(1) - spot_fit(2)*window_vec(ft1);
                    swap_array(iter,ft1) = swap_fragment - swap_fit(1) - swap_fit(2)*window_vec(ft1);
                    virtual_array(iter,ft1) = virtual_fragment - virtual_fit(1) - virtual_fit(2)*window_vec(ft1);
                    fluo_array(iter,ft1) = fluo_fragment - fluo_fit(1) - fluo_fit(2)*window_vec(ft1);
                end
                % increment
                iter = iter + 1;            
            end    
        end
        % record
        fluo_activity_array(i,:) = nanmean(fluo_array)';    
        spot_protein_array(i,:) = nanmean(spot_array)';                   
        swap_protein_array(i,:) = nanmean(swap_array)';               
        virtual_protein_array(i,:) = nanmean(virtual_array)';    
    end
    results_struct(ti).fluo_activity_array = fluo_activity_array;
    results_struct(ti).spot_protein_array = spot_protein_array;
    results_struct(ti).swap_protein_array = swap_protein_array;
    results_struct(ti).virtual_protein_array = virtual_protein_array;
    results_struct(ti).time_range = time_range;
    results_struct(ti).window_vec = window_vec;
    results_struct(ti).fluo_diff_deciles = fluo_diff_quantile_vec;
end            

%%% Make figures
close all
% % make rise series figs
% yw = [234 194 100]/256; % yellow
% bl = [115 143 193]/256; % blue
% rd = [213 108 85]/256; % red
% gr = [191 213 151]/256; % green
% br = [207 178 147]/256; % brown
% save
% save([dataPath 'fluo_input_output_results.mat'],'results_struct')