% Script to probe relationship between input protein concentration and
% output transcriptional response
clear 
close all
% define ID variables
K = 3;
w = 7;
project = 'Dl_Venus_snaBAC_mCherry_Leica_hp';
dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
dataPath = [dropboxFolder '\ProcessedEnrichmentData\' project '\'];
figPath = [dropboxFolder '\LocalEnrichmentFigures\' project '\hmm_input_output_K' num2str(K) '_w' num2str(w) '\'];

nTraces = 50; % number of individual traces to select for plotting
window_size = 10; % number of lags over which to track protein/fluo dynamics
n_boots = 100;
n_ref_hist_bins = 10;
min_time = 10;
% extract protein, gene, fluorophore info
underscores = strfind(project,'_');
protein_name = project(1:underscores(1)-1);
protein_fluor = project(underscores(1)+1:underscores(2)-1);
gene_name = project(underscores(2)+1:underscores(3)-1);
if numel(underscores) == 3
    ind = numel(project);
else
    ind = underscores(4)-1;
end
gene_fluor = project(underscores(3)+1:end);


% load data set
load([dataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'])
%%
% first make figures to ensure that hmmm results have been properly
% concatenated with protein data
for i = 1:numel(master_struct)
    subProject = master_struct(i).project;
    subID = master_struct(i).ID;
    qcPath = [figPath '\' subID '_qc_' subProject '\'];
    mkdir(qcPath);
    hmm_input_output = master_struct(i).hmm_input_output;
    s_index = 1:numel(hmm_input_output);
    rng(123);
    plot_indices = randsample(s_index,min([20,numel(s_index)]),false);
    for j = 1:numel(plot_indices)
        % MCP channel checks
        mcp_check = hmm_input_output(plot_indices(j)).mcp_check;
        fluo_check = hmm_input_output(plot_indices(j)).fluo_check;
        fluo = hmm_input_output(plot_indices(j)).fluo;
        time = hmm_input_output(plot_indices(j)).time;
        r_vec = sum(hmm_input_output(plot_indices(j)).r_mat,2);
        % make figure
        qc_fig = figure('Visible','off');
        hold on
        plot(time,fluo / nanmean(fluo))
        plot(time,fluo_check / nanmean(fluo_check));
        plot(time,mcp_check / nanmean(mcp_check))
        plot(time,r_vec / nanmean(r_vec))
        legend('fluo (HMM)', 'fluo (data)','raw mcp','activity state (HMM)')
        xlabel('time')
        ylabel([gene_name ' activity (au)'])
        saveas(qc_fig,[qcPath 'mcp_check_nc_' sprintf('%03d',plot_indices(j)) '.png'])
        
        % Protein Channel checks
        spot_protein = hmm_input_output(plot_indices(j)).spot_protein;
        null_protein = hmm_input_output(plot_indices(j)).null_protein;
        mf_protein = hmm_input_output(plot_indices(j)).mf_protein;
        % make figure
        qc_fig = figure('Visible','off');
        hold on
        plot(time,spot_protein)
        plot(time,null_protein)
        plot(time,mf_protein)
        legend('protein (spot)', 'protein (control spot)','protein (mf control)')
        xlabel('time')
        ylabel([protein_name ' - ' protein_fluor ' (au)'])
        saveas(qc_fig,[qcPath 'protein_check_nc_' sprintf('%03d',plot_indices(j)) '.png'])
    end
end

% Make single trace input-output plots
for i = 1:numel(master_struct)
    subProject = master_struct(i).project;
    subID = master_struct(i).ID;
    tracePath = [figPath '\' subID '_single_trace_' subProject '\'];
    mkdir(tracePath);
    hmm_input_output = master_struct(i).hmm_input_output;
    s_index = 1:numel(hmm_input_output);
    rng(321);
    plot_indices = randsample(s_index,min([nTraces,numel(s_index)]),false);
    for j = 1:numel(plot_indices)
        % MCP channel checks
        time = hmm_input_output(plot_indices(j)).time;
        r_vec = sum(hmm_input_output(plot_indices(j)).r_mat,2);
        spot_protein = hmm_input_output(plot_indices(j)).spot_protein;
        null_protein = hmm_input_output(plot_indices(j)).null_protein;
        delta_protein = spot_protein - null_protein;
        % make figure
        trace_fig = figure('Visible','off');
        hold on
        
        yyaxis left
        plot(time,r_vec)
        ylabel(['instantaneous ' gene_name ' activity (au)'])
        
        yyaxis right
        plot(time,delta_protein);
        ylabel(['absolute ' protein_name ' enrichment(au)'])
                
        ax = gca;
        ax.YAxis(1).Color = 'black';
        ax.YAxis(2).Color = 'black';
        
        legend('transcriptional activity', 'local protein concentration')
        xlabel('time')
        ylabel([gene_name ' activity (au)'])
        saveas(trace_fig,[tracePath 'input_output_nc' sprintf('%03d',plot_indices(j)) '.png'])               
    end
end
%%
close all
% ignore possibility of second project for noq 
hmm_input_output = master_struct(1).hmm_input_output;
Tres = nanmedian(diff([hmm_input_output.time]));
% Calcualte average cross-correlation. Use score for randomly assigned
% % trace pairs as a control
% xc_index_vec = 1:numel(hmm_input_output);
% xc_shuffle = randsample(xc_index_vec,numel(xc_index_vec),false);
% 
% xc_trend_mat = NaN(numel(xc_shuffle),2*window_size + 1);
% xc_control_mat = NaN(numel(xc_shuffle),2*window_size + 1);
% xc_trend_ct_mat = NaN(numel(xc_shuffle),2*window_size + 1);
% xc_control_ct_mat = NaN(numel(xc_shuffle),2*window_size + 1);
% 
% for i = 1:numel(xc_shuffle)
%     % compare true pairs first
%     delta_pt = hmm_input_output(i).spot_protein;%-hmm_input_output(i).mf_protein;
%     r_vec = sum(hmm_input_output(i).r_mat,2)';   
%     if numel(delta_pt) > window_size + 1
%         xc_trend_mat(i,:) = xcov(delta_pt,r_vec,window_size);
%         xc_trend_ct_mat(i,:) = [numel(delta_pt)-window_size:numel(delta_pt) fliplr(numel(delta_pt)-window_size:numel(delta_pt)-1)];
%     end
%     
%     % draw random activity vector and take cross-covariance
%     r_vec_rand = sum(hmm_input_output(xc_shuffle(i)).r_mat,2)';
%     n_max = min([numel(r_vec_rand) numel(delta_pt)]);
%     if n_max >= window_size + 1
%         xc_control_mat(i,:) = xcov(delta_pt(1:n_max),(r_vec_rand(1:n_max)),window_size);
%         xc_control_ct_mat(i,:) = [n_max-window_size:n_max fliplr(n_max-window_size:n_max-1)];
%     end
% end
% 
% trend_boots = NaN(n_boots,size(xc_control_mat,2));
% control_boots = NaN(n_boots,size(xc_control_mat,2));
% 
% trend_index_vec = find(~isnan(xc_trend_mat(:,1)));
% control_index_vec = find(~isnan(xc_trend_mat(:,1)));
% for i = 1:n_boots
%     trend_ids = randsample(trend_index_vec,numel(trend_index_vec),true);
%     trend_boots(i,:) = nansum(xc_trend_mat(trend_ids,:).*xc_trend_ct_mat(trend_ids,:)) ./ nansum(xc_trend_ct_mat(trend_ids,:));
%     control_ids = randsample(control_index_vec,numel(control_index_vec),true);
%     control_boots(i,:) = nansum(xc_control_mat(control_ids,:).*xc_control_ct_mat(control_ids,:)) ./ nansum(xc_control_ct_mat(control_ids,:));
% end
% xc_trend_mean = nanmean(trend_boots);
% xc_trend_ste = nanstd(trend_boots);
% 
% xc_control_mean = nanmean(control_boots);
% xc_control_ste = nanstd(control_boots);
% 
% % sketch method for extracting time resolution 
% 
% lag_time_vec = (-window_size:0)*Tres;
% 
% % make figure
% xc_fig = figure;
% hold on
% e1 = errorbar([lag_time_vec -fliplr(lag_time_vec(1:end-1))],xc_trend_mean,xc_trend_ste);
% e1.CapSize = 0;
% e2 = errorbar([lag_time_vec -fliplr(lag_time_vec(1:end-1))],xc_control_mean,xc_control_ste);
% e2.CapSize = 0;
% legend('active locus','random control')
% xlabel('offset (seconds)')
% ylabel('cross-covariance')
% title(['cross-covariance between ' gene_name ' activity and ' protein_name ' enrichment'])
% saveas(xc_fig,[figPath 'cross-covariance.png'])

% Generate effective transcription rate vector
occ_vec = nanmean(vertcat(hmm_input_output.z_mat));
r_mean = nanmean(vertcat(hmm_input_output.r));
if K==3
    r_mean = [r_mean(1) (r_mean(2).*occ_vec(2)+r_mean(3).*occ_vec(3))/(occ_vec(2)+occ_vec(3))];
end
for i = 1:numel(hmm_input_output)
    z_mat = hmm_input_output(i).z_mat;
    r = hmm_input_output(i).r;    
    if K == 3
        z_mat(:,2) = sum(z_mat(:,2:3),2);
        z_mat = z_mat(:,1:2);
    end
    r_mat = z_mat .* r_mean;
    hmm_input_output(i).r_vec = sum(r_mat,2);
    hmm_input_output(i).z_mat2 = z_mat;
    [~,z_vec] = max(z_mat,[],2);
    hmm_input_output(i).z_vec = z_vec;
end


% calculate prduction rate extrema

r_unit = r_mean(2);
% create indexing and storage vectors
timeBins = linspace(min_time,nanmax([hmm_input_output.time]/60),n_ref_hist_bins+1);
dT = median(diff(timeBins));
ptBins = linspace(0,prctile([hmm_input_output.null_protein],95),n_ref_hist_bins+1);
dP = median(diff(ptBins));

trend_id_cell = {'burst start','burst stop','protein peak','protein dip','consistency check'};

metric_cell = {'delta_pt','delta_pt','r_vec','r_vec','delta_pt'}; % 

measure_cell = {'burst_rise_times','burst_fall_times','blip_times','dip_times','rand_ids'};
% condition_cell = {'delta_r>.5*r_unit','delta_r<-.5*r_unit','bulk_delta_r>.2*r_unit',...
%     'bulk_delta_r<-.2*r_unit','r_sum>mean_steps*.8*r_unit','r_sum<mean_steps*.2*r_unit',...
%     'delta_smooth>delta_high','delta_smooth<delta_low',...
%     'ismember(1:numel(delta_r),randsample(1:numel(delta_r),2,true))'};

trend_index_cell = cell(1,numel(trend_id_cell));
trend_subindex_cell = cell(1,numel(trend_id_cell));
trend_linear_index_cell = cell(1,numel(trend_id_cell));

control_index_cell = cell(1,numel(trend_id_cell));
control_subindex_cell = cell(1,numel(trend_id_cell));

% calculate high and low thresholds for enrichment 
pt_prom = prctile([hmm_input_output.spot_protein],80)-prctile([hmm_input_output.spot_protein],20);

trend_weight_array = zeros(numel(timeBins)-1,numel(timeBins)-1,numel(trend_id_cell));
% generate longform ref indices as we go
ind_vec_full = [];
sub_ind_vec_full = [];
iter = 1;
% collect rise and fall information
for i = 1:numel(hmm_input_output)
    time_vec = hmm_input_output(i).time/60;
    frame_vec = 1:numel(time_vec);
    % Analyze MCP Channel first
    % find all peaks and troughs
    z_vec = hmm_input_output(i).z_vec;    
    [~,burst_starts,burst_widths] = findpeaks(z_vec,frame_vec);
    [~,burst_stops,trough_widths] = findpeaks(-z_vec,frame_vec);
    burst_peaks = burst_starts + burst_widths / 2;
    burst_troughs = burst_stops + trough_widths / 2;
    % find BIG peaks and troughs
    big_burst_peaks = burst_peaks(burst_widths > round(1/Tres*60));
    big_burst_troughs = burst_troughs(trough_widths > round(.7/Tres*60));
    % identify transitions between (uninterrupted) big troughs and big bursts
    burst_rise_times = [];
    for j = 1:numel(big_burst_troughs)
        next_start = burst_starts(find(burst_peaks>big_burst_troughs(j),1));
        next_peak = burst_peaks(find(burst_peaks>big_burst_troughs(j),1));
        next_big_peak = big_burst_peaks(find(big_burst_peaks>big_burst_troughs(j),1));
        if next_big_peak == next_peak
            burst_rise_times = [burst_rise_times next_start];
        end
    end
    burst_fall_times = [];
    for j = 1:numel(big_burst_peaks)
        next_stop = burst_stops(find(burst_troughs>big_burst_peaks(j),1));
        next_trough = burst_troughs(find(burst_troughs>big_burst_peaks(j),1));
        next_big_trough = big_burst_troughs(find(big_burst_troughs>big_burst_peaks(j),1));
        if next_big_trough == next_trough
            burst_fall_times = [burst_fall_times next_stop];
        end
    end
    % Find large protein blips and dips    
    null_vec = hmm_input_output(i).null_protein;
    spot_vec = hmm_input_output(i).spot_protein;
    spot_vec_filtered = imgaussfilt(spot_vec,2);

    [~,blip_times,blip_widths] = findpeaks(spot_vec_filtered,frame_vec,'MinPeakWidth',60*.5/Tres,'MinPeakProminence',.4*pt_prom); % NL: eye-balled atm
    [~,dip_times,dip_widths] = findpeaks(-spot_vec_filtered,frame_vec,'MinPeakWidth',60*.5/Tres,'MinPeakProminence',.4*pt_prom); % NL: eye-balled atm
    
    ind_vec_full = [ind_vec_full repelem(i,numel(time_vec))];
    sub_ind_vec_full = [sub_ind_vec_full frame_vec];
    linear_index_full = iter:iter+numel(time_vec)-1;
    early_times = find(time_vec<min_time);
    iter = iter + numel(time_vec);
    
    % random control ids to ensure method consistency
    rand_ids = randsample(1:numel(time_vec),2,true);
    
    % store ids corresponding to each condition
    for j = 1:numel(trend_id_cell)  
        
        id_vec = trend_index_cell{j};
        sub_id_vec = trend_subindex_cell{j};        
        lin_vec = trend_linear_index_cell{j};
        
        % record
        trend_index_cell{j} = [id_vec repelem(i,numel(eval(measure_cell{j})))];
        trend_subindex_cell{j} = [sub_id_vec eval(measure_cell{j})];
        trend_linear_index_cell{j} = [lin_vec linear_index_full(eval(measure_cell{j}))];
        % record weights
        N = histcounts2(time_vec(eval(measure_cell{j})),null_vec(eval(measure_cell{j})),timeBins,ptBins);
        trend_weight_array(:,:,j) = trend_weight_array(:,:,j) + N;        
    end
end

% now calculate sampling weights for control
trend_weight_array = trend_weight_array ./ sum(sum(trend_weight_array,1),2);
control_weights = NaN(size(trend_weight_array));
% make ref vectors
time_ref = [hmm_input_output.time]/60;
time_ref(time_ref<min_time) = NaN;
null_ref = [hmm_input_output.null_protein]+1e-6;

for i = 1:numel(trend_id_cell)
    % require that control samples fall outside of set of interest
    ind_vec = trend_index_cell{i};
    subind_vec = trend_subindex_cell{i};
    lin_vec = trend_linear_index_cell{i};
   
    stopDist = bwdist([0 diff(ind_vec)]==1);
    startDist = bwdist([diff(ind_vec) 0]==1);
    ft = lin_vec|startDist<3|startDist<3;
    
    t_temp = time_ref;
    t_temp(ft) = NaN;
    null_temp = null_ref;
    null_temp(ft) = NaN;
    
    base_weights = histcounts2(t_temp,null_temp,timeBins,ptBins);%,'Normalization','probability')';
    base_weights = base_weights / sum(base_weights(:));
    control_weights(:,:,i) = trend_weight_array(:,:,i) ./ base_weights;
end
control_weights(isnan(control_weights)|isinf(control_weights)) = 0;

% map weights back to observations
samp_wt_array = NaN(numel(sub_ind_vec_full),numel(trend_id_cell));
for i = 1:numel(trend_id_cell)
    ind_vec = trend_index_cell{i};
    subind_vec = trend_subindex_cell{i};
    lin_vec = trend_linear_index_cell{i};
    
    stopDist = bwdist([0 diff(ind_vec)]==1);
    startDist = bwdist([diff(ind_vec) 0]==1);
    ft = lin_vec|startDist<3|startDist<3;
    
    wt_slice = control_weights(:,:,i);
    
    x = ceil(null_ref/dP);
    x(x>n_ref_hist_bins) = n_ref_hist_bins;
    y = ceil(time_ref/dT);
    y(y>n_ref_hist_bins) = n_ref_hist_bins;
    
    wt_indices = sub2ind(size(wt_slice),y,x);
    nan_filter = isnan(wt_indices);
    % make sure we don't draw trend data points for control
    samp_wt_array(~nan_filter,i) = wt_slice(wt_indices(~nan_filter));
    samp_wt_array(ft,i) = 0;
    samp_wt_array(isnan(null_ref)|isnan(time_ref),i) = 0;
end

% draw sample indices
sample_weights_array = NaN(size(control_weights));
wt_qc_path = [figPath '/control_weight_checks/'];
mkdir(wt_qc_path);

for i = 1:numel(trend_id_cell)    
    rng(432); % seed for consistency
    ids = randsample(1:numel(null_ref),numel(trend_index_cell{i}),true,samp_wt_array(:,i));
    control_index_cell{i} = ind_vec_full(ids);
    control_subindex_cell{i} = sub_ind_vec_full(ids);
    % make figures to check that distribution of controls match trend
    check_weights = histcounts2(time_ref(ids),null_ref(ids),timeBins,ptBins);
    check_weights = check_weights / sum(check_weights(:));
    sample_weights_array(:,:,i) = check_weights;
    trend_weights = trend_weight_array(:,:,i);
    lb = prctile(trend_weights(:),5);
    ub = prctile(trend_weights(:),95);
    
    check_fig = figure('Visible','off','Position',[0 0 1024 512]);
    subplot(1,2,1)
    imagesc(trend_weights)
    colorbar
    caxis([lb ub])
    title(['Actual sample: ' trend_id_cell{i}])
    xlabel('mf protein concentration')
    ylabel('time')
    
    subplot(1,2,2)
    imagesc(check_weights)
    colorbar
    caxis([lb ub])
    title(['Control sample: ' trend_id_cell{i}])
    xlabel('mf protein concentration')
    ylabel('time')
    saveas(check_fig,[wt_qc_path trend_id_cell{i} '.png'])
end



% initialize data arrays tpo store lag pt info
lag_trend_array_cell = cell(1,numel(trend_id_cell));
lag_control_array_cell = cell(1,numel(trend_id_cell));
for i = 1:numel(trend_id_cell)
    lag_trend_array_cell{i} = NaN(numel(trend_index_cell{i}),2*window_size+1);
    lag_control_array_cell{i} = NaN(numel(control_index_cell{i}),2*window_size+1);
end
% iterate through ids for each event type and store adjacent enrichment
% signatures
lag_index = -window_size:window_size;
ctrl_cell = {'trend','control'};
% pull trend samples
for n = 1:numel(ctrl_cell)
    ctrl_id =ctrl_cell{n};
    for i = 1:numel(trend_id_cell)
        eval(['lag_array = lag_' ctrl_id '_array_cell{i};']);
        eval(['id_vec = ' ctrl_id '_index_cell{i};']);
        eval(['sub_id_vec = ' ctrl_id '_subindex_cell{i};']);
        for j = 1:numel(id_vec)
            % protein
            spot_pt = hmm_input_output(id_vec(j)).spot_protein;
            null_pt = hmm_input_output(id_vec(j)).null_protein;
            delta_pt = spot_pt-null_pt;
            % transcription
            r_vec = hmm_input_output(id_vec(j)).r_vec;
             % index
            index_vec = lag_index + sub_id_vec(j);
            ft = index_vec>=1 & index_vec<=numel(null_pt);
            % lag     
            eval(['lag_array(j,ft) = ' metric_cell{i} '(index_vec(ft));']);
        end
         eval(['lag_' ctrl_id '_array_cell{i} =  lag_array;']);
    end
end

trend_mean_array = NaN(2*window_size+1,numel(trend_id_cell));
trend_ste_array = NaN(2*window_size+1,numel(trend_id_cell));
control_mean_array = NaN(2*window_size+1,numel(trend_id_cell));
control_ste_array = NaN(2*window_size+1,numel(trend_id_cell));

% calculate bootstrap average and ste    
for i = 1:numel(trend_id_cell)
    trend_lag_array = lag_trend_array_cell{i};
    control_lag_array = lag_control_array_cell{i};
    % convenience v ec and bootstrap arrays
    index_vec = 1:size(trend_lag_array,1);
    trend_boot = NaN(n_boots,2*window_size+1);
    control_boot = NaN(n_boots,2*window_size+1);
    for j = 1:n_boots
        s_ids = randsample(index_vec,numel(index_vec),true);
        trend_boot(j,:) = nanmean(trend_lag_array(s_ids,:));
        control_boot(j,:) = nanmean(control_lag_array(s_ids,:));
    end
    trend_mean_array(:,i) = nanmean(trend_boot);
    trend_ste_array(:,i) = nanstd(trend_boot);
    
    control_mean_array(:,i) = nanmean(control_boot);
    control_ste_array(:,i) = nanstd(control_boot);
    
    % generate figure
    lag_fig = figure;
    hold on
    e1 = errorbar(lag_index*Tres,trend_mean_array(:,i),trend_ste_array(:,i));
    e1.CapSize = 0;
    e2 = errorbar(lag_index*Tres,control_mean_array(:,i),control_ste_array(:,i));
    e2.CapSize = 0;
    legend([trend_id_cell{i} ' events'],'random control')
    xlabel('offset (seconds)')
    ylabel(['absoulte ' protein_name ' enrichment (au)'])
    title([trend_id_cell{i} ' events'])
    saveas(lag_fig,[figPath trend_id_cell{i} '.png'])
end
    
