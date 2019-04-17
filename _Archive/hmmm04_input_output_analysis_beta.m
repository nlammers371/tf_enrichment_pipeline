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
mkdir(figPath)

nTraces = 50; % number of individual traces to select for plotting
window_size = 10; % number of lags over which to track protein/fluo dynamics
nBoots = 100;
n_ref_hist_bins = 10;
min_time = 10;
make_trace_plots = 0;
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

% first make figures to ensure that hmmm results have been properly
% concatenated with protein data
if make_trace_plots
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
end

close all
% ignore possibility of second project for noq 
hmm_input_output = master_struct(1).hmm_input_output;
Tres = nanmedian(diff([hmm_input_output.time]));

% Generate effective transcription rate vector and de-trend protein data
occ_vec = nanmean(vertcat(hmm_input_output.z_mat));
r_mean = nanmean(vertcat(hmm_input_output.r));
if K==3
    r_mean = [r_mean(1) (r_mean(2).*occ_vec(2)+r_mean(3).*occ_vec(3))/(occ_vec(2)+occ_vec(3))];
end
for i = 1:numel(hmm_input_output)
    % generate discrete state and continuous production vectors
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
% obtain time-averaged fluo trend
t_vec = round([hmm_input_output.time]/60);
f_vec = [hmm_input_output.fluo];
r_vec = vertcat(hmm_input_output.r_vec);

time_index = 1:50;
fluo_trend = NaN(size(time_index));
r_trend = NaN(size(time_index));
for t = 1:numel(time_index)
    fluo_trend(t) = nanmean(f_vec(t_vec==time_index(t)));
    r_trend(t) = nanmean(r_vec(t_vec==time_index(t)));
end
filter_size = 5;
% detrend input and output time series
for i = 1:numel(hmm_input_output)
    % detrend spot protein
    spot_vec = hmm_input_output(i).spot_protein;
    frame_vec = 1:numel(spot_vec);
    [p,~,mu] = polyfit(frame_vec(~isnan(spot_vec)),spot_vec(~isnan(spot_vec)),2);
    spot_trend = polyval(p,(1:numel(spot_vec)),[],mu);
    hmm_input_output(i).spot_protein_dt = spot_vec - spot_trend;
    hmm_input_output(i).spot_protein_dt_smooth = imgaussfilt(hmm_input_output(i).spot_protein_dt,filter_size);
    % detrend control 
    null_vec = hmm_input_output(i).null_protein;
    [p,~,mu] = polyfit(frame_vec(~isnan(null_vec)),null_vec(~isnan(null_vec)),2);
    null_trend = polyval(p,(1:numel(null_vec)),[],mu);
    hmm_input_output(i).null_protein_dt = null_vec - null_trend;
    % add smoothed version of control 
    hmm_input_output(i).null_protein_sm = imgaussfilt(hmm_input_output(i).null_protein_all,filter_size);
    hmm_input_output(i).delta_protein_sm = spot_vec - hmm_input_output(i).null_protein_sm;
    % detrend fluo
    fluo = hmm_input_output(i).fluo;
    activity = hmm_input_output(i).r_vec;
    time = round(hmm_input_output(i).time/60);
    for j = 1:numel(fluo)
        fluo(j) = fluo(j) - fluo_trend(time(j)==time_index);
        activity(j) = activity(j) - r_trend(time(j)==time_index);
    end
    hmm_input_output(i).fluo_dt = fluo;
    hmm_input_output(i).r_vec_dt = activity;
end


%%% Simple(ish) things first...
prctile_vec = 0:10:100;
spot_prctile_vec = NaN(size(prctile_vec));
null_prctile_vec = NaN(size(prctile_vec));
spot_protein = [hmm_input_output.spot_protein];
null_protein_smooth = [hmm_input_output.null_protein_sm];
fluo_vec = [hmm_input_output.fluo];
activity_vec = vertcat(hmm_input_output.r_vec)';

for i = 1:numel(prctile_vec)
    spot_prctile_vec(i) = prctile(spot_protein,prctile_vec(i));
    null_prctile_vec(i) = prctile(null_protein_smooth,prctile_vec(i));
end
[bin_counts, ~,~,binSpot,binNull] = histcounts2(spot_protein,null_protein_smooth,spot_prctile_vec,null_prctile_vec);
norm_counts = bin_counts ./ repmat(sum(bin_counts),numel(prctile_vec)-1,1);

resample_weights =  repmat(sum(norm_counts,2),1,numel(prctile_vec)-1) ./ norm_counts;

% now draw fluo, delta, and mf samples for each delta bin
null_boot_array = NaN(nBoots,numel(prctile_vec)-1);
spot_boot_array = NaN(nBoots,numel(prctile_vec)-1);
fluo_boot_array = NaN(nBoots,numel(prctile_vec)-1);
activity_boot_array = NaN(nBoots,numel(prctile_vec)-1);

for i = 1:numel(prctile_vec)-1
    bin_indices = find(binSpot==i);
    resample_wt_vec = resample_weights(:,i);
    sample_weights = resample_wt_vec(binNull(bin_indices));
    for n = 1:nBoots
        s_ids = randsample(bin_indices,numel(bin_indices),true,sample_weights);
        null_boot_array(n,i) = nanmean(null_protein_smooth(s_ids));
        spot_boot_array(n,i) = nanmean(spot_protein(s_ids));
        fluo_boot_array(n,i) = nanmean(fluo_vec(s_ids));
        activity_boot_array(n,i) = nanmean(activity_vec(s_ids));
    end
end

fluo_mean = nanmean(fluo_boot_array);
fluo_ste= nanstd(fluo_boot_array);

activity_mean = nanmean(activity_boot_array);
activity_ste= nanstd(activity_boot_array);

mf_mean = nanmean(null_boot_array);
mf_ste= nanstd(null_boot_array);

spot_mean = nanmean(spot_boot_array);
spot_ste= nanstd(spot_boot_array);

aspirational_fig = figure('Visible','off');
hold on
errorbar(spot_mean,activity_mean,activity_ste);
ylabel([gene_name ' activity']);
yyaxis right
plot(spot_mean,spot_mean)
hold on
plot(spot_mean,mf_mean)
ylabel([protein_name ' concentration (au)'])

legend([gene_name ' activity'], [protein_name ' (local)'], [protein_name ' (background)'],'Location','southeast')
xlabel(['local ' protein_name ' concentration (au)'])
saveas(aspirational_fig, [figPath 'aspirational_fig.png'])

%% Examine sna activity in vicinity of anamalous high and low pt points

spot_protein_dt_smooth = [hmm_input_output.spot_protein_dt_smooth];
null_protein_smooth = [hmm_input_output.null_protein_sm];

% calculate necessary sample weights
bulk_pt_high = prctile(spot_protein_dt_smooth,95);
bulk_pt_low = prctile(spot_protein_dt_smooth,5);
bulk_pt_bins = [min(spot_protein_dt_smooth),bulk_pt_low,bulk_pt_high,max(spot_protein_dt_smooth)];
% get distribution of bkg concentrations
[bin_counts, ~,~,binSpot,binNull] = histcounts2(spot_protein_dt_smooth,...
    null_protein_smooth,bulk_pt_bins,null_prctile_vec);
bin_probs = bin_counts ./ sum(bin_counts,2);
bin_weights = repmat(bin_probs(2,:),3,1) ./ bin_probs;

% take weighted sampels of activity signatures in vicinity of high, low,
% and middling pt points
binSpot(binSpot==0) = NaN;
binNull(binNull==0) = NaN;
nSamples = nansum(binSpot==1);
% initialize array
fluo_cell = cell(1,3);
local_protein_cell = cell(1,3);
null_protein_cell = cell(1,3);
activity_cell = cell(1,3);

% pull time series samples from vicinity of protein events
ind_vec = 1:3;
name_cell = {'high','other','low'};
ref_vec = -window_size:window_size;%1:(2*window_size+1);

for i = 1:numel(hmm_input_output)
    time_vec = hmm_input_output(i).time/60;
    frame_vec = 1:numel(time_vec);
    fluo_vec = hmm_input_output(i).fluo;
    activity_vec = hmm_input_output(i).r_vec;
    % find protein peaks and troughs
    pt_sp_dt_sm = hmm_input_output(i).spot_protein_dt_smooth;   
    pt_sp = hmm_input_output(i).spot_protein;   
    pt_nn_sm = hmm_input_output(i).null_protein_sm;   
    nan_filter = ~isnan(pt_nn_sm)&~isnan(pt_sp_dt_sm)&frame_vec>3&frame_vec<numel(frame_vec)-1;
    high_ids = find(nan_filter&pt_sp_dt_sm>bulk_pt_high);
    low_ids = find(nan_filter&pt_sp_dt_sm<bulk_pt_low);
    other_ids = find(nan_filter&~ismember(1:numel(pt_sp_dt_sm),[high_ids low_ids]));
    % sample proiten and activity traces
    for j = 1:numel(name_cell)
        eval(['ids =' name_cell{j} '_ids;'])
        fluo_temp = NaN(numel(ids),2*window_size+1);
        act_temp = NaN(numel(ids),2*window_size+1);
        loc_temp = NaN(numel(ids),2*window_size+1);
        bkg_temp = NaN(numel(ids),2*window_size+1);
        for k = 1:numel(ids)
            raw_ind = ref_vec+ids(k);
            ft_vec1 = raw_ind > 0 & raw_ind <= numel(fluo_vec);
            ft_vec2 = ref_vec(ft_vec1) + ids(k);
            % record
            fluo_temp(k,ft_vec1) = fluo_vec(ft_vec2);
            act_temp(k,ft_vec1) = activity_vec(ft_vec2);
            loc_temp(k,ft_vec1) = pt_sp_dt_sm(ft_vec2);
            bkg_temp(k,ft_vec1) = pt_nn_sm(ft_vec2);
            if isnan(bkg_temp(k,window_size+1))
                error('wtf')
            end
        end
        % add to main cell structures
        fluo_cell{j} = vertcat(fluo_cell{j}, fluo_temp);
        activity_cell{j} = vertcat(activity_cell{j}, act_temp);
        local_protein_cell{j} = vertcat(local_protein_cell{j}, loc_temp);
        null_protein_cell{j} = vertcat(null_protein_cell{j}, bkg_temp);
    end
end
%%
% now estimate bootstrap-averaged time trends
fluo_mean_trend_mat = NaN(2*window_size+1,3);
fluo_ste_trend_mat = NaN(2*window_size+1,3);

activity_mean_trend_mat = NaN(2*window_size+1,3);
activity_ste_trend_mat = NaN(2*window_size+1,3);

spot_protein_mean_trend_mat = NaN(2*window_size+1,3);
spot_protein_ste_trend_mat = NaN(2*window_size+1,3);

null_protein_mean_trend_mat = NaN(2*window_size+1,3);
null_protein_ste_trend_mat = NaN(2*window_size+1,3);

for i = 1:3
    % extract data
    fluo_raw = fluo_cell{i};
    act_raw = activity_cell{i};
    loc_raw = local_protein_cell{i};
    bkg_raw = null_protein_cell{i};    
    % initialize bootstrap arrays
    fluo_boot = NaN(nBoots,2*window_size+1);
    act_boot = NaN(nBoots,2*window_size+1);
    loc_boot = NaN(nBoots,2*window_size+1);
    bkg_boot = NaN(nBoots,2*window_size+1);
    % calculate weights
    bkg_vec = bkg_raw(:,window_size+1);
    [~,~,nullIDs] = histcounts(bkg_vec,null_prctile_vec);
    resample_wt_vec = bin_weights(i,:);
    resample_weights = resample_wt_vec(nullIDs);
    % conduct bootstrap sampling
    boot_index = 1:numel(bkg_vec);
    for j = 1:nBoots
        s_ids = randsample(boot_index,nSamples,true,resample_weights);
        fluo_boot(j,:) = nanmean(fluo_raw(s_ids,:));
        act_boot(j,:) = nanmean(act_raw(s_ids,:));
        loc_boot(j,:) = nanmean(loc_raw(s_ids,:));
        bkg_boot(j,:) = nanmean(bkg_raw(s_ids,:));
    end
    fluo_mean_trend_mat(:,i) = nanmean(fluo_boot);
    fluo_ste_trend_mat(:,i) = nanstd(fluo_boot);
    
    activity_mean_trend_mat(:,i) = nanmean(act_boot);
    activity_ste_trend_mat(:,i) = nanstd(act_boot);
    
    spot_protein_mean_trend_mat(:,i) = nanmean(loc_boot);
    spot_protein_ste_trend_mat(:,i) = nanstd(loc_boot);
    
    null_protein_mean_trend_mat(:,i) = nanmean(bkg_boot);
    null_protein_ste_trend_mat(:,i) = nanstd(bkg_boot);
end      




%%
trend_id_cell = {'burst start','burst stop','protein peak','protein dip','consistency check'};

metric_cell = {'spot_pt','spot_pt','r_vec','r_vec','spot_pt'}; % 

measure_cell = {'burst_rise_times','burst_fall_times','blip_times','dip_times','rand_ids'};

trend_index_cell = cell(1,numel(trend_id_cell));
trend_subindex_cell = cell(1,numel(trend_id_cell));
trend_linear_index_cell = cell(1,numel(trend_id_cell));


% calculate high and low thresholds for enrichment 
pt_prom = prctile([hmm_input_output.spot_protein_dt],75)-prctile([hmm_input_output.spot_protein_dt],25);

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
    spot_vec = hmm_input_output(i).spot_protein_dt;
    spot_vec_filtered = imgaussfilt(spot_vec,2);

    [~,blip_times,blip_widths] = findpeaks(spot_vec_filtered,frame_vec,'MinPeakWidth',60*.5/Tres,'MinPeakProminence',pt_prom); % NL: eye-balled atm
    [~,dip_times,dip_widths] = findpeaks(-spot_vec_filtered,frame_vec,'MinPeakWidth',60*.5/Tres,'MinPeakProminence',pt_prom); % NL: eye-balled atm
    
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
    end
end

% initialize data arrays tpo store lag pt info
lag_trend_array_cell = cell(1,numel(trend_id_cell));
for i = 1:numel(trend_id_cell)
    lag_trend_array_cell{i} = NaN(numel(trend_index_cell{i}),2*window_size+1);
end
% iterate through ids for each event type and store adjacent enrichment
% signatures
lag_index = -window_size:window_size;
ctrl_cell = {'trend'};
% pull trend samples
for n = 1:numel(ctrl_cell)
    ctrl_id =ctrl_cell{n};
    for i = 1:numel(trend_id_cell)
        eval(['lag_array = lag_' ctrl_id '_array_cell{i};']);
        eval(['id_vec = ' ctrl_id '_index_cell{i};']);
        eval(['sub_id_vec = ' ctrl_id '_subindex_cell{i};']);
        for j = 1:numel(id_vec)
            % protein
            spot_pt = hmm_input_output(id_vec(j)).spot_protein_dt;
            % transcription
            r_vec = hmm_input_output(id_vec(j)).r_vec;
            % index
            index_vec = lag_index + sub_id_vec(j);
            ft = index_vec>=1 & index_vec<=numel(spot_pt);
            % lag     
            eval(['lag_array(j,ft) = ' metric_cell{i} '(index_vec(ft));']);
        end
         eval(['lag_' ctrl_id '_array_cell{i} =  lag_array;']);
    end
end

trend_mean_array = NaN(2*window_size+1,numel(trend_id_cell));
trend_ste_array = NaN(2*window_size+1,numel(trend_id_cell));

% calculate bootstrap average and ste    
for i = 1:numel(trend_id_cell)
    trend_lag_array = lag_trend_array_cell{i};    
    % convenience v ec and bootstrap arrays
    index_vec = 1:size(trend_lag_array,1);
    trend_boot = NaN(nBoots,2*window_size+1);
    for j = 1:nBoots
        s_ids = randsample(index_vec,numel(index_vec),true);
        trend_boot(j,:) = nanmean(trend_lag_array(s_ids,:));
    end
    trend_mean_array(:,i) = nanmean(trend_boot);
    trend_ste_array(:,i) = nanstd(trend_boot);   
    
    % generate figure
    lag_fig = figure;
    hold on
    e1 = errorbar(lag_index*Tres,trend_mean_array(:,i),trend_ste_array(:,i));
    e1.CapSize = 0;
    legend([trend_id_cell{i} ' events'])
    xlabel('offset (seconds)')
%     ylabel(['absoulte ' protein_name ' enrichment (au)'])
    title([trend_id_cell{i} ' events'])
    saveas(lag_fig,[figPath trend_id_cell{i} '.png'])
end
    
