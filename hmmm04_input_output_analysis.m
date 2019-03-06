% Script to probe relationship between input protein concentration and
% output transcriptional response
clear 
close all
% define ID variables
project = 'Dl_Venus_snaBAC_mCherry_Leica_hp';
dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
dataPath = [dropboxFolder '\ProcessedEnrichmentData\' project '\'];
figPath = [dropboxFolder '\LocalEnrichmentFigures\' project '\hmm_input_output\'];
K = 2;
w = 7;
nTraces = 50; % number of individual traces to select for plotting
window_size = 15; % number of lags over which to track protein/fluo dynamics
n_boots = 100;
n_ref_hist_bins = 10;
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
% ignore possibility of second project for noq 
hmm_input_output = master_struct(1).hmm_input_output;

% Calcualte average cross-correlation. Use score for randomly assigned
% trace pairs as a control
xc_index_vec = 1:numel(hmm_input_output);
xc_shuffle = randsample(xc_index_vec,numel(xc_index_vec),false);

xc_trend_mat = NaN(numel(xc_shuffle),2*window_size + 1);
xc_control_mat = NaN(numel(xc_shuffle),2*window_size + 1);
xc_trend_ct_mat = NaN(numel(xc_shuffle),2*window_size + 1);
xc_control_ct_mat = NaN(numel(xc_shuffle),2*window_size + 1);

for i = 1:numel(xc_shuffle)
    % compare true pairs first
    delta_pt = hmm_input_output(i).spot_protein;%-hmm_input_output(i).mf_protein;
    r_vec = sum(hmm_input_output(i).r_mat,2)';   
    if numel(delta_pt) > window_size + 1
        xc_trend_mat(i,:) = xcov(delta_pt,r_vec,window_size);
        xc_trend_ct_mat(i,:) = [numel(delta_pt)-window_size:numel(delta_pt) fliplr(numel(delta_pt)-window_size:numel(delta_pt)-1)];
    end
    
    % draw random activity vector and take cross-covariance
    r_vec_rand = sum(hmm_input_output(xc_shuffle(i)).r_mat,2)';
    n_max = min([numel(r_vec_rand) numel(delta_pt)]);
    if n_max >= window_size + 1
        xc_control_mat(i,:) = xcov(delta_pt(1:n_max),(r_vec_rand(1:n_max)),window_size);
        xc_control_ct_mat(i,:) = [n_max-window_size:n_max fliplr(n_max-window_size:n_max-1)];
    end
end

trend_boots = NaN(n_boots,size(xc_control_mat,2));
control_boots = NaN(n_boots,size(xc_control_mat,2));

trend_index_vec = find(~isnan(xc_trend_mat(:,1)));
control_index_vec = find(~isnan(xc_trend_mat(:,1)));
for i = 1:n_boots
    trend_ids = randsample(trend_index_vec,numel(trend_index_vec),true);
    trend_boots(i,:) = nansum(xc_trend_mat(trend_ids,:).*xc_trend_ct_mat(trend_ids,:)) ./ nansum(xc_trend_ct_mat(trend_ids,:));
    control_ids = randsample(control_index_vec,numel(control_index_vec),true);
    control_boots(i,:) = nansum(xc_control_mat(control_ids,:).*xc_control_ct_mat(control_ids,:)) ./ nansum(xc_control_ct_mat(control_ids,:));
end
xc_trend_mean = nanmean(trend_boots);
xc_trend_ste = nanstd(trend_boots);

xc_control_mean = nanmean(control_boots);
xc_control_ste = nanstd(control_boots);

% sketch method for extracting time resolution 
Tres = nanmedian(diff([hmm_input_output.time]));
lag_time_vec = (-window_size:0)*Tres;

% make figure
xc_fig = figure;
hold on
e1 = errorbar([lag_time_vec -fliplr(lag_time_vec(1:end-1))],xc_trend_mean,xc_trend_ste);
e1.CapSize = 0;
e2 = errorbar([lag_time_vec -fliplr(lag_time_vec(1:end-1))],xc_control_mean,xc_control_ste);
e2.CapSize = 0;
legend('active locus','random control')
xlabel('offset (seconds)')
ylabel('cross-covariance')
title(['cross-covariance between ' gene_name ' activity and ' protein_name ' enrichment'])
saveas(xc_fig,[figPath 'cross-covariance.png'])


% calculate prduction rate extrema
r_unit = nanmean(vertcat(hmm_input_output.r));
r_unit = r_unit(2);
mean_steps = 3;
% create indexing and storage vectors
timeBins = linspace(0,nanmax([hmm_input_output.time]),n_ref_hist_bins+1);
ptBins = linspace(0,nanmax([hmm_input_output.mf_protein]),n_ref_hist_bins+1);

trend_id_cell = {'burst start','burst stop','bulk rise','bulk fall','sustained high','sustained low'};
condition_cell = {'delta_r>.5*r_unit','delta_r<-.5*r_unit','bulk_delta_r>.2*r_unit',...
    'bulk_delta_r<-.2*r_unit','r_sum>mean_steps*.8*r_unit','r_sum<mean_steps*.2*r_unit'};

trend_index_cell = cell(1,numel(trend_id_cell));
trend_subindex_cell = cell(1,numel(trend_id_cell));

control_index_cell = cell(1,numel(trend_id_cell));
control_subindex_cell = cell(1,numel(trend_id_cell));

trend_weights = zeros(numel(timeBins)-1,numel(timeBins)-1,numel(trend_id_cell));
% generate longform ref indices as we go
ind_vec_full = [];
sub_ind_vec_full = [];
% collect rise and fall information
for i = 1:numel(hmm_input_output)
    % get isntantaneous activity state a take derivative
    r_vec = sum(hmm_input_output(i).r_mat,2);
    r_sum = conv(r_vec,ones(1,mean_steps),'same');
    delta_r = diff(r_vec);
    bulk_delta_r = conv(delta_r,ones(1,mean_steps),'same');
    % time info
    time_vec = round(hmm_input_output(i).time/60);
    mf_vec = hmm_input_output(i).mf_protein;
    
    ind_vec_full = [ind_vec_full repelem(i,numel(time_vec))];
    sub_ind_vec_full = [sub_ind_vec_full 1:numel(time_vec)];
    
    % find time steps taht meet individual conditions
    for j = 1:numel(condition_cell)  
        id_vec = trend_index_cell{j};
        sub_id_vec = trend_subindex_cell{j};        
        ids = eval(['find(' condition_cell{j} ')']);
        % record
        trend_index_cell{j} = [id_vec repelem(i,numel(ids))];
        trend_subindex_cell{j} = [sub_id_vec reshape(ids,1,[])];
        % record weights
        N = histcounts2(time_vec(ids),mf_vec(ids),timeBins,ptBins);
        trend_weights(:,:,j) = trend_weights(:,:,j) + N;        
    end
end

% now calculate sampling weights for control
trend_weights = trend_weights ./ sum(trend_weights(:));
sampling_weights = NaN(size(trend_weights));
% make ref vectors
time_ref = [hmm_input_output.time]/60;
mf_vec = [hmm_input_output.mf_protein];

for i = 1:numel(condition_cell)
    % require that control samples fall outside of set of interest
    ind_vec = trend_index_cell{i};
    subind_vec = trend_subindex_cell{i};
    t_temp = time_ref;
    t_temp(ismember(ind_vec_full,ind_vec)&ismember(sub_ind_vec_full,subind_vec)) = NaN;
    base_weights = histcounts(t_temp,[0 timeBins],'Normalization','probability')';
    sampling_weights(:,i) = trend_weights(:,i) ./ base_weights;
end
sampling_weights(isnan(sampling_weights)|isinf(sampling_weights)) = 0;

samp_wt_array = NaN(numel(sub_ind_vec_full),numel(condition_cell));
for i = 1:numel(condition_cell)
    for j = 1:numel(timeBins)
        samp_wt_array(round(time_ref)==timeBins(j),i) = sampling_weights(j,i);
    end    
end

% draw sample indices
for i = 1:numel(condition_cell)    
    rng(432); % seed for consistency
    ids = randsample(1:numel(ind_vec_full),numel(trend_index_cell{i}),true,samp_wt_array(:,i));
    control_index_cell{i} = ind_vec_full(ids);
    control_subindex_cell{i} = sub_ind_vec_full(ids);
end

% initialize data arrays tpo store lag pt info
lag_trend_array_cell = cell(1,numel(trend_id_cell));
lag_control_array_cell = cell(1,numel(trend_id_cell));
for i = 1:numel(trend_id_cell)
    lag_trend_array_cell{i} = NaN(numel(trend_index_cell{i}),window_size+1);
    lag_control_array_cell{i} = NaN(numel(control_index_cell{i}),window_size+1);
end
% iterate through ids for each event type and store preceding ewnrichment
% signatures
lag_index = 0:window_size;
% pull trend samples
for i = 1:numel(trend_id_cell)
    lag_array = lag_trend_array_cell{i};
    id_vec = trend_index_cell{i};
    sub_id_vec = trend_subindex_cell{i};
    for j = 1:numel(id_vec)
        % protein
        sp_pt = hmm_input_output(id_vec(j)).spot_protein;
        mf_pt = hmm_input_output(id_vec(j)).null_protein;
        delta = sp_pt-mf_pt;
        % lag
        delta = delta(max([1,sub_id_vec(j)-window_size]):sub_id_vec(j));
        lag_array(j,window_size+1-numel(delta)+1:end) = delta;
    end
    lag_trend_array_cell{i} = lag_array;
end

% pull control samples
for i = 1:numel(trend_id_cell)
    lag_array = lag_control_array_cell{i};
    id_vec = control_index_cell{i};
    sub_id_vec = control_subindex_cell{i};
    for j = 1:numel(id_vec)
        % protein
        delta_pt = hmm_input_output(id_vec(j)).spot_protein;
        mf_pt = hmm_input_output(id_vec(j)).mf_protein;
        delta = delta_pt-mf_pt;
        % lag
        delta = delta(max(1,sub_id_vec(j)-window_size):sub_id_vec(j));
        lag_array(j,window_size+1-numel(delta)+1:end) = delta;
    end
    lag_control_array_cell{i} = lag_array;
end

trend_mean_array = NaN(window_size+1,numel(trend_id_cell));
trend_ste_array = NaN(window_size+1,numel(trend_id_cell));
control_mean_array = NaN(window_size+1,numel(trend_id_cell));
control_ste_array = NaN(window_size+1,numel(trend_id_cell));

% calculate bootstrap average and ste    
for i = 1:numel(trend_id_cell)
    trend_lag_array = lag_trend_array_cell{i};
    control_lag_array = lag_control_array_cell{i};
    % convenience v ec and bootstrap arrays
    index_vec = 1:size(trend_lag_array,1);
    trend_boot = NaN(n_boots,window_size+1);
    control_boot = NaN(n_boots,window_size+1);
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
    e1 = errorbar(lag_time_vec,trend_mean_array(:,i),trend_ste_array(:,i));
    e1.CapSize = 0;
    e2 = errorbar(lag_time_vec,control_mean_array(:,i),control_ste_array(:,i));
    e2.CapSize = 0;
    legend([trend_id_cell{i} ' events'],'random control')
    xlabel('seconds before event')
    ylabel(['absoulte ' protein_name ' enrichment (au)'])
    title([trend_id_cell{i} ' events'])
    saveas(lag_fig,[figPath trend_id_cell{i} '.png'])
end
    
