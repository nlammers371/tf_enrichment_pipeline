% Script to probe relationship between input protein concentration and
% output transcriptional response
clear 
close all
% define ID variables
project = 'Dl_Venus_snaBAC_MCPmCherry_Zoom25x_minBleaching_test';
dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
dataPath = [dropboxFolder '\ProcessedEnrichmentData\' project '\'];
figPath = [dropboxFolder '\ProcessedEnrichmentFigures\' project '\hmm_input_output\'];
K = 2;
w = 6;
nTraces = 50; % number of individual traces to select for plotting
n_lags = 15; % number of lags over which to track protein/fluo dynamics
n_boots = 100;
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
        qc_fig = figure;
        hold on
        plot(time,fluo / nanmean(fluo))
        plot(time,fluo_check / nanmean(fluo_check));
        plot(time,mcp_check / nanmean(mcp_check))
        plot(time,r_vec / nanmean(r_vec))
        legend('fluo (HMM)', 'fluo (data)','raw mcp','activity state (HMM)')
        xlabel('time')
        ylabel([gene_name ' activity (au)'])
        saveas(qc_fig,[qcPath 'mcp_check_nc_' sprintf('%03d',plot_indieces(j)) '.png'])
        
        % Protein Channel checks
        spot_protein = hmm_input_output(plot_indices(j)).spot_protein;
        null_protein = hmm_input_output(plot_indices(j)).null_protein;
        mf_protein = hmm_input_output(plot_indices(j)).mf_protein;
        % make figure
        qc_fig = figure;
        hold on
        plot(time,spot_protein)
        plot(time,null_protein)
        plot(time,mf_protein)
        legend('protein (spot)', 'protein (control spot)','protein (mf control)')
        xlabel('time')
        ylabel([protein_name ' - ' protein_fluor ' (au)'])
        saveas(qc_fig,[qcPath 'protein_check_nc_' sprintf('%03d',plot_indieces(j)) '.png'])
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
        trace_fig = figure;
        hold on
        
        yyaxis left
        plot(time,r_vec)
        ylabel(['instantaneous ' gene_name ' activity (au)'])
        
        yyaxis right
        plot(time,delta_protein);
        ylabel(['instantaneous ' protein_name ' concentration (au)'])
                
        ax = gca;
        ax.YAxis(1).Color = 'black';
        ax.YAxis(2).Color = 'black';
        
        legend('transcriptional activity', 'local protein concentration')
        xlabel('time')
        ylabel([gene_name ' activity (au)'])
        saveas(trace_fig,[qcPath 'input_output_nc' sprintf('%03d',plot_indieces(j)) '.png'])               
    end
end

% Calcualte average cross-correlation. Use score for randomly assigned
% trace pairs as a control
xc_index_vec = 1:numel(hmm_input_output);
xc_shuffle = randsample(xc_index_vec,numel(xcov(sppt,flpt,n_xcorr_lags,'unbiased')),false);

xc_trend_mat = NaN(numel(xc_shuffle),2*n_lags + 1);
xc_control_mat = NaN(numel(xc_shuffle),2*n_lags + 1);

for i = 1:numel(xc_shuffle)
    % compare true pairs first
    delta_pt = hmm_input_output(i).spot_protein-hmm_input_output(i).mf_protein;
    r_vec = sum(hmm_input_output(i).r_mat,2)';   
    if numel(delta_pt) > n_lags + 1
        xc_trend_mat(i,:) = xcov(delta_pt,r_vec,n_lags,'unbiased');
    end
    
    % draw random activity vector and take cross-covariance
    r_vec_rand = sum(hmm_input_output(xc_shuffle(i)).r_mat,2)';
    n_max = min([numel(r_vec_rand) numel(delta_pt)]);
    if n_max >= n_lags + 1
        xc_control_mat(i,:) = xcov(delta_pt(1:n_max),r_vec_rand(1:n_max),n_lags,'unbiased');
    end
end

trend_boots = NaN(n_boots,size(xc_control_mat,2));
control_boots = NaN(n_boots,size(xc_control_mat,2));

trend_index_vec = find(~isnan(xc_trend_mat(:,1)));
control_index_vec = find(~isnan(xc_trend_mat(:,1)));
for i = 1:n_boots
    trend_boots(i,:) = nanmean(xc_trend_mat(randsample(trend_index_vec,numel(trend_index_vec),true),:));
    control_boots(i,:) = nanmean(xc_control_mat(randsample(control_index_vec,numel(control_index_vec),true),:));
end
xc_trend_mean = nanmean(trend_boots);
xc_trend_ste = nanstd(trend_boots);

xc_control_mean = nanmean(control_boots);
xc_control_ste = nanstd(control_boots);

% sketch method for extracting time resolution 
Tres = nanmedian([hmm_input_output.time]);
lag_time_time_vec = (-n_lags:0)*Tres;

% make figure
xc_fig = figure;
hold on
e1 = errorbar(lag_time_vec,trend_mean_array(:,i),trend_ste_array(:,i));
e1.CapSize = 0;
e2 = errorbar(lag_time_vec,control_mean_array(:,i),control_ste_array(:,i));
e2.CapSize = 0;
legend('active locus','random control')
xlabel('offset (seconds)')
ylabel('cross-covariance')
title(['cross-covariance between ' gene_name ' activity and ' protein_name ' enrichment'])
saveas(xc_fig,[figPath 'cross-covariance.png'])


% calculate prduction rate extrema
r_unit = nanmean([hmm_input_output.r_inf ; ]);
r_unit = r_unit(2);
mean_steps = 5;
% create indexing and storage vectors
time_index = 1:nanmax(round([hmm_input_output.time]/60));
trend_id_cell = {'burst start','burst stop','bulk rise','bulk fall','sustained high','sustained low'};
condition_cell = {'delta_z>0','delta_z<0','bulk_delta_r>r_unit',...
    'bulk_delta_r<-r_unit','r_sum>mean_steps*.8*r_unit','r_sum<mean_steps*.2*r_unit'};

trend_index_cell = cell(1,numel(trend_id_cell));
trend_subindex_cell = cell(1,numel(trend_id_cell));

control_index_cell = cell(1,numel(trend_id_cell));
control_subindex_cell = cell(1,numel(trend_id_cell));

time_weights = zeros(numel(time_index),numel(trend_id_cell));
% generate longform ref indices as we go
ind_vec_full = [];
sub_ind_vec_full = [];
% collect rise and fall information
for i = 1:numel(hmm_input_output)
    % calculate most likely state and take derivative
    [~, z_vec] = max(hmm_input_output(i).z_mat,[],2);
    delta_z = [0 diff(z_vec)];
    % get isntantaneous activity state a take derivative
    r_vec = sum(hmm_input_output(i).r_mat,2);
    r_sum = conv(r_vec,ones(1,mean_steps),'same');
    delta_r = diff(r_vec);
    bulk_delta_r = conv(delta_r,ones(1,mean_steps),'same');
    % time info
    time_vec = round(hmm_input_output(i).time/60);
    
    ind_vec_full = [ind_vec_full repelem(i,numel(time_vec))];
    sub_ind_vec_full = [sub_ind_vec_full 1:numel(time_vec)];
    
    for j = 1:numel(condition_cell)  
        id_vec = trend_index_cell{j};
        sub_id_vec = trend_subindex_cell{j};
        ids = eval(['find(' condition_cell{j} ')']);
        trend_index_cell{j} = [id_vec repelem(i,numel(ids))];
        trend_subindex_cell{j} = [sub_id_vec ids];
        time_weights(:,j) = time_weights(:,j) + histc(time_vec,time_index);
    end
end

% calculate sampling weights for control
time_weights = time_weights ./ sum(time_weights);
time_ref = round([hmm_input_output.time]/60);
base_weights = histcounts(round([hmm_input_output.time]/60),time_index,'Normalization','probability');
sampling_weights = time_weights ./ base_weights;
samp_wt_vec = NaN(size(sub_ind_vec_full));
for i = 1:numel(time_index)
    samp_wt_vec(time_ref==time_index(i)) = sampling_weights(i);
end
% draw sample indices
for i = 1:numel(condition_cell)    
    rng(432); % seed for consistency
    ids = randsample(1:numel(ind_vec_full),numel(trend_index_cell{i}),true,samp_wt_vec);
    control_index_cell{i} = ind_vec_full(ids);
    control_subindex_cell{i} = sub_ind_vec_full(ids);
end

% initialize data arrays tpo store lag pt info
lag_trend_array_cell = cell(1,numel(trend_id_cell));
lag_control_array_cell = cell(1,numel(trend_id_cell));
for i = 1:numel(trend_id_cell)
    lag_trend_array_cell{i} = NaN(numel(trend_index_cell{i}),n_lags+1);
    lag_control_array_cell{i} = NaN(numel(control_index_cell{i}),n_lags+1);
end
% iterate through ids for each event type and store preceding ewnrichment
% signatures
lag_index = 0:n_lags;
% pull trend samples
for i = 1:numel(trend_id_cell)
    lag_array = lag_trend_array_cell{i};
    id_vec = trend_index_cell{j};
    sub_id_vec = trend_subindex_cell{j};
    for j = 1:numel(id_vec)
        % protein
        delta_pt = hmm_input_output(id_vec(j)).spot_protein;
        mf_pt = hmm_input_output(id_vec(j)).mf_protein;
        delta = delta_pt-mf_pt;
        % lag
        delta = delta(max(1,sub_id_vec(j)-n_lags),sub_id_vec(j));
        lag_array(j,n_lags+1-numel(delta)+1:end) = delta;
    end
    lag_trend_array_cell{i} = lag_array;
end

% pull control samples
for i = 1:numel(trend_id_cell)
    lag_array = lag_control_array_cell{i};
    id_vec = control_index_cell{j};
    sub_id_vec = control_subindex_cell{j};
    for j = 1:numel(id_vec)
        % protein
        delta_pt = hmm_input_output(id_vec(j)).spot_protein;
        mf_pt = hmm_input_output(id_vec(j)).mf_protein;
        delta = delta_pt-mf_pt;
        % lag
        delta = delta(max(1,sub_id_vec(j)-n_lags),sub_id_vec(j));
        lag_array(j,n_lags+1-numel(delta)+1:end) = delta;
    end
    lag_control_array_cell{i} = lag_array;
end

trend_mean_array = NaN(n_lags+1,numel(trend_ids_cell));
trend_ste_array = NaN(n_lags+1,numel(trend_ids_cell));
control_mean_array = NaN(n_lags+1,numel(trend_ids_cell));
control_ste_array = NaN(n_lags+1,numel(trend_ids_cell));

% calculate bootstrap average and ste    
for i = 1:numel(trend_id_cell)
    trend_lag_array = lag_trend_array_cell{i};
    control_lag_array = lag_control_array_cell{i};
    % convenience v ec and bootstrap arrays
    index_vec = 1:size(lag_trend_array,1);
    trend_boot = NaN(n_boots,n_lags+1);
    control_boot = NaN(n_boots,n_lags+1);
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
    
