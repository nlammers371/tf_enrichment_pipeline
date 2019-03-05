% Script to probe relationship between input protein concentration and
% output transcriptional response
clear 
close all
% define ID variables
project = 'Dl_Venus_snaBAC_MCPmCherry_Zoom25x_minBleaching_test';
dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
dataPath = [dropboxFolder '\ProcessedEnrichmentData\' project '\'];
figPath = [dropboxFolder '\ProcessedEnrichmentFigures\' project '\input_output\'];
K = 2;
w = 6;
nTraces = 50; % number of individual traces to select for plotting

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

n_lags = 10;
n_bins = 10;
% load data set
% make burst vectors
rise_vec_full = [];
fall_vec_full = [];
r_vec_full = [];

particle_vec_full = [];

for i = 1:numel(hmm_input_output)
    rise_vec_full = [rise_vec_full [0 reshape(hmm_input_output(i).zz_mat(2,1,:),1,[])]];
    fall_vec_full = [fall_vec_full [0 reshape(hmm_input_output(i).zz_mat(1,2,:),1,[])]];
    r_vec_full = [r_vec_full reshape(sum(hmm_input_output(i).r_mat,2),1,[])];
    particle_vec_full = [particle_vec_full repelem(hmm_input_output(i).ParticleID,numel(hmm_input_output(i).time))];
end
time_vec_full = [hmm_input_output.time];
mf_vec = [hmm_input_output.mf_protein];
fluo_vec_full = [hmm_input_output.fluo];
delta_vec_full = [hmm_input_output.spot_protein] - [hmm_input_output.null_protein];

null_filter = isnan(delta_vec_full) | isnan(mf_vec);
time_vec_full = time_vec_full(~null_filter);
mf_vec = mf_vec(~null_filter);
delta_vec_full = delta_vec_full(~null_filter);
fluo_vec_full = fluo_vec_full(~null_filter);
rise_vec_full = rise_vec_full(~null_filter);
fall_vec_full = fall_vec_full(~null_filter);
r_vec_full = r_vec_full(~null_filter);
particle_vec_full = particle_vec_full(~null_filter);

set_vec_full = floor(particle_vec_full);

% identify outliers
low_ids = find(r_vec_full<prctile(r_vec_full,20));
high_ids = find(r_vec_full>prctile(r_vec_full,80));
rise_ids = find(rise_vec_full>.6);
fall_ids = find(fall_vec_full>.6);

% calculate sampling weights that allow us to draw control distributions
% that mimic mf and fluo values for subsets of interest
index_vec = 1:numel(fluo_vec_full);
% generate ID vector assigning each observation to a bin in 2D array
% bin_id_vec = NaN(size(fluo_vec));

% for i = 1:n_bins
%     for j = 1:n_bins
%         ft = fluo_vec < Xedges(i+1) & fluo_vec >= Xedges(i) & mf_vec < Yedges(j+1) & mf_vec >= Yedges(j);
%         bin_id_vec(ft) = sub2ind([n_bins,n_bins],j,i);
%     end
% end
% now calculate resampling weights fior each scenario
low_ft = ismember(index_vec,low_ids);
[baseHist,Xedges,Yedges,binX,binY] = histcounts2(fluo_vec_full(~low_ft),mf_vec(~low_ft),n_bins,'Normalization','probability');
bin_id_vec = sub2ind([n_bins,n_bins],binX,binY);
baseHist = baseHist + 1e-6;

lowHist = histcounts2(fluo_vec_full(low_ids),mf_vec(low_ids),Xedges,Yedges,'Normalization','probability');
low_wt = lowHist ./ baseHist;
low_wt_vec = low_wt(bin_id_vec);
low_ctrl_ids = randsample(index_vec(~low_ft),numel(low_ids),true,low_wt_vec);

% high
high_ft = ismember(index_vec,high_ids);
[baseHist,Xedges,Yedges,binX,binY] = histcounts2(fluo_vec_full(~high_ft),mf_vec(~high_ft),n_bins,'Normalization','probability');
bin_id_vec = sub2ind([n_bins,n_bins],binX,binY);
baseHist = baseHist + 1e-6;

highHist = histcounts2(fluo_vec_full(high_ids),mf_vec(high_ids),Xedges,Yedges,'Normalization','probability');
high_wt = highHist ./ baseHist;
high_wt_vec = high_wt(bin_id_vec);
high_ctrl_ids = randsample(index_vec(~high_ft),numel(high_ids),true,high_wt_vec);

% rise
rise_ft = ismember(index_vec,rise_ids);
[baseHist,Xedges,Yedges,binX,binY] = histcounts2(fluo_vec_full(~rise_ft),mf_vec(~rise_ft),n_bins,'Normalization','probability');
bin_id_vec = sub2ind([n_bins,n_bins],binX,binY);
baseHist = baseHist + 1e-6;

riseHist = histcounts2(fluo_vec_full(rise_ids),mf_vec(rise_ids),Xedges,Yedges,'Normalization','probability');
rise_wt = riseHist ./ baseHist;
rise_wt_vec = rise_wt(bin_id_vec);
rise_ctrl_ids = randsample(index_vec(~rise_ft),numel(rise_ids),true,rise_wt_vec);


% fall 
fall_ft = ismember(index_vec,fall_ids);
[baseHist,Xedges,Yedges,binX,binY] = histcounts2(fluo_vec_full(~fall_ft),mf_vec(~fall_ft),n_bins,'Normalization','probability');
bin_id_vec = sub2ind([n_bins,n_bins],binX,binY);
baseHist = baseHist + 1e-6;

fallHist = histcounts2(fluo_vec_full(fall_ids),mf_vec(fall_ids),Xedges,Yedges,'Normalization','probability');
fall_wt = fallHist ./ baseHist;
fall_wt_vec = fall_wt(bin_id_vec);
fall_ctrl_ids = randsample(index_vec(~fall_ft),numel(fall_ids),true,fall_wt_vec);

% low
low_array = NaN(numel(low_ids),n_lags+1);
low_ctrl_array = NaN(numel(low_ids),n_lags+1);
for i = 1:numel(low_ids)
    % sample
    ParticleID = particle_vec_full(low_ids(i));
    pt_vec = [NaN(1,max(0,n_lags-low_ids(i)+1)) particle_vec_full(max(1,low_ids(i)-n_lags):low_ids(i))];
    dt_vec = [NaN(1,max(0,n_lags-low_ctrl_ids(i)+1)) delta_vec_full(max(1,low_ctrl_ids(i)-n_lags):low_ctrl_ids(i))];
    low_array(i,pt_vec==ParticleID) = dt_vec(pt_vec==ParticleID);
    % control
    ParticleID = particle_vec_full(low_ctrl_ids(i));
    pt_vec = [NaN(1,max(0,n_lags-low_ctrl_ids(i)+1)) particle_vec_full(max(1,low_ctrl_ids(i)-n_lags):low_ctrl_ids(i))];
    dt_vec = [NaN(1,max(0,n_lags-low_ctrl_ids(i)+1)) delta_vec_full(max(1,low_ctrl_ids(i)-n_lags):low_ctrl_ids(i))];
    low_ctrl_array(i,pt_vec==ParticleID) = dt_vec(pt_vec==ParticleID);
end


high_array = NaN(numel(high_ids),n_lags+1);
high_ctrl_array = NaN(numel(high_ids),n_lags+1);
for i = 1:numel(high_ids)
    % sample
    ParticleID = particle_vec_full(high_ids(i));
    pt_vec = [NaN(1,max(0,n_lags-high_ids(i)+1)) particle_vec_full(max(1,high_ids(i)-n_lags):high_ids(i))];
    dt_vec = [NaN(1,max(0,n_lags-high_ctrl_ids(i)+1)) delta_vec_full(max(1,high_ctrl_ids(i)-n_lags):high_ctrl_ids(i))];
    high_array(i,pt_vec==ParticleID) = dt_vec(pt_vec==ParticleID);
    % control
    ParticleID = particle_vec_full(high_ctrl_ids(i));
    pt_vec = [NaN(1,max(0,n_lags-high_ctrl_ids(i)+1)) particle_vec_full(max(1,high_ctrl_ids(i)-n_lags):high_ctrl_ids(i))];
    dt_vec = [NaN(1,max(0,n_lags-high_ctrl_ids(i)+1)) delta_vec_full(max(1,high_ctrl_ids(i)-n_lags):high_ctrl_ids(i))];
    high_ctrl_array(i,pt_vec==ParticleID) = dt_vec(pt_vec==ParticleID);
end


rise_array = NaN(numel(rise_ids),n_lags+1);
rise_ctrl_array = NaN(numel(rise_ids),n_lags+1);

for i = 1:numel(rise_ids)
    % sample
    ParticleID = particle_vec_full(rise_ids(i));
    pt_vec = [NaN(1,max(0,n_lags-rise_ids(i)+1)) particle_vec_full(max(1,rise_ids(i)-n_lags):rise_ids(i))];
    dt_vec = [NaN(1,max(0,n_lags-rise_ctrl_ids(i)+1)) delta_vec_full(max(1,rise_ctrl_ids(i)-n_lags):rise_ctrl_ids(i))];
    rise_array(i,pt_vec==ParticleID) = dt_vec(pt_vec==ParticleID);
    % control
    ParticleID = particle_vec_full(rise_ctrl_ids(i));
    pt_vec = [NaN(1,max(0,n_lags-rise_ctrl_ids(i)+1)) particle_vec_full(max(1,rise_ctrl_ids(i)-n_lags):rise_ctrl_ids(i))];
    dt_vec = [NaN(1,max(0,n_lags-rise_ctrl_ids(i)+1)) delta_vec_full(max(1,rise_ctrl_ids(i)-n_lags):rise_ctrl_ids(i))];
    rise_ctrl_array(i,pt_vec==ParticleID) = dt_vec(pt_vec==ParticleID);
end


fall_array = NaN(numel(fall_ids),n_lags+1);
fall_ctrl_array = NaN(numel(fall_ids),n_lags+1);

for i = 1:numel(fall_ids)
    % sample
    ParticleID = particle_vec_full(fall_ids(i));
    pt_vec = [NaN(1,max(0,n_lags-fall_ids(i)+1)) particle_vec_full(max(1,fall_ids(i)-n_lags):fall_ids(i))];
    dt_vec = [NaN(1,max(0,n_lags-fall_ctrl_ids(i)+1)) delta_vec_full(max(1,fall_ctrl_ids(i)-n_lags):fall_ctrl_ids(i))];
    fall_array(i,pt_vec==ParticleID) = dt_vec(pt_vec==ParticleID);
    % control
    ParticleID = particle_vec_full(fall_ctrl_ids(i));
    pt_vec = [NaN(1,max(0,n_lags-fall_ctrl_ids(i)+1)) particle_vec_full(max(1,fall_ctrl_ids(i)-n_lags):fall_ctrl_ids(i))];
    dt_vec = [NaN(1,max(0,n_lags-fall_ctrl_ids(i)+1)) delta_vec_full(max(1,fall_ctrl_ids(i)-n_lags):fall_ctrl_ids(i))];
    fall_ctrl_array(i,pt_vec==ParticleID) = dt_vec(pt_vec==ParticleID);
end