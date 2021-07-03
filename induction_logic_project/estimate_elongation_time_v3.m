clear
close all

% set basic paths
DataRoot = 'C:\Users\nlamm\Dropbox (Personal)\InductionLogic\';
if ~exist(DataRoot)
  DataRoot = 'S:\Nick\Dropbox\InductionLogic\';
end

project = '20210430';
DataPath = [DataRoot filesep 'raw_data' filesep project filesep];

% load WT data
load([DataPath 'spot_struct.mat'])
Tres = spot_struct(1).Tres;

% set parameters for autocorr analysis
window_size = 15;
n_boots = 100;

%% calculate autocorrelation
particle_id_vec = [spot_struct.particleID];
set_id_vec = [spot_struct.setID];
set_index = unique(set_id_vec);
gene_id_vec = [spot_struct.geneID];
gene_index = unique(gene_id_vec);
fluo_vec = [spot_struct.fluo];
non_nan_ids = find(~isnan(particle_id_vec));

% initialize array
% trace_array = NaN(length(spot_struct(1).fluoInterp),length(non_nan_ids));
iter = 1;
low_thresh = 1e3;
sm_kernel_size = .5;

high_thresh_vec = [];
for s = 1:length(gene_index)
  fluo_vec = [spot_struct(gene_id_vec==gene_index(s)).fluo];
  high_thresh_vec(s) = prctile(fluo_vec,90);
end
%%
event_mat = NaN(1,2*window_size+1);
event_id_vec = NaN(1,1);
gene_id_vec2 = NaN(1,1);
genes_to_use = [];
id_vec = [];

for i = 1:length(spot_struct)
  
    trace = imgaussfilt(spot_struct(i).fluoInterp,sm_kernel_size);
    geneID = spot_struct(i).geneID;
    trace_len = length(trace);
    % identify low and high points
    low_filter = trace<=low_thresh & trace~=0;
    low_ids = find(low_filter);
    high_filter = trace>=high_thresh_vec(geneID);
    high_ids = find(high_filter);

    % iterate through trace to see if we can find any low->high or high->low
    % transitions
    id_vec = zeros(size(trace));
    id_vec(low_filter) = -1;
    id_vec(high_filter) = 1;

    event_ids = find(id_vec);

    for e = 1:length(event_ids)-1
        ind1 = event_ids(e);
        ind2 = event_ids(e+1);
        current_event = id_vec(ind1);
        next_event = id_vec(ind2);
        if ind2-ind1 <= 2*window_size && current_event~=next_event 
            % extract snip          
            snip_full = NaN(1,2*window_size+1);
            indices = ind1:min([trace_len,ind1+2*window_size]);
            snip_raw = trace(indices);
            snip_full(indices-ind1+1) = snip_raw;
            event_mat(end+1,:) = snip_full;
            event_id_vec(end+1) = next_event;
            gene_id_vec2(end+1) = geneID;
        end
    end
end

%% calculate weighted average

gene_index = unique(gene_id_vec(non_nan_ids));
autocorr_results = struct;
for g = 1:length(gene_index)
    filtered_traces = fragment_array(:,genes_to_use==gene_index(g));
    autocorr_results(g).filtered_traces = filtered_traces;
    [autocorr_results(g).wt_autocorr, autocorr_results(g).a_boot_errors, autocorr_results(g).wt_dd, autocorr_results(g).dd_boot_errors] = ...
        weighted_autocorrelation(filtered_traces, n_lags, 1,n_boots,ones(size(filtered_traces)));
end
%%
close all
bins = linspace(2e3,5e4);
figure;
hold on
histogram(autocorr_results(1).filtered_traces(:),bins,'Normalization','probability')
histogram(autocorr_results(2).filtered_traces(:),bins,'Normalization','probability')
histogram(autocorr_results(3).filtered_traces(:),bins,'Normalization','probability')
%% Looks like there's a ton of heterogeneity
slow_ids = find(rise_mat(6,:)>0.5);
fast_ids = find(rise_mat(6,:)<=0.5);

figure(1);
yyaxis left
plot(nanmean(rise_mat(:,fast_ids),2))

yyaxis right
plot(diff(nanmean(rise_mat(:,fast_ids),2),2))

