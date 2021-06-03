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
n_lags = 20;
n_boots = 100;

%% calculate autocorrelation
particle_id_vec = [spot_struct.particleID];
gene_id_vec = [spot_struct.geneID];
non_nan_ids = find(~isnan(particle_id_vec));

% initialize array
% trace_array = NaN(length(spot_struct(1).fluoInterp),length(non_nan_ids));
iter = 1;
% trace_array = vertcat(spot_struct(non_nan_ids).fluoInterp)';
sm_factor = 0.5;
min_dp = 25;
thresh = 2e3;
autocorr_mat = NaN(n_lags+1,1);
weight_mat = NaN(n_lags+1,1);
genes_to_use = [];
id_vec = [];
for i = 1:length(spot_struct)
  
  trace = spot_struct(i).fluoInterp;%imgaussfilt(spot_struct(i).fluoInterp,sm_factor);  
  
  % subdivide trace into contiguous active perods
  on_filter = trace>=thresh;
  start_id = find(on_filter,1);
  stop_id = find(on_filter,1,'last');    
  
  if stop_id - start_id + 1 >= min_dp 
    fragment = trace(start_id:stop_id);
    
    autocorr_mat(:,iter) = autocorr(fragment,n_lags);
    weight_mat(:,iter) = length(fragment):-1:length(fragment)-n_lags;
    
    fragment_array(1:length(fragment),iter) = fragment;
    id_vec(iter) = i;
    genes_to_use(iter) = gene_id_vec(i);
    iter = iter + 1;
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
slow_ids = find(autocorr_mat(6,:)>0.5);
fast_ids = find(autocorr_mat(6,:)<=0.5);

figure(1);
yyaxis left
plot(nanmean(autocorr_mat(:,fast_ids),2))

yyaxis right
plot(diff(nanmean(autocorr_mat(:,fast_ids),2),2))

