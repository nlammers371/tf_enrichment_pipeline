clear
close all

% set basic paths
DataRoot = 'C:\Users\nlamm\Dropbox (Personal)\InductionLogic\';

project = '20200807_WT';
DataPath = [DataRoot project filesep];

% load WT data
load([DataPath 'spot_struct.mat'])
Tres = spot_struct(1).Tres;

% set parameters for autocorr analysis
n_lags = 20;
n_boots = 100;

% % %% calculate autocorrelation
% % [wt_autocorr, a_boot_errors, wt_dd, dd_boot_errors] = ...
% %     weighted_autocorrelation(wt_array, n_lags, 1,n_boots,ones(size(wt_array)));

% employ mild smoothing filter to reduce noisy fluctuations
% wt_traces_filtered = wt_array;
% for i = 1:size(wt_array,2)
%   wt_traces_filtered(:,i) = imgaussfilt(wt_traces_filtered(:,i),.1);
% end

%% calculate autocorrelation
min_dp = 50;
autocorr_mat = NaN(n_lags+1,1);
weight_mat = NaN(n_lags+1,1);
iter = 1;

for i = 1:length(spot_struct)
  
  trace = spot_struct(i).fluo;  
  
  % subdivide trace into contiguous active perods
  off_ids = find(trace==0);
  contig_fragments = diff(off_ids);  
  long_fragments = find(contig_fragments > n_lags);
  
  for j = 1:length(long_fragments)    
    fragment = trace(off_ids(long_fragments(j)):off_ids(long_fragments(j)+1));
    
    autocorr_mat(:,iter) = autocorr(fragment,n_lags);
    weight_mat(:,iter) = length(fragment):-1:length(fragment)-n_lags;
    
    fragment_array(1:length(fragment),iter) = fragment;
    
    iter = iter + 1;
  end

end

%% calculate weighted average

autocorr_boot_array = NaN(n_lags+1,n_boots);
autocorr_dd_boot_array = NaN(n_lags-1,n_boots);
sample_index = 1:size(autocorr_mat,2);
for n = 1:n_boots
  boot_indices = randsample(sample_index,length(sample_index),true);
  autocorr_boot_array(:,n) = nansum(autocorr_mat(:,boot_indices).*weight_mat(:,boot_indices),2) ./ nansum(weight_mat(:,boot_indices),2);
  autocorr_dd_boot_array(:,n) = diff(autocorr_boot_array(:,n),2);
end

autocorr_mean = nanmean(autocorr_boot_array,2);
autocorr_ste = nanstd(autocorr_boot_array,[],2);
autocorr_dd_mean = nanmean(autocorr_dd_boot_array,2);
autocorr_dd_ste = nanstd(autocorr_dd_boot_array,[],2);

%% Looks like there's a ton of heterogeneity
slow_ids = find(autocorr_mat(6,:)>0.5);
fast_ids = find(autocorr_mat(6,:)<=0.5);

figure(1);
yyaxis left
plot(nanmean(autocorr_mat(:,fast_ids),2))

yyaxis right
plot(diff(nanmean(autocorr_mat(:,fast_ids),2),2))

