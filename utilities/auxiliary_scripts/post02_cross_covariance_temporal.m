% Script to probe relationship between input protein concentration and
% output transcriptional response using cross-covariance
clear 
close all
% define ID variables
K = 3;
w = 7;
project = 'Dl-Ven x hbP2P';
nBoots = 100;
% project = 'Dl_Venus_hbP2P_MCPmCherry_Zoom2_7uW14uW';
dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
% dropboxFolder = 'C:\Users\nlamm\Dropbox (Garcia Lab)\';
dataPath = [dropboxFolder '\ProcessedEnrichmentData\' project '\'];
figPath = [dropboxFolder '\LocalEnrichmentFigures\' project '\hmm_input_output_K' num2str(K) '_w' num2str(w) '\'];
mkdir(figPath)
load([dataPath 'input_output_snips.mat'])
gene_name = 'snaBAC';
protein_name = 'Dorsal';

%%% Make time-dependent cross-covariance plots
% define some colors
yw = [234 194 100]/256; % yellow
bl = [115 143 193]/256; % blue
rd = [213 108 85]/256; % red
gr = [191 213 151]/256; % green
br = [207 178 147]/256; % brown

t_vec = [input_output_snips.t_center]; % center time for each snip
early_indices = find(t_vec < 20*60);
mid_indices = find(t_vec >= 60*20 & t_vec < 60*35);
late_indices = (t_vec >= 60*35);
n_lags = floor(numel(input_output_snips(3).hmm_ctrl_xcov)/2);

% conduct bootstrap sampling
early_xcov_mat = NaN(nBoots,2*n_lags+1);
mid_xcov_mat = NaN(nBoots,2*n_lags+1);
late_xcov_mat = NaN(nBoots,2*n_lags+1);
for n = 1:nBoots
    early_boot = randsample(early_indices,numel(early_indices),true);
    mid_boot = randsample(mid_indices,numel(mid_indices),true);
    late_boot = randsample(late_indices,numel(late_indices),true);

    early_xcov_mat(n,:) = nanmean(vertcat(input_output_snips(early_boot).hmm_spot_xcov));
    mid_xcov_mat(n,:) = nanmean(vertcat(input_output_snips(mid_boot).hmm_spot_xcov));
    late_xcov_mat(n,:) = nanmean(vertcat(input_output_snips(late_boot).hmm_spot_xcov));
end
early_xcov_mean = nanmean(early_xcov_mat);
early_xcov_ste = nanstd(early_xcov_mat);
mid_xcov_mean = nanmean(mid_xcov_mat);
mid_xcov_ste = nanstd(mid_xcov_mat);
late_xcov_mean = nanmean(late_xcov_mat);
late_xcov_ste = nanstd(late_xcov_mat);    

lag_axis = (-n_lags:n_lags)*20/60;
test_fig = figure;
hold on
fill([lag_axis fliplr(lag_axis)],[early_xcov_mean+early_xcov_ste fliplr(early_xcov_mean-early_xcov_ste)],bl,'FaceAlpha',.3,'EdgeAlpha',0)
fill([lag_axis fliplr(lag_axis)],[mid_xcov_mean+mid_xcov_ste fliplr(mid_xcov_mean-mid_xcov_ste)],gr,'FaceAlpha',.3,'EdgeAlpha',0)
fill([lag_axis fliplr(lag_axis)],[late_xcov_mean+late_xcov_ste fliplr(late_xcov_mean-late_xcov_ste)],rd,'FaceAlpha',.3,'EdgeAlpha',0)

p1 = plot(lag_axis,early_xcov_mean,'Color',bl,'LineWidth',1.3);
p2 = plot(lag_axis,mid_xcov_mean,'Color',gr,'LineWidth',1.3);
p3 = plot(lag_axis,late_xcov_mean,'Color',rd,'LineWidth',1.3);
legend([p1 p2 p3],'early','middle','late')
grid on
xlabel('offset (minutes)')
ylabel(['cross-covariance (' gene_name ' x ' protein_name])
saveas(test_fig,[figPath 'temporal_xcov.png'])