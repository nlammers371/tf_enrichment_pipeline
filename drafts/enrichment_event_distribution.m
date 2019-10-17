% Script to investigate distribution of enrichment "event" sizes 
clear
close all
addpath('utilities')
% define core ID variables
projectTarget = 'Dl-Ven_snaBAC-mCh';
projectControl = 'Dl-Ven_hbP2P-mCh';
dropboxFolder = 'E:\Nick\LivemRNA\Dropbox (Personal)\';
dataPath = [dropboxFolder 'ProcessedEnrichmentData\'];
figPath = [dropboxFolder 'LocalEnrichmentFigures\enrichment_events\'];
mkdir(figPath)
% load data
w = 7;
K = 3;
load([dataPath projectTarget '\hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'],'hmm_input_output');
hmm_target = hmm_input_output;
load([dataPath projectControl '\hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'],'hmm_input_output');
hmm_control = hmm_input_output;
clear hmm_input_output

%% generate vectors of intensity changes and corresponding enrichment event size
% target
fluo_target = [hmm_target.fluo];
protein_delta_target = [hmm_target.spot_protein]-[hmm_target.null_protein];
mf_target = [hmm_target.mf_protein];

fluo_control = [hmm_control.fluo];
protein_delta_control = [hmm_control.spot_protein]-[hmm_control.null_protein];
mf_control = [hmm_control.mf_protein];

% analyze distributions
cm_hist = brewermap([],'Set2');

close all
pt_bins = linspace(-600,885,200);

overall_hist_pt = figure;
hold on
histogram(protein_delta_control,pt_bins,'Normalization','probability','FaceColor',cm_hist(3,:))
histogram(protein_delta_target,pt_bins,'Normalization','probability','FaceColor',cm_hist(2,:))
legend('enrichment (control)','enrichment (target)')
grid on
box on
set(gca,'Fontsize',14)
xlim([-500 500])
saveas(overall_hist_pt,[figPath 'full_enrichment_histograms.png'])
saveas(overall_hist_pt,[figPath 'full_enrichment_histograms.pdf'])

%%
% target
pd_target = fitdist(protein_delta_target','normal');
target_fit_line = pdf(pd_target,pt_bins');

target_hist_pt = figure;
hold on;
histogram(protein_delta_target,pt_bins,'Normalization','probability','FaceColor',cm_hist(2,:))
plot(pt_bins,target_fit_line/sum(target_fit_line),'Color','black','LineWidth',1)
grid on 
box on
xlabel('Dl enrichment')
ylabel('share')
set(gca,'Fontsize',14)
saveas(target_hist_pt,[figPath 'target_enrichment_fit.png'])
saveas(target_hist_pt,[figPath 'target_enrichment_fit.pdf'])

% control
pd_control = fitdist(protein_delta_control','normal');
control_fit_line = pdf(pd_control,pt_bins');

control_hist_pt = figure;
hold on;
histogram(protein_delta_control,pt_bins,'Normalization','probability','FaceColor',cm_hist(3,:))
plot(pt_bins,control_fit_line/sum(control_fit_line),'Color','black','LineWidth',1)
grid on 
box on
xlabel('Dl enrichment')
ylabel('share')
set(gca,'Fontsize',14)
saveas(control_hist_pt,[figPath 'control_enrichment_fit.png'])
saveas(control_hist_pt,[figPath 'control_enrichment_fit.pdf'])

%% Fit distributions to enrichment data, filtered by spot fluorescence
close all
q_vec = [.2 .4 .6 .8];
fluo_q_target = quantile(fluo_target,q_vec);
% generate separate vectors for low, middle, and high spots
pt_target_dim = protein_delta_target(fluo_target<=fluo_q_target(1));
pt_target_mid = protein_delta_target(fluo_target>fluo_q_target(2)&fluo_target<=fluo_q_target(3));
pt_target_bright = protein_delta_target(fluo_target>fluo_q_target(4));

pd_tg_dim = fitdist(pt_target_dim','normal')
pd_tg_mid = fitdist(pt_target_mid','normal')
pd_tg_bright = fitdist(pt_target_bright','normal')


cm_target = flipud(brewermap(3,'YlOrRd'));
pt_bins = linspace(-350,450);
test_fig_target = figure;
hold on
histogram(pt_target_bright,pt_bins,'Normalization','probability','FaceColor',cm_target(2,:))
% histogram(pt_delta_mid,pt_bins,'Normalization','probability','FaceColor',cm_target(2,:))
histogram(pt_target_dim,pt_bins,'Normalization','probability','FaceColor',cm_target(1,:))

fluo_q_control = quantile(fluo_delta_control,q_vec);
% generate separate vectors for low, middle, and high spots
pt_control_dim = protein_delta_control(fluo_delta_control<=fluo_q_control(1));
pt_control_mid = protein_delta_control(fluo_delta_control>fluo_q_control(2)&fluo_delta_control<=fluo_q_control(3));
pt_control_bright = protein_delta_control(fluo_delta_control>fluo_q_control(4));

pd_ct_dim = fitdist(pt_control_dim','normal')
pd_ct_mid = fitdist(pt_control_mid','normal')
pd_ct_bright = fitdist(pt_control_bright','normal')

test_fig_control = figure;
hold on
histogram(pt_control_bright,pt_bins,'Normalization','probability','FaceColor',cm_target(3,:))
% histogram(pt_delta_mid,pt_bins,'Normalization','probability','FaceColor',cm_target(2,:))
histogram(pt_control_dim,pt_bins,'Normalization','probability','FaceColor',cm_target(1,:))

%% 
mf_ub = 250;
mf_lb = 150;
pt_bins = linspace(-450,600,200);
mf_tg_ft = mf_target >= mf_lb & mf_target <= mf_ub;
mf_ct_ft = mf_control >= mf_lb & mf_control <= mf_ub;

delta_target = protein_delta_target(mf_tg_ft);
delta_control = protein_delta_control(mf_ct_ft);
% generate hypothetical distribution
mu_locus = nanmean(delta_target) - nanmean(delta_control);
sigma_locus = sqrt(nanvar(delta_target) - nanvar(delta_control));
logn_mu = log((mu_locus^2)/sqrt(sigma_locus^2+mu_locus^2))
logn_sigma = sqrt(log(sigma_locus^2/(mu_locus^2)+1))
locus_dist_sim = [normrnd(mu_locus,sigma_locus,1,numel(delta_target))];

delta_target_sim = randsample(delta_control,numel(delta_target),true) + locus_dist_sim;


hist_fig = figure;
hold on
histogram(delta_target,pt_bins,'Normalization','probability')
histogram(delta_target_sim,pt_bins,'Normalization','probability')
