% Script to attempt a systematic dissection of various factors driving
% proteinxtranscription burst coincidence
clear
close all
addpath('../utilities')
% define core ID variables
project = 'Dl-Ven_snaBAC-mCh_v3';
DropboxFolder =  'E:\Nick\LivemRNA\Dropbox (Personal)\';
[~, DataPath, FigRoot] =   header_function(DropboxFolder, project);
FigPath = [FigRoot '\' project '\burst_analyses\'];
mkdir(FigPath)
% load data
load([DataPath 'hmm_input_output_results.mat'])
w = 7;
K = 3;
load([DataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'],'hmm_input_output');

% define size of window of interest
roi_window = 6; 
window_size = 15;
start = window_size + 2;
% extract roi_vectors from wapo and locus arrays
locus_protein_vec = nansum(results_struct.spot_array_dt(:,start:start + roi_window),2);% - ...
%     nanmean(results_struct.spot_array_dt(:,start-2-roi_window:start-2),2);
% swap_protein_vec = nanmean(results_struct.swap_array_dt(:,start:start + roi_window),2);
% mf_protein_vec = nanmean(results_struct.mf_array(:,start:start + roi_window),2);
% pull other trend vectors
feature_sign_vec = results_struct.feature_sign_vec';
lag_size_vec = results_struct.lag_size_vec';
lead_size_vec = results_struct.lead_size_vec';
lag_dur_vec = results_struct.lag_dur_vec';
lead_dur_vec = results_struct.lead_dur_vec';
tr_burst_size_vec = lag_dur_vec.*lag_size_vec;
mf_protein_vec = results_struct.mf_protein_vec;
% make rise filter
rise_ft = feature_sign_vec == 1;
analysis_ft = rise_ft & lead_dur_vec>5 & ~isnan(locus_protein_vec)&~isnan(lag_size_vec);



burst_vec = (2:1:12)/3;
burst_sigma = 1;
n_boots = 100;
locus_pt_array = NaN(n_boots,numel(burst_vec));
index_vec = find(analysis_ft);
for b = 1:numel(burst_vec)
    for n = 1:n_boots
        s_ids = randsample(index_vec,numel(index_vec),true);
        boot_durs = lag_dur_vec(s_ids)/3;
        boot_protein = locus_protein_vec(s_ids);
        burst_curr = burst_vec(b);
        wt_vec = exp(-.5*((boot_durs-burst_curr)/burst_sigma).^2);
        locus_pt_array(n,b) = nansum(wt_vec.*boot_protein) ./ nansum(wt_vec);
    end
end
    
pt_mean = nanmean(locus_pt_array);
pt_ste = nanstd(locus_pt_array);
%%
dur_sz_fig = figure;
hm_cm = flipud(brewermap([],'RdYlBu'));
hold on
e = errorbar(burst_vec,pt_mean,pt_ste,'Color','black','LineWidth',1.5);
e.CapSize = 0;
scatter(burst_vec,pt_mean,75,'MarkerFaceColor',[213,108,85]/256,'MarkerEdgeColor','black')
% grid on
p = plot(0,0);
box on
xlim([.66 4]);
% ylim([40 130])
xlabel('burst duration (minutes)')
ylabel('Dorsal enrichment (au)')
set(gca,'FontSize',14)
set(gca,'Xtick',(3:3:12)/3)
% set(gca,'Ytick',(40:20:140))
StandardFigurePBoC(p,gca)
% saveas(dur_sz_fig,[FigPath 'burst_dur_surge_sz_fig.pdf'])
% saveas(dur_sz_fig,[FigPath 'burst_dur_surge_sz_fig.png'])

%%
burst_dur_ft = lag_dur_vec(analysis_ft);
burst_size_ft = lag_size_vec(analysis_ft);
locus_protein_ft = locus_protein_vec(analysis_ft);
mf_protein_ft = mf_protein_vec(analysis_ft);

mdl1 = fitlm(burst_dur_ft',burst_size_ft')

mdl2 = fitlm(burst_dur_ft',locus_protein_ft')

mdl3 = fitlm(burst_size_ft',locus_protein_ft')

mdl4 = fitlm(mf_protein_ft',locus_protein_ft)

mdl5 = fitlm([burst_dur_ft burst_size_ft mf_protein_ft'],locus_protein_ft)

