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
% w = 7;
% K = 3;
% load([DataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'],'hmm_input_output');

% define size of window of interest
roi_window = 6; 
window_size = 15;
start = window_size + 2;
% extract roi_vectors from wapo and locus arrays
locus_protein_vec = nansum(results_struct.spot_array_dt(:,start:start + roi_window),2);% - ...
feature_sign_vec = results_struct.feature_sign_vec';
lag_size_vec = results_struct.lag_size_vec';
lead_size_vec = results_struct.lead_size_vec';
lag_dur_vec = results_struct.lag_dur_vec';
lead_dur_vec = results_struct.lead_dur_vec';
% tr_burst_size_vec = lag_dur_vec.*lag_size_vec;
mf_protein_vec = results_struct.mf_protein_vec;
% make rise filter
rise_ft = feature_sign_vec == 1;
analysis_ft = rise_ft & lead_dur_vec>5 & ~isnan(locus_protein_vec)&~isnan(lag_size_vec);

% generate vector of burst amplitudes
amp_val_vec = lag_size_vec(feature_sign_vec==1);
n_bins = 15;
amp_range = linspace(prctile(amp_val_vec,5),prctile(amp_val_vec,95),n_bins);

amp_sigma = median(diff(amp_range));
n_boots = 100;
amp_locus_pt_array = NaN(n_boots,numel(amp_range));
index_vec = find(analysis_ft);
for a = 1:numel(amp_range)
    for n = 1:n_boots
        s_ids = randsample(index_vec,numel(index_vec),true);
        boot_amps = lag_size_vec(s_ids);
        boot_protein = locus_protein_vec(s_ids);
        amp_curr = amp_range(a);
        wt_vec = exp(-.5*((boot_amps-amp_curr)/amp_sigma).^2);
        amp_locus_pt_array(n,a) = nansum(wt_vec.*boot_protein) ./ nansum(wt_vec);
    end
end
    
amp_pt_mean = nanmean(amp_locus_pt_array);
amp_pt_ste = nanstd(amp_locus_pt_array);
%
burst_amp_fig = figure;
hold on
e = errorbar(amp_range*3,amp_pt_mean,amp_pt_ste,'Color','black','LineWidth',1.5);
e.CapSize = 0;
scatter(amp_range*3,amp_pt_mean,75,'MarkerFaceColor',[213,108,85]/256,'MarkerEdgeColor','black')
% grid on
p = plot(0,0);
box on
xlim([.95*amp_range(1)*3  amp_range(end)*3+amp_range(1)*3-.95*amp_range(1)*3]);
% ylim([40 130])
xlabel('burst amplitude (au/min)')
ylabel('Dorsal enrichment (au)')
set(gca,'FontSize',14)
% set(gca,'Xtick',(3:3:12)/3)
% set(gca,'Ytick',(40:20:140))
StandardFigure(p,gca)
set(gca,'Color',[228 220 209]/255) 
saveas(burst_amp_fig,[FigPath 'burst_amp_surge_sz_fig.pdf'])
saveas(burst_amp_fig,[FigPath 'burst_amp_surge_sz_fig.png'])

%%
dur_range = 2:15;
burst_sigma = 2;
min_buffer_len = 5;
max_burst_dur = 30;


% initialize data arrays
dur_locus_pt_array = NaN(n_boots,numel(dur_range));
dur_sigma = 1; % 1 time step

for i = 1:numel(dur_range)    
    bd = dur_range(i);
    for n = 1:n_boots
        s_ids = randsample(index_vec,numel(index_vec),true);
        lag_boot_vec = lag_dur_vec(s_ids);
        boot_protein = locus_protein_vec(s_ids);
        burst_wts = exp(-.5*((lag_boot_vec-bd)/dur_sigma).^2);
        dur_locus_pt_array(n,i) = nansum(boot_protein.*burst_wts)/nansum(burst_wts);        
    end
end
dur_pt_mean = nanmean(dur_locus_pt_array);
dur_pt_ste = nanstd(dur_locus_pt_array);
    
%
burst_dur_fig = figure;
hm_cm = flipud(brewermap([],'RdYlBu'));
hold on
e = errorbar(dur_range/3,dur_pt_mean,dur_pt_ste,'Color','black','LineWidth',1.5);
e.CapSize = 0;
scatter(dur_range/3,dur_pt_mean,75,'MarkerFaceColor',[213,108,85]/256,'MarkerEdgeColor','black')
% grid on
p = plot(0,0);
box on
xlim([.95*dur_range(1)/3  dur_range(end)/3+dur_range(1)/3-.95*dur_range(1)/3]);
% ylim([40 130])
xlabel('burst duration (min)')
ylabel('Dorsal enrichment (au)')
set(gca,'FontSize',14)
% set(gca,'Xtick',(3:3:12)/3)
% set(gca,'Ytick',(40:20:140))
StandardFigure(p,gca)
set(gca,'Color',[228 220 209]/255) 
saveas(burst_dur_fig,[FigPath 'burst_dur_surge_sz_fig.pdf'])
saveas(burst_dur_fig,[FigPath 'burst_dur_surge_sz_fig.png'])

%% Regressions
% generate normalized regression vectors
% duration
burst_dur_ft = lag_dur_vec(analysis_ft);
burst_dur_ft = burst_dur_ft / nanstd(burst_dur_ft);
burst_dur_ft = burst_dur_ft - nanmean(burst_dur_ft);
% size
burst_size_ft = lag_size_vec(analysis_ft);
burst_size_ft = burst_size_ft / nanstd(burst_size_ft);
burst_size_ft = burst_size_ft - nanmean(burst_size_ft);
% nuclear protein
mf_protein_ft = mf_protein_vec(analysis_ft);
mf_protein_ft = mf_protein_ft / nanstd(mf_protein_ft);
mf_protein_ft = mf_protein_ft - nanmean(mf_protein_ft);
% average protein
locus_protein_ft = locus_protein_vec(analysis_ft);

mdl1 = fitlm(burst_dur_ft',burst_size_ft')

mdl2 = fitlm(burst_dur_ft',locus_protein_ft')

mdl3 = fitlm(burst_size_ft',locus_protein_ft')

mdl4 = fitlm(mf_protein_ft',locus_protein_ft)

mdl5 = fitlm([burst_dur_ft burst_size_ft mf_protein_ft'],locus_protein_ft)

