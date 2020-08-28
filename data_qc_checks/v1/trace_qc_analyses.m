clear
close all

% paths
DropboxFolder =  'E:\Nick\LivemRNA\Dropbox (Personal)\';
DataRoot = [DropboxFolder 'ProcessedEnrichmentData\'];
FigPath = [DropboxFolder 'LocalEnrichmentFigures\PipelineOutput\DataQCComparisons\'];
mkdir(FigPath);

% projects to compare
project_cell = {'Dl-Ven_snaBAC-mCh_v4','2xDl-Ven_snaBAC-mCh_v2','Dl-Ven_snaBAC-mCh_F-F-F'};
title_cell = {'original line','2x Dorsal','Dorsal CRISPR'};

% load sets and generate relevant data structures
master_struct = struct;

for p = 1:numel(project_cell)
    project = project_cell{p};
    load([DataRoot project '/nucleus_struct.mat'])
    % record basics
    master_struct(p).project = project;
    master_struct(p).nucleus_struct = nucleus_struct;
    set_index = unique([nucleus_struct.setID]);
    % record longform metrics
    protein_vec = vertcat(nucleus_struct.raw_nc_protein)';
    fluo_vec = [nucleus_struct.fluo];
    time_vec = [nucleus_struct.time];
    f_filter = ~isnan(fluo_vec);
    master_struct(p).fluo_vec = fluo_vec(f_filter);
    master_struct(p).time_vec = time_vec(f_filter);
    master_struct(p).time_vec_raw = time_vec(~isnan(protein_vec));
    master_struct(p).protein_vec = protein_vec(f_filter);
    master_struct(p).protein_vec_full = protein_vec(~isnan(protein_vec));
    % get trace QC info
    pt_indices = find(~isnan([nucleus_struct.ParticleID]));
    df_vec = [];
    df6_vec = [];
    nan_mean_vec = [];
    qc_iter = 0;
    for i = pt_indices
        fluo = nucleus_struct(i).fluo;
        fluo = fluo(find(~isnan(fluo),1):find(~isnan(fluo),1,'last'));
        df_vec = [df_vec diff(find(~isnan(fluo)))];
        nan_mean_vec = [nan_mean_vec nanmean(isnan(fluo))];
        
        if numel(fluo) >= 6
            ndf = diff(find(isnan([NaN fluo NaN])))-1;
            df6_vec = [df6_vec sum(ndf(ndf>6))/numel(fluo)];
        end
        
        qc_iter = qc_iter + 1*(sum(~isnan(fluo))>5);
    end
    master_struct(p).df_vec = df_vec;
    master_struct(p).df6_vec = df6_vec;
    master_struct(p).nan_mean_vec = nan_mean_vec;
    master_struct(p).traces_per_set = qc_iter/numel(unique([nucleus_struct.setID]));
    master_struct(p).qc_frac = qc_iter/numel(pt_indices);
end
   
%% Make trace QC plots
gap_bins = 0:30;
close all

gap_fig = figure;
hold on
for p = 1:numel(project_cell)
    histogram(master_struct(p).df_vec,gap_bins,'Normalization','probability')
end
legend(title_cell{:})
ylabel('share')
xlabel('gap (time steps)')
set(gca,'Fontsize',14)
xlim([min(gap_bins) max(gap_bins)])
set(gca,'Yscale','Log')
saveas(gap_fig,[FigPath 'gap_histogram.png'])

% fraction missing
nan_bins = linspace(0,1,50);
nan_fig = figure;
hold on
for p = 1:numel(project_cell)
    histogram(master_struct(p).nan_mean_vec,nan_bins,'Normalization','probability')
end
legend(title_cell{:})
ylabel('share')
xlabel('fraction of missing observations')
% set(gca,'Yscale','Log')
xlim([min(nan_bins) max(nan_bins)])
set(gca,'Fontsize',14)
ylim([0 .05])
saveas(nan_fig,[FigPath 'nan_histogram.png'])

% fraction of long runs
ctg_fig = figure;
hold on
for p = 1:numel(project_cell)
    histogram(master_struct(p).df6_vec,nan_bins,'Normalization','probability')
end
legend(title_cell{:},'Location','northwest')
ylabel('share')
xlabel('fraction of long (>2 minute) runs')
% set(gca,'Yscale','Log')
set(gca,'Fontsize',14)
ylim([0 .05])
xlim([min(nan_bins) max(nan_bins)])
saveas(ctg_fig,[FigPath 'longrun_histogram.png'])

%% bar plots

frac_qc_fig = figure;
hold on
for p = 1:numel(project_cell)
    barh(p,master_struct(p).qc_frac)
end
set(gca,'Ytick',1:3)
set(gca,'Yticklabel',title_cell)
xlabel('fraction of usable traces')
set(gca,'FontSize',14)
saveas(frac_qc_fig,[FigPath 'fraction_qc_bar.png'])

frac_qc_fig = figure;
hold on
for p = 1:numel(project_cell)
    barh(p,master_struct(p).traces_per_set)
end
set(gca,'Ytick',1:3)
set(gca,'Yticklabel',title_cell)
xlabel('number of usable traces')
set(gca,'FontSize',14)
saveas(frac_qc_fig,[FigPath 'number_qc_bar.png'])
% n_qc_fig = figure;
% barh([master_struct.traces_per_set])


%% input/output
close all
time_window = [0 20];
% raw fluorescence
fluo_bins = linspace(0,400);
fluo_fig = figure;
hold on
for p = 1:numel(project_cell)
    p_time = master_struct(p).time_vec/60;
    histogram(master_struct(p).fluo_vec(p_time >= time_window(1) & p_time <= time_window(2)),fluo_bins,'Normalization','probability')
end
legend(title_cell{:})
ylabel('share')
xlabel('spot fluorescence (au)')
set(gca,'Fontsize',14)
saveas(fluo_fig,[FigPath 'raw_fluo_hist.png'])

% raw Dorsal
dorsal_bins = linspace(0,5000);
dorsal_fig = figure;
hold on
for p = 1:numel(project_cell)
    p_time = master_struct(p).time_vec/60;
    histogram(master_struct(p).protein_vec(p_time >= time_window(1) & p_time <= time_window(2)),dorsal_bins,'Normalization','probability')
end
legend(title_cell{:})
ylabel('share')
xlabel('nuclear Dorsal concentration (au)')
set(gca,'Fontsize',14)
saveas(dorsal_fig,[FigPath 'raw_dl_hist.png'])


% set-wise dorsal histograms
color_cell = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250]};
for p = 1:numel(project_cell)
    dl_fig = figure;
    hold on 
    p_time = master_struct(p).time_vec/60;
    p_time_raw = master_struct(p).time_vec_raw/60;
    histogram(master_struct(p).protein_vec_full(p_time_raw >= time_window(1) & p_time_raw...
        <= time_window(2)),dorsal_bins,'FaceColor',color_cell{p}/1.25,'EdgeAlpha',0)
    ylabel('share')
    xlabel('nuclear Dorsal concentration (au)')
    set(gca,'Fontsize',14)
    saveas(dl_fig,[FigPath 'raw_dl_hist_' title_cell{p} '.png'])
    
    histogram(master_struct(p).protein_vec(p_time >= time_window(1) & p_time...
        <= time_window(2)),dorsal_bins,'FaceColor',color_cell{p}/.93,'EdgeAlpha',0)
    legend('all nuclei','nuclei with spots')
    saveas(dl_fig,[FigPath 'combined_dl_hist_' title_cell{p} '.png'])
end
    
    

% average input/output
dorsal_vec = [master_struct.protein_vec];
dorsal_bins = linspace(1000,5500);
dl_sigma = median(diff(dorsal_bins));
nBoots = 100;
% initialize arrays 
fluo_mean_array = NaN(numel(dorsal_bins),numel(project_cell),nBoots);
for p = 1:numel(project_cell)
    p_dorsal = master_struct(p).protein_vec;
    p_fluo = master_struct(p).fluo_vec;
    p_time = master_struct(p).time_vec/60;
    index_vec = find(p_time >= time_window(1) & p_time <= time_window(2));
    for n = 1:nBoots
        boot_indices = randsample(index_vec,numel(index_vec),true);
        b_dorsal = p_dorsal(boot_indices);
        b_fluo = p_fluo(boot_indices);
        for d = 1:numel(dorsal_bins)
            diff_vec = b_dorsal-dorsal_bins(d);
            if sum(abs(diff_vec) <= dl_sigma) >= 50
                wt_vec = exp(-.5*(diff_vec/dl_sigma).^2);
                fluo_mean_array(d,p,n) = nansum(wt_vec.*b_fluo) / nansum(wt_vec);
            end
        end
    end
end



mean_fluo = nanmean(fluo_mean_array,3);
ste_fluo = nanstd(fluo_mean_array,[],3);

in_out_fig = figure;
hold on
for p = 1:numel(project_cell)
    e = errorbar(dorsal_bins,mean_fluo(:,p),ste_fluo(:,p),'LineWidth',1.5);
    e.CapSize = 0;
end
e = errorbar(2*dorsal_bins,mean_fluo(:,1),ste_fluo(:,1),'--','Color',[0, 0.4470, 0.7410]);
e.CapSize = 0;
xlabel('Nuclear Dorsal concentration (au)')
ylabel('mean {\it snail} spot intensity (au)')
set(gca,'Fontsize',14);
xlim([800 5000])
grid on
legend([title_cell{:} {'orignal line (2x)'}],'Location','northwest')
saveas(in_out_fig,[FigPath 'input_output_comparison.png'])