clear
close all
addpath('../utilities')
% specify paths
DropboxFolder =  'E:\Nick\LivemRNA\Dropbox (Personal)\';

project = 'Dl-Ven_snaBAC-mCh_v2';
master_struct = struct;
% load data
[~, DataPath, FigRoot] =   header_function(DropboxFolder, project);
load([DataPath 'nucleus_struct.mat']);
load([DataPath 'nucleus_struct_protein.mat']);
load([DropboxFolder 'ProcessedEnrichmentData\bleedthrough_controls\bleedthrough_model.mat'],'lm_full')
FigPath = [FigRoot '/mcp_saturation/'];
mkdir(FigPath)

%% generate indexing vectors

offset_vec = [nucleus_struct.fluoOffset];
fluo_vec = [nucleus_struct.fluo];
time_vec = [nucleus_struct.time];
nan_ft = ~isnan(fluo_vec);
offset_vec = offset_vec(nan_ft);
fluo_vec = fluo_vec(nan_ft);
time_vec = time_vec(nan_ft);    
dorsal_vec = [nucleus_struct_protein.mf_null_protein_vec];
set_vec = NaN(size(fluo_vec));
iter = 1;
for i = 1:numel(nucleus_struct)
    fluo = nucleus_struct(i).fluo;
    N = sum(~isnan(fluo));
    set_vec(iter:iter+N-1) = repelem(nucleus_struct(i).setID,N);
    iter = iter + N;
end

%% look at median fluo for different dorsal cohorts across full range of observed MCP levels

% first, generate adjusted mcp vector that accounts for Dl bleedthrough
venus_offset = lm_full.Coefficients.Estimate(1);
venus_slope = lm_full.Coefficients.Estimate(2);
offset_vec_adjusted = offset_vec - venus_offset - dorsal_vec*venus_slope;

% filter for time and remove set with incorrect time stamps
time_ft = time_vec >= 4*60 & time_vec<=20*60;
analysis_ft =  time_ft & set_vec~=4;

% assign to bins by according to dorsal concentration
dl_bins = 0:.2:1;
dorsal_ids = discretize(dorsal_vec,quantile(dorsal_vec(analysis_ft),dl_bins));

% get average fluoresence as a function of MCP for each dorsal cohort
nBoots = 100;
minDP = 50;
pct = 99;
offset_index = linspace(prctile(offset_vec_adjusted(analysis_ft),1),prctile(offset_vec_adjusted(analysis_ft),99),15);
offset_ids = discretize(offset_vec_adjusted,offset_index);
fluo_mean_array = NaN(numel(offset_index)-1,numel(dl_bins)-1,nBoots);
fluo_max_array = NaN(numel(offset_index)-1,numel(dl_bins)-1,nBoots);

% iterate
for d = 1:numel(dl_bins)-1
    for o = 1:numel(offset_index)-1
        index_vec = find(analysis_ft&offset_ids==o&dorsal_ids==d);
        if numel(index_vec) > minDP
            for n = 1:nBoots
                boot_ids = randsample(index_vec,numel(index_vec),true);
                fluo_mean_array(o,d,n) = nanmean(fluo_vec(boot_ids));
                fluo_max_array(o,d,n) = prctile(fluo_vec(boot_ids),100);
            end
        end
    end
end
    
% calculate bootstrap average and standard error
fluo_mean_mat = nanmean(fluo_mean_array,3);
fluo_mean_mat_se = nanstd(fluo_mean_array,[],3);

fluo_max_mat = nanmean(fluo_max_array,3);
fluo_max_mat_se = nanstd(fluo_max_array,[],3);


%% % make figures
offset_plot = offset_index(1:end-1) + mean(diff(offset_index))/2;
close all
mean_bound = figure;
cmap = brewermap(size(fluo_mean_mat,2),'Blues');
hold on
for s = 1:size(fluo_mean_mat,2)
    e = errorbar(offset_plot,fluo_mean_mat(:,s),fluo_mean_mat_se(:,s),'CapSize',0,'Color','black');
    scatter(offset_plot,fluo_mean_mat(:,s),'MarkerFaceColor',cmap(s,:),'MarkerEdgeColor','black');
end

e.CapSize = 0;
xlabel('MCP offset')
ylabel('average spot fluorescence')
set(gca,'Fontsize',14)
grid on
box on
% ylim([50 200])
xlim([.5 1.11])
saveas(mean_bound,[FigPath 'average_mcp_sat_dl_cohorts.png'])   


max_bound = figure;
cmap = brewermap(size(fluo_mean_mat,2),'Blues');
hold on
for s = 1:size(fluo_mean_mat,2)
    e = errorbar(offset_plot,fluo_max_mat(:,s),fluo_max_mat_se(:,s),'CapSize',0,'Color','black');
    scatter(offset_plot,fluo_max_mat(:,s),'MarkerFaceColor',cmap(s,:),'MarkerEdgeColor','black');
end

e.CapSize = 0;
xlabel('MCP offset')
ylabel('maximum spot fluorescence')
set(gca,'Fontsize',14)
grid on
box on
% ylim([50 200])
xlim([.5 1.11])
saveas(mean_bound,[FigPath 'max_mcp_sat_dl_cohorts.png'])   


%%
set_index = unique(set_vec);
NaN_big = NaN(numel(set_vec),numel(unique(set_vec)));
label_vec = {'OG1 (sna)','OG2 (sna)','OG3 (sna)','OG4 (sna)','OG5 (sna)','OG6 (sna)','OG7 (sna)'};
for i = 1:numel(set_index)
    offsets = offset_vec_adjusted(time_ft & set_vec == i);
    NaN_big(1:numel(offsets),i) = offsets;
end

% close all
box_fig = figure;
hold on
bplot(NaN_big);
% plot(0:numel(label_vec)+1,repelem(boundary,2+numel(label_vec)),'--','LineWidth',1.5,'Color','black')
set(gca,'xtick',1:1:numel(label_vec))
set(gca,'xticklabels',label_vec)
ylabel('MCP offset')
% xlabel('embryo')
set(gca,'Fontsize',10)
xtickangle(-45)
grid on
% saveas(box_fig,[FigPath 'trouble_box_plot.png'])   