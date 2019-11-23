clear
close all
addpath('../utilities')
% specify paths
DropboxFolder =  'E:\Nick\LivemRNA\Dropbox (Personal)\';

project_cell = {'Dl-Ven_hbP2P-mCh_v2','Bcd-GFP_hbP2P-mCh_v2'};
master_struct = struct;
for  p = 1:numel(project_cell)
    [~, DataPath, FigRoot] =   header_function(DropboxFolder, project_cell{p});
    load([DataPath 'nucleus_struct.mat']);
    master_struct(p).nucleus_struct = nucleus_struct;
    master_struct(p).project = project_cell{p};
end
clear nucleus_struct;
FigPath = [FigRoot '/data_qc_checks/'];
mkdir(FigPath)

% generate indexing vectors
offset_vec = [];
fluo_vec = [];
time_vec = [];
set_vec = [];
ap_vec = [];

for m = 1:numel(master_struct)   
    nucleus_struct = master_struct(m).nucleus_struct;
    offset_temp = [nucleus_struct.fluoOffset];
    fluo_temp = [nucleus_struct.fluo];
    time_temp = [nucleus_struct.time];
    ap_temp = [nucleus_struct.APPosParticle];
    nan_ft = ~isnan(fluo_temp);
    offset_temp = offset_temp(nan_ft);
    fluo_temp = fluo_temp(nan_ft);
    time_temp = time_temp(nan_ft);   
    ap_temp = ap_temp(nan_ft);   
    set_temp = NaN(size(fluo_temp));
    iter = 1;
    for i = 1:numel(nucleus_struct)
        fluo = nucleus_struct(i).fluo;
        N = sum(~isnan(fluo));
        set_temp(iter:iter+N-1) = (m-1)*10 + repelem(nucleus_struct(i).setID,N);
        iter = iter + N;
    end
    % record
    offset_vec = [offset_vec offset_temp];
    ap_vec = [ap_vec ap_temp];
    fluo_vec = [fluo_vec fluo_temp];
    time_vec = [time_vec time_temp];
    set_vec = [set_vec set_temp];
end
clear nucleus_struct;
%%%
close all
%%% take bootstrap estimates of 99th percentile fluo and mean offset
nBoots = 100;
pct = 95;
set_index = unique(set_vec);
offset_array = NaN(nBoots,numel(set_index));
fluo_array = NaN(nBoots,numel(set_index));
time_ft = time_vec <= 1500 & time_vec >= 300;
% ap_ft = ap_vec <= 35 & time_vec >= 300;
% iterate    
for s = 1:numel(set_index)    
    set_ft = set_vec == set_index(s) & time_ft;% & mf_ft & time_ft;    
    index_vec = find(set_ft);
    N = numel(index_vec);
    n8_ids = floor((1:N)/N *100) == pct;
    for n = 1:nBoots
        boot_ids = randsample(index_vec,numel(index_vec),true);     
        fluo_boot = fluo_vec(boot_ids);
        offset_boot = offset_vec(boot_ids);
        
        [fluo_sorted, f_rank] = sort(fluo_boot);
        offset_sorted = offset_boot(f_rank);            
        offset_array(n,s) = nanmean(offset_sorted(n8_ids));
        fluo_array(n,s) = nanmean(fluo_sorted(n8_ids));
    end
end
% calculate average and standard error
offset_mean = nanmean(offset_array);
offset_ste = nanstd(offset_array);

fluo_mean = nanmean(fluo_array);
fluo_ste = nanstd(fluo_array);
    
%%% Make figure
% sna_filter = set_index<10|set_index>19;
close all
fig = figure;
cmap2 = brewermap([],'Set2');
hold on
errorbar(offset_mean,fluo_mean,-fluo_ste,fluo_ste,-offset_ste,offset_ste,'o','Color','black')
s1 = scatter(offset_mean(set_index<10),fluo_mean(set_index<10),'MarkerFaceColor',cmap2(2,:),'MarkerEdgeAlpha',0);
s2 = scatter(offset_mean(set_index>10),fluo_mean(set_index>10),'MarkerFaceColor',cmap2(3,:),'MarkerEdgeAlpha',0);
legend([s1 s2], 'OG line (homo)', 'OG line (het)','Location','southeast')
grid on
box on
xlabel('MCP offset')
ylabel('max spot fluorescence')
set(gca,'Fontsize',14)
ylim([0 280])
saveas(fig,[FigPath 'mcp_vs_fluo_hb_reporter.png'])    

%% look at max achieved fluo across full range of observed MCP levels
use_vec = true(size(set_vec));
offset_index = linspace(prctile(offset_vec(use_vec),1),prctile(offset_vec(use_vec),99),25);
off_window = 1;%*median(diff(offset_index));
fluo_max_array = NaN(nBoots,numel(offset_index));

for o = 1:numel(offset_index)
    ol = offset_index(max(1,o-off_window));
    ou = offset_index(min(numel(offset_index),o+off_window));
    off_ids = find(offset_vec>=ol & offset_vec <ou & use_vec);
    for n = 1:nBoots
        boot_ids = randsample(off_ids,numel(off_ids),true);
        fluo_max_array(n,o) = prctile(fluo_vec(boot_ids),pct);
    end
end

fluo_max_vec = nanmean(fluo_max_array);
fluo_max_vec_se = nanstd(fluo_max_array);

% make figure
boundary = 0.825;
close all
max_bound = figure;
hold on
e = errorbar(offset_index,fluo_max_vec,fluo_max_vec_se,'Color','black','LineWidth',1);
scatter(offset_index,fluo_max_vec,20,'MarkerFaceColor',cmap2(3,:),'MarkerEdgeColor','black')
plot([boundary boundary],[0 500],'--','Color','black')
e.CapSize = 0;
xlabel('MCP offset')
ylabel('maximum spot fluorescence')
set(gca,'Fontsize',14)
grid on
box on
% ylim([0 500])
% xlim([.3 1.5])
saveas(max_bound,[FigPath 'mcp_max_boundary.png'])   


%%
NaN_big = NaN(numel(set_vec),numel(unique(set_vec)));
label_vec = {'OG1 (sna)','OG2 (sna)','OG3 (sna)','OG4 (sna)','OG5 (sna)','OG6 (sna)','OG7 (sna)',...
    'OG1 (hb)','OG2 (hb)','OG3 (hb)','OG4 (hb)','eNosHom1','eNosHom2','eNosHom3','eNosHet1'};
for i = 1:numel(set_index)
    offsets = offset_vec(set_vec==set_index(i));
    NaN_big(1:numel(offsets),i) = offsets;
end

% close all
box_fig = figure;
hold on
bplot(NaN_big);
plot(0:numel(label_vec)+1,repelem(boundary,2+numel(label_vec)),'--','LineWidth',1.5,'Color','black')
set(gca,'xtick',1:1:numel(label_vec))
set(gca,'xticklabels',label_vec)
ylabel('MCP offset')
% xlabel('embryo')
set(gca,'Fontsize',10)
xtickangle(-45)
grid on
saveas(box_fig,[FigPath 'trouble_box_plot.png'])   