clear
close all
addpath('../utilities')
% specify paths
DropboxFolder =  'E:\Nick\LivemRNA\Dropbox (Personal)\';

project_cell = {'Dl-Ven_snaBAC-mCh','Dl-Ven_hbP2P-mCh'};
master_struct = struct;
for  p = 1:numel(project_cell)
    [~, DataPath, FigRoot] =   header_function(DropboxFolder, project_cell{p});
    load([DataPath 'nucleus_struct.mat']);
    master_struct(p).nucleus_struct = nucleus_struct;
    load([DataPath 'nucleus_struct_protein.mat']);
    master_struct(p).nucleus_struct_protein = nucleus_struct_protein;    
    master_struct(p).project = project_cell{p};
end
FigPath = [FigRoot '/data_qc_checks/'];
mkdir(FigPath)
clear nucleus_struct;
clear nucleus_struct_protein;
%% generate indexing vectors
offset_vec = [];
mcp_vec = [];
dorsal_vec = [];
fluo_vec = [];
time_vec = [];
set_vec = [];

for m = 1%:numel(master_struct)   
    % from nucleus_struct
    nucleus_struct = master_struct(m).nucleus_struct;
    offset_temp = [nucleus_struct.fluoOffset];
    fluo_temp = [nucleus_struct.fluo];
    time_temp = [nucleus_struct.time];
    nan_ft = ~isnan(fluo_temp);
    offset_temp = offset_temp(nan_ft);
    fluo_temp = fluo_temp(nan_ft);
    time_temp = time_temp(nan_ft);    
    % from protein strict
    nucleus_struct_protein = master_struct(m).nucleus_struct_protein;
    dorsal_temp = [nucleus_struct_protein.mf_null_protein_vec];
    mcp_temp = [nucleus_struct_protein.edge_null_mcp_vec];
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
    fluo_vec = [fluo_vec fluo_temp];
    time_vec = [time_vec time_temp];
    set_vec = [set_vec set_temp];
    mcp_vec = [mcp_vec mcp_temp];
    dorsal_vec = [dorsal_vec dorsal_temp];
end
%% simple cross-correlations
close all
% MCP sample vs offset
nan_ft = ~isnan(mcp_vec); 
p_mcp_off = polyfit(offset_vec(nan_ft),mcp_vec(nan_ft),1);
yfit_mcp_off = p_mcp_off(1)*offset_vec+p_mcp_off(2);
hold on;

dl_off_fig = figure;
hold on
scatter(offset_vec,mcp_vec)
plot(offset_vec,yfit_mcp_off,'r-.');
xlabel('MCP (random site)')
ylabel('MCP (inferred offset)')

% offset vs. Dorsal
nan_ft = ~isnan(dorsal_vec); 
p_dl_off = polyfit(dorsal_vec(nan_ft),offset_vec(nan_ft),1);
yfit_dl_off = p_dl_off(1)*dorsal_vec+p_dl_off(2);
hold on;

dl_off_fig = figure;
hold on
scatter(dorsal_vec,offset_vec)
plot(dorsal_vec,yfit_dl_off,'r-.');
xlabel('Dl level')
ylabel('MCP (inferred offset)')


% mcp vs. Dorsal
nan_ft = ~isnan(mcp_vec); 
p_dl_mcp = polyfit(dorsal_vec(nan_ft),mcp_vec(nan_ft),1);
yfit_dl_mcp = p_dl_mcp(1)*dorsal_vec+p_dl_mcp(2);
hold on;

dl_mcp_fig = figure;
hold on
scatter(dorsal_vec,mcp_vec)
plot(dorsal_vec,yfit_dl_mcp,'r-.');
xlabel('Dl level')
ylabel('MCP (random site)')


time_mcp_fig = figure;
hold on
scatter(time_vec,offset_vec)

%%

lm_time = fitlm([time_vec'],offset_vec')

lm_off = fitlm([time_vec' dorsal_vec'],offset_vec')



%%
%%%
close all
%%% take bootstrap estimates of 99th percentile fluo and mean offset
nBoots = 100;
pct = 100;
set_index = unique(set_vec);
offset_array = NaN(nBoots,numel(set_index));
fluo_array = NaN(nBoots,numel(set_index));
time_ft = time_vec <= 1200 & time_vec >= 300;
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
%% calculate average and standard error
offset_mean = nanmean(offset_array);
offset_ste = nanstd(offset_array);

fluo_mean = nanmean(fluo_array);
fluo_ste = nanstd(fluo_array);
    
%%% Make figure
sna_filter = set_index<10|set_index>19;
close all
fig = figure;
cmap2 = brewermap([],'Set2');
hold on
errorbar(offset_mean(sna_filter),fluo_mean(sna_filter),-fluo_ste(sna_filter),...
    fluo_ste(sna_filter),-offset_ste(sna_filter),offset_ste(sna_filter),'o','Color','black')
s1 = scatter(offset_mean(set_index<10),fluo_mean(set_index<10),'MarkerFaceColor',cmap2(2,:),'MarkerEdgeAlpha',0);
s2 = scatter(offset_mean(set_index>20&set_index<24),fluo_mean(set_index>20&set_index<24),'MarkerFaceColor',cmap2(3,:),'MarkerEdgeAlpha',0);
s3 = scatter(offset_mean(set_index==24),fluo_mean(set_index==24),'MarkerFaceColor',cmap2(4,:),'MarkerEdgeAlpha',0);
legend([s1 s2 s3], 'OG line', 'eNos2x (homo)', 'eNos2x (het)','Location','southeast')
grid on
box on
xlabel('MCP offset')
ylabel('max spot fluorescence')
set(gca,'Fontsize',14)
ylim([0 500])
saveas(fig,[FigPath 'mcp_vs_fluo.png'])    

%% look at max achieved fluo across full range of observed MCP levels
exclude_vec = ~(10<set_vec & set_vec < 20);
offset_index = linspace(prctile(offset_vec(exclude_vec),.1),prctile(offset_vec(exclude_vec),99),50);
off_window = 2;%*median(diff(offset_index));
fluo_max_array = NaN(nBoots,numel(offset_index));

for o = 1:numel(offset_index)
    ol = offset_index(max(1,o-off_window));
    ou = offset_index(min(numel(offset_index),o+off_window));
    off_ids = find(offset_vec>=ol & offset_vec <ou & exclude_vec);
    for n = 1:nBoots
        boot_ids = randsample(off_ids,numel(off_ids),true);
        fluo_max_array(n,o) = nanmax(fluo_vec(boot_ids));
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
ylim([0 500])
xlim([.3 1.5])
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
