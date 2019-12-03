clear
close all
addpath('../utilities')
% specify paths
DropboxFolder =  'E:\Nick\LivemRNA\Dropbox (Personal)\';

project_cell = {'Dl-Ven_snaBAC-mCh_v2','Bcd-GFP_snaBAC-mCh_v2'};
master_struct = struct;
for  p = 1:numel(project_cell)
    [~, DataPath, FigRoot] =   header_function(DropboxFolder, project_cell{p});
    load([DataPath 'nucleus_struct.mat']);
    master_struct(p).nucleus_struct = nucleus_struct;
    master_struct(p).project = project_cell{p};
end
FigPath = [FigRoot '/mcp_saturation/'];
mkdir(FigPath)

%%% generate indexing vectors
offset_vec = [];
fluo_vec = [];
time_vec = [];
set_vec = [];

for m = 1:numel(master_struct)   
    nucleus_struct = master_struct(m).nucleus_struct;
    offset_temp = [nucleus_struct.fluoOffset];
    fluo_temp = [nucleus_struct.fluo];
    time_temp = [nucleus_struct.time];
    nan_ft = ~isnan(fluo_temp);
    offset_temp = offset_temp(nan_ft);
    fluo_temp = fluo_temp(nan_ft);
    time_temp = time_temp(nan_ft);    
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
end
clear nucleus_struct;
%%%
close all
%%% take bootstrap estimates of 99th percentile fluo and mean offset
nBoots = 100;
pct = 99;
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
sna_filter = true(size(set_index));
close all
fig = figure;
cmap2 = brewermap([],'Set2');
hold on
errorbar(offset_mean(sna_filter),fluo_mean(sna_filter),-fluo_ste(sna_filter),...
    fluo_ste(sna_filter),-offset_ste(sna_filter),offset_ste(sna_filter),'o','Color','black')

s1 = scatter(offset_mean(set_index<10),fluo_mean(set_index<10),'MarkerFaceColor',cmap2(2,:),'MarkerEdgeAlpha',0);
s2 = scatter(offset_mean(set_index>9),fluo_mean(set_index>9),'MarkerFaceColor',cmap2(3,:),'MarkerEdgeAlpha',0);
legend([s1 s2], 'OG line', 'OG Trans-Het','Location','southeast')
grid on
box on
xlabel('MCP offset')
ylabel('spot intensity (99th percentile)')
set(gca,'Fontsize',14)
ylim([0 500])
saveas(fig,[FigPath 'sna_mcp_vs_fluo.png'])    

