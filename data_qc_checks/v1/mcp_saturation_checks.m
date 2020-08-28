clear
close all
addpath('../utilities')
% specify paths
DropboxFolder =  'E:\Nick\LivemRNA\Dropbox (Personal)\';
project = 'Dl-Ven_snaBAC-mCh';
[~, DataPath, FigRoot] =   header_function(DropboxFolder, project);
FigPath = [FigRoot project '/data_qc_checks/'];
mkdir(FigPath)
% load data
load([DataPath 'nucleus_struct_protein.mat']);
load([DataPath 'nucleus_struct.mat']);
%% generate indexing vectors
offset_vec = [nucleus_struct.fluoOffset];
fluo_vec = [nucleus_struct.fluo];
nan_ft = ~isnan(fluo_vec);
offset_vec = offset_vec(nan_ft);
fluo_vec = fluo_vec(nan_ft);
time_vec = [nucleus_struct_protein.time];
dorsal_vec = [nucleus_struct_protein.mf_null_protein_vec];
set_vec = NaN(size(fluo_vec));
iter = 1;
for i = 1:numel(nucleus_struct)
    fluo = nucleus_struct(i).fluo;
    N = sum(~isnan(fluo));
    set_vec(iter:iter+N-1) = repelem(nucleus_struct(i).setID,N);
    iter = iter + N;
end
% set limits for [Dl] and time
mcp_time_bounds = [5 20]*60;
mf_bounds = [prctile(dorsal_vec,55),prctile(dorsal_vec,65)];


% look at time averages for offset in each dataset
time_ref = 1:50;
set_index = unique(set_vec);
offset_time = NaN(numel(time_ref),numel(set_index));
for s = 1:numel(set_index)
    for t = 1:numel(time_ref)
        offset_time(t,s) = nanmean(offset_vec(round(time_vec/60)==time_ref(t)&set_index(s)==set_vec));
    end
end


%
close all
%%% take bootstrap estimates of 95th percentile fluo and mean offset
nBoots = 100;
pct = 98;
offset_array = NaN(nBoots,numel(set_index));
fluo_array = NaN(nBoots,numel(set_index));

% filters
mf_ft = dorsal_vec >= mf_bounds(1) & dorsal_vec <= mf_bounds(2);
time_ft = time_vec >= mcp_time_bounds(1) & time_vec <= mcp_time_bounds(2);
% iterate    
for s = 1:numel(set_index)    
    set_ft = set_vec == set_index(s) & mf_ft & time_ft;    
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
mean_offset_mean = nanmean(offset_array);
mean_offset_ste = nanstd(offset_array);

max_fluo_mean = nanmean(fluo_array);
max_fluo_ste = nanstd(fluo_array);
    
%%% Make figure
fig = figure;
cmap2 = brewermap([],'Set2');
hold on
errorbar(mean_offset_mean,max_fluo_mean,-max_fluo_ste,max_fluo_ste,-mean_offset_ste,mean_offset_ste,'o','Color','black')
scatter(mean_offset_mean,max_fluo_mean,'MarkerFaceColor',cmap2(2,:));
grid on
box on
xlabel('MCP offset (95th percentile )')
ylabel('spot fluorescence (95th percentile )')
set(gca,'Fontsize',14)
% ylim([0 320])
saveas(fig,[FigPath 'mcp_vs_fluo_OG.png'])    

%% look at MCP offset vs Dl levels...is there a trend?
mf_lb = prctile(dorsal_vec,1);
mf_ub = prctile(dorsal_vec,90);
mf_index = linspace(mf_lb,mf_ub,50);

% generate protein vs time grid of offset values
offset_grid = NaN(numel(mf_index)-1,numel(time_ref)-1);
for t = 5:numel(time_ref)-1
    T1 = time_ref(max(t-1,1))*60;
    T2 = time_ref(min(numel(time_ref),t+2))*60;    
    t_ft = time_vec >= T1 & time_vec < T2;
    for m = 1:numel(mf_index)-1
        D1 = mf_index(max(m-1,1));
        D2 = mf_index(min(m+2,numel(mf_index)));
        m_ft = dorsal_vec >= D1 & dorsal_vec < D2;
        offset_grid(m,t) = nanmean(offset_vec(m_ft & t_ft));
    end
end

offset_hm = figure;
% cmap1 = flipud(brewermap([],'RdYlBu'));
% colormap(cmap1)
imagesc(offset_grid)
caxis([.7, 1.2])
xlabel('time into nc14')
ylabel('Dl-Venus intensity (au)')
h = colorbar;
ylabel(h,'MCP offset (au)')
set(gca,'Ytick',5:5:50,'yticklabels',round(mf_index(5:5:50)/10)*10)
set(gca,'xtick',5:5:50);%,'xticklabels',10:5:50)
set(gca,'Fontsize',14)
xlim([5 49])
ylim([10 49])
saveas(offset_hm,[FigPath 'mcp_vs_dl_and_time.png'])   

offset_plot = figure;
hold on
e = errorbar(mf_index(1:end-1),nanmean(offset_grid,2),nanstd(offset_grid,[],2),'Color','black');
e.CapSize = 0;
scatter(mf_index(1:end-1),nanmean(offset_grid,2),'MarkerFaceColor','black','MarkerEdgeAlpha',0)
ylabel('MCP offset (au)')
xlabel('Dl-Venus intensity (au)')
set(gca,'Fontsize',14)
grid on
saveas(offset_plot,[FigPath 'mcp_vs_dl.png'])   