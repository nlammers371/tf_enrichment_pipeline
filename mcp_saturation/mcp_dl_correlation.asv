clear
close all
addpath('../utilities')
% specify paths
DropboxFolder =  'E:\Nick\LivemRNA\Dropbox (Personal)\';

project_cell = {'Dl-Ven_snaBAC-mCh_v2','Dl-Ven_hbP2P-mCh_v2'};
master_struct = struct;
for  p = 1:numel(project_cell)
    [~, DataPath, FigRoot] =   header_function(DropboxFolder, project_cell{p});
%     load([DataPath 'nucleus_struct.mat']);
%     master_struct(p).nucleus_struct = nucleus_struct;
    load([DataPath 'nucleus_struct_protein.mat']);
    master_struct(p).nucleus_struct_protein = nucleus_struct_protein;    
    master_struct(p).project = project_cell{p};
end
FigPath = [FigRoot '/mcp_offset_analyses/'];
mkdir(FigPath)
clear nucleus_struct;
clear nucleus_struct_protein;

%%% Generate indexing vectors
offset_vec = [];
mcp_vec = [];
dorsal_vec = [];
fluo_vec = [];
time_vec = [];
set_vec = [];

for m = 1:numel(master_struct)   
    % from nucleus_struct
    nucleus_struct_protein = master_struct(m).nucleus_struct_protein;
    offset_temp = [nucleus_struct_protein.fluoOffset];
    fluo_temp = [nucleus_struct_protein.fluo];
    time_temp = [nucleus_struct_protein.time];
%     nan_ft = ~isnan(fluo_temp);
%     offset_temp = offset_temp(nan_ft);
%     fluo_temp = fluo_temp(nan_ft);
%     time_temp = time_temp(nan_ft);       
    dorsal_temp = [nucleus_struct_protein.mf_null_protein_vec];
%     mcp_temp = [nucleus_struct_protein.edge_null_mcp_vec];
    set_temp = NaN(size(fluo_temp));
    iter = 1;
    for i = 1:numel(nucleus_struct_protein)
        fluo = nucleus_struct_protein(i).fluo;
        N = sum(~isnan(fluo));
        set_temp(iter:iter+N-1) = (m-1)*10 + repelem(nucleus_struct_protein(i).setID,N);
        iter = iter + N;
    end
    % record
    offset_vec = [offset_vec offset_temp];
    fluo_vec = [fluo_vec fluo_temp];
    time_vec = [time_vec time_temp];
    set_vec = [set_vec set_temp];
%     mcp_vec = [mcp_vec mcp_temp];
    dorsal_vec = [dorsal_vec dorsal_temp];
end
set_filter = ~ismember(set_vec,[4,13]);
dorsal_vec = dorsal_vec(set_filter);
set_vec = set_vec(set_filter);
offset_vec = offset_vec(set_filter);
time_vec = time_vec(set_filter);
fluo_vec = fluo_vec(set_filter);

%%% Simple cross-correlations
cmap = brewermap([],'Set3');
close all
% MCP sample vs offset
% nan_ft = ~isnan(mcp_vec); 
% p_mcp_off = polyfit(offset_vec(nan_ft),mcp_vec(nan_ft),1);
% yfit_mcp_off = p_mcp_off(1)*offset_vec+p_mcp_off(2);
% hold on;
% 
% mcp_off_fig = figure;
% hold on
% scatter(offset_vec,mcp_vec_plot)
% plot(offset_vec,yfit_mcp_off,'-','Color','black');
% xlabel('MCP (random site)')
% ylabel('MCP (inferred offset)')
% grid on
% box on
% set(gca,'Fontsize',14)
% ylim([0 450])
% saveas(mcp_off_fig, [FigPath 'offset_vs_mcp.png'])

%%% offset vs. Dorsal
nan_ft = ~isnan(dorsal_vec); 
p_dl_off = polyfit(dorsal_vec(nan_ft),offset_vec(nan_ft),1);
yfit_dl_off = p_dl_off(1)*dorsal_vec+p_dl_off(2);
hold on;

dl_off_fig = figure;
colormap(cmap)
hold on
scatter(dorsal_vec,offset_vec,10,set_vec)
plot(dorsal_vec,yfit_dl_off,'');
xlabel('Dl level')
ylabel('MCP (inferred offset)')
grid on
box on
set(gca,'Fontsize',14)
saveas(dl_off_fig,[FigPath 'offset_vs_dorsal.png'])

% mcp vs. Dorsal
nan_ft = ~isnan(mcp_vec); 
p_dl_mcp = polyfit(dorsal_vec(nan_ft),mcp_vec(nan_ft),1);
yfit_dl_mcp = p_dl_mcp(1)*dorsal_vec+p_dl_mcp(2);
hold on;

% dl_mcp_fig = figure;
% colormap(cmap)
% hold on
% scatter(dorsal_vec,mcp_vec_plot,10,set_vec)
% plot(dorsal_vec,yfit_dl_mcp,'r-.');
% xlabel('Dl level')
% ylabel('MCP (random site)')
% ylim([0 450])
% grid on
% box on
% set(gca,'Fontsize',14)
% saveas(dl_mcp_fig,[FigPath 'mcp_vs_dorsal.png'])


time_mcp_fig = figure;
colormap(cmap)
hold on
scatter(time_vec/60,offset_vec,10,set_vec)
grid on
box on
xlabel('time (minutes)')
ylabel('MCP (inferred offset)')
set(gca,'Fontsize',14)
saveas(time_mcp_fig,[FigPath 'offset_vs_time.png'])


time_dl_fig = figure;
colormap(cmap)
hold on
scatter(time_vec/60,dorsal_vec,10,set_vec)
grid on
box on
xlabel('time (minutes)')
ylabel('dorsal concentration')
set(gca,'Fontsize',14)
saveas(time_dl_fig,[FigPath 'dorsal_vs_time.png'])

%%% Estimate degree to which time can account for correlation
% get average offset as a function of dorsal
dorsal_index = linspace(prctile(dorsal_vec,1),prctile(dorsal_vec,99),25);
dl_window = 1;
nBoots = 100;
offset_array = NaN(nBoots,numel(dorsal_index));
time_array = NaN(nBoots,numel(dorsal_index));
for d = 1:numel(dorsal_index)
    dl = dorsal_index(max(1,d-1));
    dh = dorsal_index(min(numel(dorsal_index),d+1));
    index_vec = find(dorsal_vec<dh & dorsal_vec >= dl);
    for n = 1:nBoots
        s_ids = randsample(index_vec,numel(index_vec),true);
        offset_array(n,d) = nanmean(offset_vec(s_ids));
        time_array(n,d) = nanmean(time_vec(s_ids));
    end
end

offset_mean_vec = nanmean(offset_array);
offset_se_vec = nanstd(offset_array);

%%% permorm linear regressions

% generate table structure for regression
regression_tbl = array2table([time_vec' dorsal_vec' ...
    offset_vec'],'VariableNames',{'time','dorsal','offset'});

regression_tbl.setID = categorical(set_vec)';

lm_time = fitlm(regression_tbl,'offset~time');
lm_set = fitlm(regression_tbl,'offset~time*setID');
lm_dl = fitlm(regression_tbl,'offset~time + dorsal');

%%% generate model predictions
time_only_pd = NaN(size(dorsal_index));
time_set_pd = NaN(size(dorsal_index));
time_dl_pd = NaN(size(dorsal_index));

for d = 1:numel(dorsal_index)
    dl = dorsal_index(max(1,d-1));
    dh = dorsal_index(min(numel(dorsal_index),d+1));
    
    index_ft = dorsal_vec<dh & dorsal_vec >= dl;
    temp_tbl = array2table([time_vec(index_ft)' dorsal_vec(index_ft)' ...
                offset_vec(index_ft)'],'VariableNames',{'time','dorsal','offset'});
    temp_tbl.setID = categorical(set_vec(index_ft))';
    % generate predictions
    tp = predict(lm_time,temp_tbl);
    dp = predict(lm_dl,temp_tbl);
    sp = predict(lm_set,temp_tbl);
 
    time_only_pd(d) = nanmean(tp);
    time_set_pd(d) = nanmean(sp);
    time_dl_pd(d) = nanmean(dp);
end


%%% Make figure
time_exp_fig = figure;
hold on
colormap(cmap);
% actual trend
% errorbar(dorsal_index, offset_mean_vec,offset_se_vec,'o','Color','black')
e = scatter(dorsal_index, offset_mean_vec,15,'MarkerEdgeColor','black','MarkerfaceColor','black');
% prediction with time only
p1 = plot(dorsal_index,time_only_pd,'Color',cmap(5,:),'LineWidth',1.5);
p2 = plot(dorsal_index,time_dl_pd,'Color',cmap(6,:),'LineWidth',1.5);
p3 = plot(dorsal_index,time_set_pd,'Color',cmap(7,:),'LineWidth',1.5);
legend('raw data','fit (time only)','fit (time + dorsal)','fit (time x embryo effects)','Location','northwest')
xlabel('Dorsal concentration')
ylabel('MCP (inferred offset')
set(gca,'Fontsize',14)
grid on
box on
% xlim([100,377])
saveas(time_exp_fig,[FigPath 'time_contribution.png'])
saveas(time_exp_fig,[FigPath 'time_contribution.pdf'])

%%% Test hypothesis that embryo-sepcfic effects can explain apparent Dl vs. MCP correlation
lm_mcp_set = fitlm(regression_tbl,'offset~time*setID');
lm_dl_set = fitlm(regression_tbl,'dorsal~time*setID');

dl_vars1 = [0 ; lm_dl_set.Coefficients.Estimate(3:7)];
mcp_vars1 = [0 ; lm_mcp_set.Coefficients.Estimate(3:7)];

dl_vars2 = lm_dl_set.Coefficients.Estimate(8:10);
mcp_vars2 = lm_mcp_set.Coefficients.Estimate(8:10);

embryo_effects = figure;
hold on
colormap(cmap);
scatter(dl_vars1,mcp_vars1)
scatter(dl_vars2,mcp_vars2)
xlabel('Dorsal effects')
ylabel('MCP effects')
grid on
box on
set(gca,'Fontsize',14)
saveas(embryo_effects,[FigPath 'embryo_effects.png'])

%% These results suggest that, within a single set, temporal variability alone 
 % should be able to account for apparent correlation...
 
% generate table structure for regression
setID = 1;
set_tbl = regression_tbl(regression_tbl.setID=='1',:);

lm_time = fitlm(set_tbl,'offset~time')
lm_dl_time = fitlm(set_tbl,'offset~time+dorsal')

lm_dl = fitlm(regression_tbl,'offset~time + dorsal');