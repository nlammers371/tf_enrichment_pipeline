% Script conduct locus enrichment analyses
function main04_make_exploratory_figs(project)

close all

% extract protein, gene, fluorophore info
underscores = strfind(project,'_');
protein_name = project(1:underscores(1)-1);
protein_fluor = project(underscores(1)+1:underscores(2)-1);
gene_name = project(underscores(2)+1:underscores(3)-1);
if numel(underscores) == 3
    ind = numel(project);
else
    ind = underscores(4)-1;
end
gene_fluor = project(underscores(3)+1:end);
% parameters for plots
dist_lim = .6; % min distance from edge permitted (um)
n_boots = 100; % number of bootstraps to use to estimate standard error

% Load analysis data
ReadPath = ['../../dat/' project '/'];
load([ReadPath 'nucleus_struct_protein.mat']);

% Set write path
FigPath = ['../../fig/' project '/'];
mkdir(FigPath)

%%%% Start with some simple bulk analyses
PixelSize = nucleus_struct_ctrl(1).PixelSize;

dist_vec = [nucleus_struct_ctrl.edgeDistSpot]*PixelSize;
ctrl_vec = [nucleus_struct_ctrl.qc_flag_vec]; 
review_vec = [nucleus_struct_ctrl.qc_flag_vec]; 
ctrl_vec(review_vec==0) = 0;

pt_null_vec = [nucleus_struct_ctrl.pt_null];
pt_spot_vec = [nucleus_struct_ctrl.pt_spot];

%%%% Filter for samples meeting qc standards
qc_ft = dist_vec >= dist_lim & ctrl_vec == 1;
dist_vec_ctrl = dist_vec(qc_ft);
pt_null_vec_ctrl = pt_null_vec(qc_ft);
pt_spot_vec_ctrl = pt_spot_vec(qc_ft);

% filter for data points with appropriate control

[~,p] = kstest2(pt_spot_vec_ctrl,pt_null_vec_ctrl);
lb = prctile([pt_null_vec_ctrl pt_spot_vec_ctrl],1);
ub = prctile([pt_null_vec_ctrl pt_spot_vec_ctrl],99);
bins = linspace(lb,ub,100);

spot_ct = histc(pt_spot_vec_ctrl,bins);
null_ct = histc(pt_null_vec_ctrl,bins);
spot_ct = spot_ct / sum(spot_ct);
null_ct = null_ct / sum(null_ct);

hist_fig = figure;
cm = jet(128);
hold on
% hist plots
yyaxis left
b1 = bar(bins,null_ct,1,'FaceAlpha',.5,'FaceColor',cm(30,:));
b2 = bar(bins,spot_ct,1,'FaceAlpha',.5,'FaceColor',cm(100,:));
ylabel('share')
ax = gca;
ax.YColor = 'black';
% cumulative plots
yyaxis right
plot(bins,cumsum(null_ct),'-','LineWidth',1.5,'Color',cm(30,:));
plot(bins,cumsum(spot_ct),'-','LineWidth',1.5,'Color',cm(100,:));
ylabel('cumulative share')
ax = gca;
ax.YColor = 'black';

legend('control',['active locus (' gene_name ')'])
xlabel('concentration (au)')
title(strvcat(['Distribution of ' protein_name '-' protein_fluor ' Concentrations'],...
       [ '    min dist: ' num2str(dist_lim) ' \mu m, p=' num2str(p)]))
grid on
saveas(hist_fig,[FigPath 'hist_plots_' gene_name '_' protein_name '_dist_' num2str(dist_lim) '.png']);    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% compare protein snippets %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nan_ft = ~isnan([nucleus_struct_ctrl.pt_spot]);
snip_ctrl_vec = snip_struct_ctrl.ctrl_flags_final&dist_vec(nan_ft)>dist_lim;

pt_null_mean = nanmean(snip_struct_ctrl.pt_snippet_null(:,:,snip_ctrl_vec==1),3);
pt_spot_mean = nanmean(snip_struct_ctrl.pt_snippet_spot(:,:,snip_ctrl_vec==1),3);
pt_diff_mean = nanmean(snip_struct_ctrl.pt_snippet_spot(:,:,snip_ctrl_vec==1)-snip_struct_ctrl.pt_snippet_null(:,:,snip_ctrl_vec==1),3);

ub = prctile([reshape(pt_null_mean,1,[]) reshape(pt_spot_mean,1,[])],99);
lb = prctile([reshape(pt_null_mean,1,[]) reshape(pt_spot_mean,1,[])],1);

pt_snippet_spot_fig = figure;
colormap(jet(128))
imagesc(pt_spot_mean)
% title(['Active Locus (' gene_name ')'])
caxis([lb ub])
h = colorbar;
ylabel('\mum','FontSize',12)
xlabel('\mum','FontSize',12)
ylabel(h,[protein_name '-' protein_fluor ' concentration (au)'],'FontSize',12)
set(gca,'xtick',3:5:size(pt_null_mean,1),'xticklabel',round(((3:5:size(pt_null_mean,1))-round(size(pt_null_mean,1)/2))*PixelSize,2))
set(gca,'ytick',3:5:size(pt_null_mean,1),'yticklabel',round(((3:5:size(pt_null_mean,1))-round(size(pt_null_mean,1)/2))*PixelSize,2))
saveas(pt_snippet_spot_fig,[FigPath 'mean_pt_snippet_spot.png']);    

pt_snippet_null_fig = figure;
colormap(jet(128))
imagesc(pt_null_mean)
% title('Control')
caxis([lb ub])
h = colorbar;
ylabel('\mum','FontSize',12)
xlabel('\mum','FontSize',12)
ylabel(h,[protein_name '-' protein_fluor ' concentration (au)'],'FontSize',12)
set(gca,'xtick',3:5:size(pt_null_mean,1),'xticklabel',round(((3:5:size(pt_null_mean,1))-round(size(pt_null_mean,1)/2))*PixelSize,2))
set(gca,'ytick',3:5:size(pt_null_mean,1),'yticklabel',round(((3:5:size(pt_null_mean,1))-round(size(pt_null_mean,1)/2))*PixelSize,2))
saveas(pt_snippet_null_fig,[FigPath 'mean_pt_snippet_null.png']);    

% Make diff image
pt_rel_snip = pt_spot_mean ./ pt_null_mean;
diff_snippet_fig = figure;
colormap(jet(128))
imagesc(pt_rel_snip)
% title(['Difference Between Active and Control Snips (' gene_name ')'])
h = colorbar;
ylabel('\mum','FontSize',12)
xlabel('\mum','FontSize',12)
ylabel(h,'fold enrichment (au)','FontSize',12)
set(gca,'xtick',3:5:size(pt_null_mean,1),'xticklabel',round(((3:5:size(pt_null_mean,1))-round(size(pt_null_mean,1)/2))*PixelSize,2))
set(gca,'ytick',3:5:size(pt_null_mean,1),'yticklabel',round(((3:5:size(pt_null_mean,1))-round(size(pt_null_mean,1)/2))*PixelSize,2))
saveas(diff_snippet_fig,[FigPath 'mean_pt_snippet_diff.png']);    

% compare mcp snippets
null_mean_fluo = nanmean(snip_struct_ctrl.fluo_snippet_null(:,:,snip_ctrl_vec==1),3);
spot_mean_fluo = nanmean(snip_struct_ctrl.fluo_snippet_spot(:,:,snip_ctrl_vec==1),3);
ub = prctile([reshape(null_mean_fluo,1,[]) reshape(spot_mean_fluo,1,[])],99);
lb = prctile([reshape(null_mean_fluo,1,[]) reshape(spot_mean_fluo,1,[])],1);

fluo_snippet_spot_fig = figure;
colormap(jet(128))
imagesc(spot_mean_fluo)
% title(['Active Locus (' gene_name ')'])
caxis([lb ub])
h = colorbar;
ylabel('\mum','FontSize',12)
xlabel('\mum','FontSize',12)
ylabel(h,[gene_name ' expression (au)'],'FontSize',12)
set(gca,'xtick',3:5:size(pt_null_mean,1),'xticklabel',round(((3:5:size(pt_null_mean,1))-round(size(pt_null_mean,1)/2))*PixelSize,2))
set(gca,'ytick',3:5:size(pt_null_mean,1),'yticklabel',round(((3:5:size(pt_null_mean,1))-round(size(pt_null_mean,1)/2))*PixelSize,2))
saveas(fluo_snippet_spot_fig,[FigPath 'mean_fluo_snippet_spot.png']);    

fluo_snippet_null_fig = figure;
colormap(jet(128))
imagesc(null_mean_fluo)
% title('Control')
caxis([lb ub])
h = colorbar;
ylabel('\mum','FontSize',12)
xlabel('\mum','FontSize',12)
ylabel(h,[gene_name ' expression (au)'],'FontSize',12)
set(gca,'xtick',3:5:size(pt_null_mean,1),'xticklabel',round(((3:5:size(pt_null_mean,1))-round(size(pt_null_mean,1)/2))*PixelSize,2))
set(gca,'ytick',3:5:size(pt_null_mean,1),'yticklabel',round(((3:5:size(pt_null_mean,1))-round(size(pt_null_mean,1)/2))*PixelSize,2))
saveas(fluo_snippet_null_fig,[FigPath 'mean_fluo_snippet_null.png']);    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Investigate Edge Artifact %%%%%%%%%%%%%%%%%%%%%%%%
% dist_sigma = .1;
% dist_vec_bah = dist_vec(ctrl_vec==1);
% delta_vec_ctrl =  pt_spot_vec_ctrl - pt_null_vec_ctrl;
% dist_index = 0:.1:2;
% 
% dist_mat = NaN(numel(dist_index),n_boots);
% dist_dv_nn = NaN(numel(dist_index),n_boots);
% dist_dv_sp = NaN(numel(dist_index),n_boots);
% 
% for n = 1:n_boots
%     s_ids = randsample(1:numel(delta_vec_ctrl),numel(delta_vec_ctrl),true);
%     dv_samp = delta_vec_ctrl(s_ids);    
%     nn_samp = pt_null_vec_ctrl(s_ids);
%     sp_samp = pt_spot_vec_ctrl(s_ids);
%     dist_samp = dist_vec_blah(s_ids);        
%     for t = 1:numel(dist_index)
%         d_weights = exp(-.5*((dist_samp-dist_index(t))/dist_sigma).^2);
%         dist_mat(t,n) = nansum(dv_samp.*d_weights) / nansum(d_weights);
%         dist_dv_nn(t,n) = (nansum(nn_samp.*d_weights) / nansum(d_weights));
%         dist_dv_sp(t,n) = (nansum(sp_samp.*d_weights) / nansum(d_weights));
%     end
% end
% % calcualte mean and standard error
% dist_mean = nanmean(dist_mat,2);
% dist_ste = nanstd(dist_mat,[],2);
% 
% dist_dv_nn_mean = nanmean(dist_dv_nn,2);
% dist_dv_nn_ste = nanstd(dist_dv_nn,[],2);
% 
% dist_dv_sp_mean = nanmean(dist_dv_sp,2);
% dist_dv_sp_ste = nanstd(dist_dv_sp,[],2);
% 
% % make dist-dependent fold enrichment figure
% fold_fig = figure;
% e = errorbar(dist_index,100* dist_mean ./ dist_dv_nn_mean,100* dist_ste ./ dist_dv_nn_mean);
% e.CapSize = 0;
% xlabel('distance from edge (\mu m)')
% ylabel('apparent % enrichment')
% title(['Pointwise Percent Enrixhment vs. Distance from Nucleus Boundary'])
% grid on
% y_max = 100 * nanmax(dist_mean ./ dist_dv_nn_mean);
% y_min = 100 * nanmin(dist_mean ./ dist_dv_nn_mean);
% ylim([min([0 , y_min-.1*y_min*sign(y_min)]) 1.1*y_max])
% saveas(fold_fig, [FigPath 'edge_artifact_plot.png'])

% plot protein profiles as function edge distance
% fold_fig = figure;
% hold on
% e = errorbar(dist_index, dist_dv_nn_mean, dist_dv_nn_ste);
% e.CapSize = 0;
% e = errorbar(dist_index, dist_dv_sp_mean, dist_dv_sp_ste);
% e.CapSize = 0;
% legend('control', [gene_name ' locus'], 'Location', 'best')
% xlabel('distance from edge (\mu m)')
% ylabel('apparent % enrichment')
% title('Pointwise Percent Enrixhment vs. Distance from Nucleus Boundary')
% grid on
% y_max = nanmax(dist_dv_sp_mean);
% ylim([0 1.1*y_max])
% saveas(fold_fig, [FigPath 'edge_dist_pt_profiles.png'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Probe mean difference between control and spot samples over time %%%%%%
delta_vec_ctrl = pt_spot_vec_ctrl - pt_null_vec_ctrl;

%%% make histogram of differences
delta_hist = figure;
histogram(delta_vec_ctrl)
xlabel('pointwise difference (au)')
ylabel('counts')
title('Distribution of Pointwise Differences Between Locus and Control')
grid on
saveas(delta_hist,[FigPath 'diff_hist.png'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% plot relative enrichment at locus as a function of time %%%%%%%%%%%
time_vec = [nucleus_struct_ctrl.time];
time_vec_dist = time_vec(ctrl_vec==1&dist_vec>dist_lim);

% take bootstrap samples
time_mean = 0:20:round(max(time_vec_dist));
time_delta_mat = NaN(numel(time_mean),n_boots);
sp_dv_time_mat = NaN(numel(time_mean),n_boots);
nn_dv_time_mat = NaN(numel(time_mean),n_boots);
t_window = 120; % seconds
for n = 1:n_boots
    s_ids = randsample(1:numel(delta_vec_ctrl),numel(delta_vec_ctrl),true);
    dv_samp = delta_vec_ctrl(s_ids);    
    nn_samp = pt_null_vec_ctrl(s_ids);
    sp_samp = pt_spot_vec_ctrl(s_ids);
    time_samp = time_vec_dist(s_ids);
    for t = 1:numel(time_mean)
        if time_mean(t) < min(time_vec_dist) || time_mean(t) > max(time_vec_dist)
            continue
        end
        t_weights = exp(-.5*((time_samp-time_mean(t))/t_window).^2);
        time_delta_mat(t,n) = nansum(dv_samp.*t_weights)/nansum(t_weights);        
        nn_dv_time_mat(t,n) =  nansum(nn_samp.*t_weights)/nansum(t_weights);
        sp_dv_time_mat(t,n) = nansum(sp_samp.*t_weights)/nansum(t_weights);
    end
end
% calculate average and standard error
time_delta_mean = nanmean(time_delta_mat,2);
time_delta_ste = nanstd(time_delta_mat,[],2);
nn_time_mean = nanmean(nn_dv_time_mat,2);
nn_time_ste = nanstd(nn_dv_time_mat,[],2);
sp_time_mean = nanmean(sp_dv_time_mat,2);
sp_time_ste = nanstd(sp_dv_time_mat,[],2);


fold_fig = figure;
e = errorbar(time_mean/60,100 * time_delta_mean ./ nn_time_mean,100 * time_delta_ste ./ nn_time_mean);
e.CapSize = 0;
xlabel('minutes')
ylabel('% enrichment at locus')
title(['Percent Enrichment of ' protein_name ' at ' gene_name ' Loci'])
grid on
y_max = nanmax(100 * time_delta_mean ./ nn_time_mean);
y_min = nanmin(100 * time_delta_mean ./ nn_time_mean);
ylim([y_min-.1*y_min*sign(y_min) 1.1*y_max])
saveas(fold_fig, [FigPath 'time_fold_enrichment.png'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Make Radial Profile Figure %%%%%%%%%%%%%%%%%%%%%%%%%%
% generate position reference matrices

r_sigma = .1;
[yref, xref] = meshgrid(1:size(pt_diff_mean,1),1:size(pt_diff_mean,2));
r_grid = sqrt((xref-round(size(pt_diff_mean,1)/2)).^2+(yref-round(size(pt_diff_mean,1)/2)).^2)*PixelSize;
r_index = round(0:.1:2,1);
r_vec = NaN(size(r_index));
% extract_snips
pt_snips_spot = snip_struct_ctrl.pt_snippet_spot(:,:,snip_ctrl_vec==1);
pt_snips_null = snip_struct_ctrl.pt_snippet_null(:,:,snip_ctrl_vec==1);
index_vec = 1:size(pt_snips_null,3);
% initialize profile arrays
r_spot_mat = NaN(n_boots,numel(r_index));
r_null_mat = NaN(n_boots,numel(r_index));

for n = 1:n_boots
    s_ids = randsample(index_vec,numel(index_vec),true);
    pt_spot = nanmean(pt_snips_spot(:,:,s_ids),3);
    pt_null = nanmean(pt_snips_null(:,:,s_ids),3);
    for r = 1:numel(r_index)
        r_weights = exp(-.5*((r_grid-r_index(r))/r_sigma).^2);
        r_spot_mat(n,r) = nansum(pt_spot(:).*r_weights(:)) ./ nansum(r_weights(:));
        r_null_mat(n,r) = nansum(pt_null(:).*r_weights(:)) ./ nansum(r_weights(:));
    end
end

r_spot_mean = nanmean(r_spot_mat);
r_spot_ste = nanstd(r_spot_mat);

r_null_mean = nanmean(r_null_mat);
r_null_ste = nanstd(r_null_mat);

% make figure
cm = jet(128);
r_fig = figure;
hold on
e = errorbar(r_index,r_null_mean / r_null_mean(1),r_null_ste / r_null_mean(1),'Color','black','LineWidth',1.75);
e.CapSize = 0;
e = errorbar(r_index,r_spot_mean / r_null_mean(1),r_spot_ste / r_null_mean(1),'Color',cm(35,:),'LineWidth',1.75);
e.CapSize = 0;
xlabel('radius (\mu m)')
ylabel('relative enrichment')
legend('control','active locus')
% grid on

saveas(r_fig, [FigPath 'radial_enrichment.png'])
