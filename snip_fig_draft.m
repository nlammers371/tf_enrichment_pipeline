% Script conduct locus enrichment analyses
clear
close all
%%% id variables
% pt_string = 'Kr-GFP';
% project = 'kr_eve2_reporter';
protein_name = 'Bcd-GFP';
gene_name = 'eve stripe 2';
% gene_name = 'hb';
% gene_name = 'snail';
% project = 'bcd_sna_zoom3x_pinhole06';
project = 'Bcd_eve2_reporter';
% project = 'bcd_hb_zoom3_LowPower';
% parameters for plots
dist_lim = .6;
% project = 'bcd_hb_pos_control';
%%% Load analysis data
ReadPath = ['../dat/' project '/'];
load([ReadPath 'snip_struct_ctrl.mat'])
load([ReadPath 'nucleus_struct_ctrl.mat']);
PixelSize = nucleus_struct_ctrl(1).PixelSize;
% pix_cutoffs = dist_cutoffs / PixelSize;
% Set write path
FigPath = ['../fig/' project '/snip_figures/'];
mkdir(FigPath)

% indexing vectors
dist_vec = [nucleus_struct_ctrl.edgeDistSpot]*PixelSize; 

%%%%%%%%%%%%%%%%%%% compare protein snippets %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
snip_ctrl_vec = snip_struct_ctrl.ctrl_flags_final&dist_vec(~isnan([nucleus_struct_ctrl.pt_spot]))>dist_lim;
null_mean = nanmean(snip_struct_ctrl.pt_snippet_null(:,:,snip_ctrl_vec==1),3);
spot_mean = nanmean(snip_struct_ctrl.pt_snippet_spot(:,:,snip_ctrl_vec==1),3);
pt_diff_snip = nanmean(snip_struct_ctrl.pt_snippet_spot(:,:,snip_ctrl_vec==1)-snip_struct_ctrl.pt_snippet_null(:,:,snip_ctrl_vec==1),3);

ub = prctile([reshape(null_mean,1,[]) reshape(spot_mean,1,[])],99);
% lb = prctile([reshape(null_mean,1,[]) reshape(spot_mean,1,[])],1);

% interpolate/smooth as needed
if ~strcmp(gene_name,'eve stripe 2')
    null_mean = imgaussfilt(null_mean,1);
    spot_mean = imgaussfilt(spot_mean,1);
    pt_diff_snip = imgaussfilt(pt_diff_snip,1);    
    
else
    null_mean = imresize(null_mean,25/size(null_mean,1));
    spot_mean = imresize(spot_mean,25/size(spot_mean,1));
    pt_diff_snip = imresize(pt_diff_snip,25/size(pt_diff_snip,1));
end
mx = max([null_mean(:)' spot_mean(:)']);

null_mean = null_mean / mx;
spot_mean = spot_mean / mx;
diff_pct = 100*((spot_mean - null_mean) ./ null_mean);


null_snippet_fig = figure;
imagesc(null_mean)
colormap(jet(128))
title(['Control (' gene_name ')'])
caxis([.5 1])
h = colorbar;
ylabel('\mum','FontSize',12)
xlabel('\mum','FontSize',12)
ylabel(h,'concentration (au)','FontSize',12)
set(gca,'xtick',3:5:25,'xticklabel',round(((3:5:25)-13)*PixelSize,2))
set(gca,'ytick',3:5:25,'yticklabel',round(((3:5:25)-13)*PixelSize,2))
saveas(null_snippet_fig,[FigPath 'mean_null_snippet_smooth.png']);    

spot_snippet_fig = figure;
colormap(jet(128))
imagesc(spot_mean);
title(['Active ' gene_name ' Locus'])
caxis([.5 1])
h = colorbar;
ylabel('\mum','FontSize',12)
xlabel('\mum','FontSize',12)
ylabel(h,'concentration (au)','FontSize',12)
set(gca,'xtick',3:5:25,'xticklabel',round(((3:5:25)-13)*PixelSize,2))
set(gca,'ytick',3:5:25,'yticklabel',round(((3:5:25)-13)*PixelSize,2))

saveas(spot_snippet_fig,[FigPath 'mean_spot_snippet_smooth.png']);    

% Make diff image
% ub = prctile([pt_diff_snip(:)],99);
% lb = prctile([pt_diff_snip(:)],1);

diff_snippet_fig = figure;
colormap(jet(128))
imagesc(diff_pct)
caxis([-5 15])
h = colorbar;
title(['Percent Bcd Enrichment at active ' gene_name ' locus'])
ylabel('\mum','FontSize',12)
xlabel('\mum','FontSize',12)
ylabel(h,'% enrichment (au)','FontSize',12)
set(gca,'xtick',3:5:25,'xticklabel',round(((3:5:25)-13)*PixelSize,2))
set(gca,'ytick',3:5:25,'yticklabel',round(((3:5:25)-13)*PixelSize,2))
saveas(diff_snippet_fig,[FigPath 'diff_pt_snippet_smooth.png']);    
% %%
% % compare mcp snippets
% null_mean_fluo = nanmean(snip_struct_ctrl.fluo_snippet_null(:,:,snip_ctrl_vec==1),3);
% spot_mean_fluo = nanmean(snip_struct_ctrl.fluo_snippet_spot(:,:,snip_ctrl_vec==1),3);
% ub = prctile([reshape(null_mean_fluo,1,[]) reshape(spot_mean_fluo,1,[])],99);
% lb = prctile([reshape(null_mean_fluo,1,[]) reshape(spot_mean_fluo,1,[])],1);
% 
% mcp_snippet_fig = figure('Position',[0 0 1024 512]);
% colormap(jet(128))
% hold on
% subplot(1,2,1)
% 
% imagesc(null_mean_fluo)
% title('control')
% caxis([lb ub])
% colorbar
% 
% subplot(1,2,2)
% imagesc(spot_mean_fluo)
% title(['active locus (' gene_name ')'])
% caxis([lb ub])
% colorbar
% saveas(mcp_snippet_fig,[FigPath 'mean_fluo_snippet.png']);  
