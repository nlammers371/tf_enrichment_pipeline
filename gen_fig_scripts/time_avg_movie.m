% Script to make movie of time-averaged local enrichment
clear
close all
% project = 'Bcd_GFP_hb_mCherry_Zoom2x';
project = 'Bcd_GFP_hb_mCherry_Zoom3x_LowPower';
% extract protein, gene, fluorophore info
underscores = strfind(project,'_');
protein_name = project(1:underscores(1)-1);
protein_fluor = project(underscores(1)+1:underscores(2)-1);
gene_name = project(underscores(2)+1:underscores(3)-1);
gene_fluor = project(underscores(3)+1:underscores(4)-1);
% load snips
ReadPath = ['../../../dat/' project '/'];
load([ReadPath 'snip_struct_ctrl.mat'])
load([ReadPath 'nucleus_struct_ctrl.mat']);
PixelSize = nucleus_struct_ctrl(1).PixelSize;
ctrl_vec = snip_struct_ctrl.ctrl_flags_final;
% make save path
WritePath = ['../../../fig/' project '/time_avg_frames/'];
mkdir(WritePath)
% extract snip stacks
pt_spot_snips = snip_struct_ctrl.pt_snippet_spot(:,:,ctrl_vec==1);
% set movie parameters
frame_inc = 25;
n_frames = ceil(size(pt_spot_snips,3)/frame_inc);
index_vec = 1:size(pt_spot_snips,3);
% shuffle
pt_spot_snips = pt_spot_snips(:,:,randsample(index_vec,numel(index_vec),false));
full_mean = nanmean(pt_spot_snips,3);
% generate position reference matrices
[yref, xref] = meshgrid(1:size(full_mean,1),1:size(full_mean,2));
r_grid = sqrt((xref-round(size(full_mean,1)/2)).^2+(yref-round(size(full_mean,1)/2)).^2)*PixelSize;
r_index = 0:.1:2;
% loop through frames
parfor n = 1:n_frames
    mean_frame = nanmean(pt_spot_snips(:,:,1:min([frame_inc*n,numel(index_vec)])),3);
    temp_fig = figure('Visible','off');
    colormap(jet(128));
%     imagesc(imgaussfilt(mean_frame,1))
    imagesc(mean_frame)
%     caxis([lb ub])
%     h = colorbar;
% 	zlim([1600 2400])
    set(gca,'xtick',3:5:size(full_mean,1),'xticklabel',round(((3:5:size(full_mean,1))-round(size(full_mean,1)/2))*PixelSize,2))
    set(gca,'ytick',3:5:size(full_mean,1),'yticklabel',round(((3:5:size(full_mean,1))-round(size(full_mean,1)/2))*PixelSize,2))
    yl = get(gca,'ylim');
    xl = get(gca,'xlim');
%     title([num2str(n*frame_inc) ' frames'])
    text(xl(1)+.25,yl(2) - 1,[num2str(n*frame_inc) ' frames'],...
        'Color','black','BackgroundColor','white')    
    saveas(temp_fig,[WritePath 'temp_avg_' sprintf('%03d',n) '.tif'])            
    close all
    % generate radial profile plots
    r_vec = NaN(size(r_index));
    for i = 1:numel(r_index)
        r_vec(i) = nanmean(mean_frame(round(r_grid,1)==round(r_index(i),1)));
    end
    r_fig = figure('Visible','off');
    plot(r_index,r_vec,'LineWidth',1.5);
    grid on
    xlabel('\mum')
    ylabel('au')
    ylim([500 650])
    yl = get(gca,'ylim');
    xl = get(gca,'xlim');
    text(xl(1)+.02,yl(1)+6,[num2str(n*frame_inc) ' frames'],...
        'Color','black','BackgroundColor','white')       
%     error('asfa')
    saveas(r_fig,[WritePath 'r_avg_' sprintf('%03d',n) '.tif'])            
    
    %%% Make single-frame versions
    single_frame = pt_spot_snips(:,:,min([numel(index_vec),n*frame_inc]));
    temp_fig = figure('Visible','off');
    colormap(jet(128));
%     imagesc(imgaussfilt(mean_frame,1))
    imagesc(single_frame)
%     caxis([lb ub])
%     h = colorbar;
% 	zlim([1600 2400])
    set(gca,'xtick',3:5:size(full_mean,1),'xticklabel',round(((3:5:size(full_mean,1))-round(size(full_mean,1)/2))*PixelSize,2))
    set(gca,'ytick',3:5:size(full_mean,1),'yticklabel',round(((3:5:size(full_mean,1))-round(size(full_mean,1)/2))*PixelSize,2))
%     yl = get(gca,'ylim');
%     xl = get(gca,'xlim');
%     title([num2str(n*frame_inc) ' frames'])
%     text(xl(1)+.25,yl(2) - 1,[num2str(n*frame_inc) ' frames'],...
%         'Color','black','BackgroundColor','white')    
    saveas(temp_fig,[WritePath 'single_frame_' sprintf('%03d',n) '.tif'])            
    close all
    % generate radial profile plots
    r_vec = NaN(size(r_index));
    for i = 1:numel(r_index)
        r_vec(i) = nanmean(single_frame(round(r_grid,1)==round(r_index(i),1)));
    end
    r_fig = figure('Visible','off');
    plot(r_index,r_vec,'LineWidth',1.5);
    grid on
    xlabel('\mum')
    ylabel('au')
    ylim([200 800])
%     yl = get(gca,'ylim');
%     xl = get(gca,'xlim');
%     text(xl(1)+.02,yl(1)+6,[num2str(n*frame_inc) ' frames'],...
%         'Color','black','BackgroundColor','white')       
%     error('asfa')
    saveas(r_fig,[WritePath 'r_single_frame_' sprintf('%03d',n) '.tif'])            
end