clear
close all
addpath('utilities')
dropboxFolder = 'E:\Nick\LivemRNA\Dropbox (Personal)\';
dataPath = [dropboxFolder 'ProcessedEnrichmentData\'];
figPath = [dropboxFolder 'LocalEnrichmentFigures\fluo_input_output\'];
mkdir(figPath)
% add path to utilities
addpath('../utilities')
% load data set
load([dataPath 'fluo_in_out.mat'],'fluo_io_struct')
%%%
VoxelSize = fluo_io_struct.voxel_size;
PixelSize = sqrt(VoxelSize/.5);

% make fluo heatmap plots

% extract snips
fluo_q1_snip = fluo_io_struct.fluo_q1_snip;
fluo_q4_snip = fluo_io_struct.fluo_q4_snip;
snip_size = size(fluo_q1_snip,1);

% determine bounds
fluo_lb = prctile([fluo_q1_snip(:)'  fluo_q4_snip(:)'],1);
fluo_ub = prctile([fluo_q1_snip(:)'  fluo_q4_snip(:)'],99);

% define tick strings
xtick_string = "set(gca,'xtick',1:5:snip_size,'xticklabel',round(((1:5:snip_size)-round(snip_size/2))*PixelSize,2))";
ytick_string = "set(gca,'ytick',1:5:snip_size,'yticklabel',round(((1:5:snip_size)-round(snip_size/2))*PixelSize,2))";

% Make titles
spot_mcp_title_dim = 'snail MCP-mCherry (bottom quintile)';
spot_mcp_title_bright = 'snail MCP-mCherry (top quintile)';
spot_mcp_Clabel = 'snail intensity (AU)';

% define heatmap
magenta = [162 60 150]/256/.65;
gray = magenta/3;
cm_magenta = interp1([0,1],vertcat(gray,magenta),linspace(0,1,128));
% dim spot
makeHeatmapPlots(fluo_q1_snip, 1, spot_mcp_title_dim, spot_mcp_Clabel,cm_magenta,PixelSize,fluo_lb,fluo_ub)
makeHeatmapPlots(fluo_q4_snip, 1, spot_mcp_title_bright, spot_mcp_Clabel,cm_magenta,PixelSize,fluo_lb,fluo_ub)

%%% Now protein
% extract snips
protein_q1_snip = fluo_io_struct.protein_q1_snip / VoxelSize;
protein_q4_snip = fluo_io_struct.protein_q4_snip / VoxelSize;

% determine bounds
protein_lb = 0;%prctile([protein_q1_snip(:)'  protein_q4_snip(:)'],1);
protein_ub = floor(max([protein_q1_snip(:)'  protein_q4_snip(:)'])/10)*10;

% Make titles
spot_protein_title_dim = 'Dorsal-Venus at Active snail locus (bottom quintile)';
spot_protein_title_bright = 'Dorsal-Venus at Active snail locus (top quintile)';
spot_protein_Clabel = 'Venus intensity (AU)';

bad_cm = flipud(brewermap([],'RdYlBu'));
% define heatmap
green_top = [12 118 60]/256/.65;
green_bottom = green_top/3;
% cm_green = interp1([0,1],vertcat(green_bottom,green_top),linspace(0,1,128));
% dim spot
makeHeatmapPlots(protein_q1_snip, 1, spot_protein_title_dim, spot_protein_Clabel,bad_cm,PixelSize,protein_lb,protein_ub)
makeHeatmapPlots(protein_q4_snip, 1, spot_protein_title_bright, spot_protein_Clabel,bad_cm,PixelSize,protein_lb,protein_ub)


%%% Now bar plots
protein_target_mean = fluo_io_struct.protein_target_mean;
protein_target_ste = fluo_io_struct.protein_target_ste;
protein_control_mean = fluo_io_struct.protein_control_mean;
protein_control_ste = fluo_io_struct.protein_control_ste;

% calculate y limits
ymax = 1.1*max(protein_target_mean);
ymin = 1.5*min(protein_control_mean);

target_fig = figure;
hold on
b = bar([1 3 5],protein_target_mean);
b.FaceColor = bad_cm(210,:);
errorbar([1 3 5],protein_target_mean,protein_target_ste,'LineStyle','none','Color','black','LineWidth',1.5,'CapSize',20)
box on
xlabel('quintile')
ylabel('Dorsal enrichment (AU)')
set(gca,'xtick',[1 3 5],'FontSize',14);
ylim([ymin ymax])

control_fig = figure;
hold on
b = bar([1 3 5],protein_control_mean);
b.FaceColor = bad_cm(50,:);
errorbar([1 3 5],protein_control_mean,protein_control_ste,'LineStyle','none','Color','black','LineWidth',1.5,'CapSize',20)
box on
xlabel('quintile')
ylabel('Dorsal enrichment (AU)')
set(gca,'xtick',[1 3 5],'FontSize',14);
ylim([ymin ymax])



