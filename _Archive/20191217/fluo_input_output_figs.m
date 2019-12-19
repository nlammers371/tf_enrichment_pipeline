clear 
% close all
addpath('utilities')
% set ID variables
targetProject = 'Dl-Ven_snaBAC-mCh';
controlProject = 'Dl-Ven_hbP2P-mCh';
DropboxFolder = 'E:\Nick\LivemRNA\Dropbox (Personal)\';
[~, DataPathTarget, FigureRoot] =   header_function(DropboxFolder, targetProject); 
[~, DataPathControl, ~] =   header_function(DropboxFolder, controlProject); 

FigPath = [FigureRoot '\' targetProject '\input_output01\'];
mkdir(FigPath)
% load data
load([DataPathTarget 'fluo_input_output.mat'])
target_results_struct = results_struct;
load([DataPathControl 'fluo_input_output.mat'])
control_results_struct = results_struct;
clear results_struct

% plot simple figure with just computational control
time_axis = target_results_struct(1).ref_vec * 20/60;
% calculate upper and lower bound vectors
target_mean = target_results_struct(1).spot_protein_mean;
br_spot_ub = target_mean + target_results_struct(1).spot_protein_ste;
br_spot_lb = target_mean - target_results_struct(1).spot_protein_ste;

virt_mean = target_results_struct(1).serial_protein_mean;
br_virt_ub = virt_mean + target_results_struct(1).serial_protein_ste;
br_virt_lb = virt_mean - target_results_struct(1).serial_protein_ste;

sna_mean = target_results_struct(1).response_mean;

% make figure
burst_trend_fig = figure;
cmap1 = brewermap([],'Set2');

% snail activity
yyaxis right
p1 = plot(time_axis,sna_mean,'--','Color','black','LineWidth',2);
ylabel('snail transcription (au)')
% set(gca,'ytick',.2:.1:1.2)
% ylim([.2 1.2])
ax = gca;
ax.YColor = 'black';

% Dorsal activity
yyaxis left
hold on
fill([time_axis fliplr(time_axis)],[br_spot_ub fliplr(br_spot_lb)],cmap1(2,:),'FaceAlpha',.5,'EdgeAlpha',0)
p2 = plot(time_axis,target_mean,'-','Color',cmap1(2,:),'LineWidth',2);

fill([time_axis fliplr(time_axis)],[br_virt_ub fliplr(br_virt_lb)],cmap1(3,:),'FaceAlpha',.5,'EdgeAlpha',0)
p3 = plot(time_axis,virt_mean,'-','Color',cmap1(3,:),'LineWidth',2);

ylabel('relative Dl concentration (au)')
% set(gca,'ytick',-20:4:20)
ax = gca;
ax.YColor = 'black';

grid on
xlabel('offset (minutes)')
legend([p1 p2 p3],'{\it snail} MS2','Dl at {\it snail} locus','Dl at control locus', 'Location','northeast')
set(gca,'Fontsize',12,'xtick',-4:2:4)
chH = get(gca,'Children');
set(gca,'Children',flipud(chH));
% save
saveas(burst_trend_fig,[FigPath 'fluo_locus_trend.tif'])
saveas(burst_trend_fig,[FigPath 'fluo_locus_trend.pdf'])


% plot with ALL the controls

% make figure
burst_trend_all_fig = figure;
cmap1 = brewermap([],'Set2');

% snail activity
yyaxis right
p1 = plot(time_axis,sna_mean,'--','Color','black','LineWidth',2);
ylabel('snail transcription (au)')
ax = gca;
ax.YColor = 'black';

% Dorsal activity
yyaxis left
hold on
fill([time_axis fliplr(time_axis)],[br_spot_ub fliplr(br_spot_lb)],cmap1(2,:),'FaceAlpha',.5,'EdgeAlpha',0)
p2 = plot(time_axis,target_mean,'-','Color',cmap1(2,:),'LineWidth',2);
% virtual control
p3 = plot(time_axis,virt_mean,'-','Color',cmap1(3,:),'LineWidth',2);
% swap control
p4 = plot(time_axis,target_results_struct(1).swap_protein_mean,'-','Color',cmap1(5,:),'LineWidth',2);
%  hbP2P control
p5 = plot(time_axis,control_results_struct(1).spot_protein_mean,'-','Color',cmap1(6,:),'LineWidth',2);

ylabel('relative Dl concentration (au)')
% set(gca,'ytick',-20:4:20)
ax = gca;
ax.YColor = 'black';

grid on
xlabel('offset (minutes)')
legend([p1 p2 p3 p4 p5],'{\it snail} MS2','Dl at {\it snail} locus','Dl at control locus',...
    'Dl at nearest neighbor locus','Dl at {\it hbP2P} locus','Location','northeast')
set(gca,'Fontsize',12,'xtick',-4:2:4)
chH = get(gca,'Children');
set(gca,'Children',flipud(chH));
% save
saveas(burst_trend_all_fig,[FigPath 'fluo_locus_trend_all_ctrl.tif'])
saveas(burst_trend_all_fig,[FigPath 'fluo_locus_trend_all_ctrl.pdf'])

% with bio control only
% make figure
burst_trend_bio_fig = figure;
cmap1 = brewermap([],'Set2');

bio_mean = control_results_struct(1).spot_protein_mean;
br_bio_ub = bio_mean + target_results_struct(1).spot_protein_ste;
br_bio_lb = bio_mean - target_results_struct(1).spot_protein_ste;
% snail activity
yyaxis right
p1 = plot(time_axis,sna_mean,'--','Color','black','LineWidth',2);
ylabel('snail transcription (au)')
ax = gca;
ax.YColor = 'black';

% Dorsal activity
yyaxis left
hold on
fill([time_axis fliplr(time_axis)],[br_spot_ub fliplr(br_spot_lb)],cmap1(2,:),'FaceAlpha',.5,'EdgeAlpha',0)
p2 = plot(time_axis,target_mean,'-','Color',cmap1(2,:),'LineWidth',2);
%  hbP2P control
fill([time_axis fliplr(time_axis)],[br_bio_ub fliplr(br_bio_lb)],cmap1(6,:),'FaceAlpha',.5,'EdgeAlpha',0)
p3 = plot(time_axis,bio_mean,'-','Color',cmap1(6,:),'LineWidth',2);


ylabel('relative Dl concentration (au)')
% set(gca,'ytick',-20:4:20)
ax = gca;
ax.YColor = 'black';

grid on
xlabel('offset (minutes)')
legend([p1 p2 p3],'{\it snail} MS2','Dl at {\it snail} locus','Dl at {\it hbP2P} locus','Location','northeast')
set(gca,'Fontsize',12,'xtick',-4:2:4)
chH = get(gca,'Children');
set(gca,'Children',flipud(chH));
% save
saveas(burst_trend_bio_fig,[FigPath 'fluo_locus_trend_bio_ctrl.tif'])
saveas(burst_trend_bio_fig,[FigPath 'fluo_locus_trend_bio_ctrl.pdf'])
