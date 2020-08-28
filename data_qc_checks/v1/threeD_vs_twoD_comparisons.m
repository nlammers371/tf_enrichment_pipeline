clear
close all
addpath('utilities')

% set paths
DropboxFolder =  'E:\Nick\LivemRNA\Dropbox (Personal)\';
project = '2xDl-Ven_snaBAC-mCh';
[~, DataPath, FigureRoot] =   header_function(DropboxFolder, project);
FigPath = [FigureRoot 'dorsal_sna_enrichment_titration/'];
mkdir(FigPath)
% load data
load([DataPath 'nucleus_struct_protein.mat'])
load([DataPath 'nucleus_struct.mat'])
%% basic plot and data qc params 
DistLim = 0.8; % min distance from edge permitted (um)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pull useful vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% distance from nucleus edge
dist_vec = [nucleus_struct_protein.spot_edge_dist_vec];
dist_ft_vec = dist_vec >= DistLim;

close all
xPos2D = [nucleus_struct.xPosParticle];
% fluo_vec = [nucleus_struct_protein.fluo];
xPos3D = [nucleus_struct.xPosParticle3D];

figure;
scatter(xPos2D,xPos3D)
grid on
xlabel('x pos (2D fit)')
ylabel('x pos (3D fit)')

yPos2D = [nucleus_struct.yPosParticle];
yPos3D = [nucleus_struct.yPosParticle3D];
figure;
scatter(yPos2D,yPos3D)
grid on
xlabel('y pos (2D fit)')
ylabel('y pos (3D fit)')

zPos2D = [nucleus_struct.zPosParticle];
zPos3D = [nucleus_struct.zPosParticle3D];
figure;
scatter(zPos2D,zPos3D)
grid on
xlabel('z pos (2D fit)')
ylabel('z pos (3D fit)')

fluo = [nucleus_struct.fluo];
fluo3D = [nucleus_struct.fluo3D];
figure;
scatter(fluo,fluo3D)
grid on
xlabel('fluo (2D fit)')
ylabel('fluo (3D fit)')

%% Examine protein fields
close all
% absolute spot enrichment
delta_protein_vec = ([nucleus_struct_protein.spot_protein_vec] - [nucleus_struct_protein.edge_null_protein_vec]);
delta_protein_dist_vec = delta_protein_vec(dist_ft_vec);
delta_protein3D_vec = ([nucleus_struct_protein.spot_protein_vec_3d] - [nucleus_struct_protein.edge_null_protein_vec_3d]);
delta_protein3D_dist_vec = delta_protein3D_vec(dist_ft_vec);

lb = prctile([delta_protein_dist_vec delta_protein3D_dist_vec],.1);
ub = prctile([delta_protein_dist_vec delta_protein3D_dist_vec],99.9);
pt_bins = linspace(lb,ub);

hist_fig = figure;
hold on
histogram(delta_protein_dist_vec,pt_bins,'Normalization','probability')
histogram(delta_protein3D_dist_vec,pt_bins,'Normalization','probability')
legend('2D','3D')
xlabel('Dorsal enrichment')
ylabel('share')
grid on
box on

