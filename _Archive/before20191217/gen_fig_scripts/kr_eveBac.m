% Script conduct locus enrichment analyses
clear
close all
project = 'Kr_GFP_eveBAC_mCherry';
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
dist_lim = .6;
n_boots = 100; % number of bootstraps to use to estimate standard error
% project = 'bcd_hb_pos_control';
%%% Load analysis data
ReadPath = ['../../../dat/' project '/'];
load([ReadPath 'snip_struct_ctrl.mat'])
load([ReadPath 'nucleus_struct_ctrl.mat']);

% Set write path
FigPath = ['../../fig/' project '/'];
mkdir(FigPath)
%%
PixelSize = nucleus_struct_ctrl(1).PixelSize;
%%%% Start with some simple bulk analyses
dist_vec = [nucleus_struct_ctrl.edgeDistSpot]*PixelSize;
ctrl_vec = [nucleus_struct_ctrl.ctrl_flags]; 
pt_null_vec = [nucleus_struct_ctrl.pt_null];
pt_spot_vec = [nucleus_struct_ctrl.pt_spot];
fluo_vec = [nucleus_struct_ctrl.fluo];
ap_vec = [nucleus_struct_ctrl.ap_vector];
time_vec = [nucleus_struct_ctrl.time];
set_vec= [];
for i = 1:numel(nucleus_struct_ctrl)
    set_vec = [set_vec repelem(nucleus_struct_ctrl(i).setID,numel(nucleus_struct_ctrl(i).fluo))];
end
% filter 
ft = ctrl_vec == 1 & dist_vec > dist_lim & ap_vec > .3 & ap_vec < .4 & time_vec > 25*60;
pt_null_vec = pt_null_vec(ft);
pt_spot_vec = pt_spot_vec(ft);
fluo_vec = fluo_vec(ft);
time_vec = time_vec(ft);
ap_vec = ap_vec(ft);
