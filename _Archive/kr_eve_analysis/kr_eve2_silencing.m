% Script to dissect Kr dynamics leading up to eve2 silencing
clear
close all
%%% id variables
Tres = 20;
pt_string = 'Kruppel';
project = 'kr_eve2_reporter';
print_traces = 1;
%%% Load analysis data
ReadPath = ['../../dat/' project '/'];
load([ReadPath 'snip_struct_ctrl.mat'])
regression_table = readtable([ReadPath project '_analysis_data.csv']);

% Set write path
FigPath = ['../../fig/' project '/'];
mkdir(FigPath)

%%% Make indexing vectors
include_vec = [1 5 6 8 10 12];
t_bounds = [15 35]*60;
min_ap = 41/100;

set_vec = regression_table.setID;
ap_vec = regression_table.AP;
t_vec = regression_table.time;

ft_vec = ap_vec>=min_ap&t_vec>=t_bounds(1)&t_vec<t_bounds(2);
analysis_table = regression_table(ft_vec,:);
%%
%%% Plot Mean Kr Concentration in Locus and Control Snippets as a Function
%%% of Time till turn-off
pt_spot_vec = analysis_table.pt_spot;
pt_null_vec = analysis_table.pt_null;
time_vec = analysis_table.time;
fluo_vec = analysis_table.fluo_int;
ncID_vec = analysis_table.NucleusID;
ncID_index = unique(ncID_vec);

% generate relative time variable
time_vec_rel = NaN(size(time_vec));
for i = 1:numel(ncID_index)
    f_vec = fluo_vec(ncID_vec==ncID_index(i));
    t_vec = time_vec(ncID_vec==ncID_index(i));
    lt = find(~isnan(f_vec),1,'last');
    if lt <= numel(f_vec)-5
        time_vec_rel(ncID_vec==ncID_index(i)) = t_vec - t_vec(lt);
    end
end
%%
% track protein and fljuo levels in 30 min prior to turn off
fluo_vec_rel = NaN(1,31);
pt_spot_vec_rel = NaN(1,31);
pt_null_vec_rel = NaN(1,31);
t_rel_index = fliplr(-30:0);
for i = 1:numel(fluo_vec_rel)
    t_ft = round(time_vec_rel/60)==t_rel_index(i);
    fluo_vec_rel(i) = nanmean(fluo_vec(t_ft));
    pt_spot_vec_rel(i) = nanmean(pt_spot_vec(t_ft));
    pt_null_vec_rel(i) = nanmean(pt_null_vec(t_ft));
end
%%
pt_spot_vec_late = pt_spot_vec(time_vec_rel<=-500);
pt_null_vec_late = pt_null_vec(time_vec_rel<=-500);