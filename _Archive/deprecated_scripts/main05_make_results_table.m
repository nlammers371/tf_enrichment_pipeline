% main05_make_results_table(project)
%
% DESCRIPTION
% Script to compile experimental traces into longform table amd generate QC
% Figs
%
% ARGUMENTS
% project: master ID variable 
%
%
% OUTPUT: analysis_data: longform table structure. Script also generates QC
%         Figures

function main05_make_results_table(project)


%-----------------------------ID Variables--------------------------------%
%%% Generate filenames and writepath
DataPath = ['../../dat/' project '/'];
FigPath = ['../../fig/' project '/qc_plots/'];
mkdir(FigPath)
%%% Load Data
load([DataPath '/nucleus_struct_ctrl.mat'])

%%% Make set Key
set_index = unique([nucleus_struct_ctrl.setID]);
source_index = {};
for i = 1:numel(set_index)
    set_struct = nucleus_struct_ctrl([nucleus_struct_ctrl.setID]==set_index(i));
    source_index = [source_index{:} {set_struct(1).source_path}];
end

% Define names of fields to include in final table structure
var_name_cell = {'NucleusID','ParticleID','ctrl_flag','setID','AP','EdgeDist','CenterDist',...
            'xParticle','yParticle','zParticle','xNC','yNC', 'xNull','yNull',...
            'frame','time','fluo_int','fluo_spot','fluo_null','on','pt_mf','pt_spot','pt_null'};
%%% initialize vectors to store output
ParticleID_vec = [];
NucleusID_vec = [];
on_vec = [];
% position info
xParticle_vec = [nucleus_struct_ctrl.xPosParticle];
yParticle_vec = [nucleus_struct_ctrl.yPosParticle];
zParticle_vec = [nucleus_struct_ctrl.brightestZs];
xNC_vec = [nucleus_struct_ctrl.xPos];
yNC_vec = [nucleus_struct_ctrl.yPos];

xNull_vec = [nucleus_struct_ctrl.xPosNull_edge];
yNull_vec = [nucleus_struct_ctrl.yPosNull_edge];

ctrl_flag_vec = [nucleus_struct_ctrl.ctrl_flags];
EdgeDist_vec = [nucleus_struct_ctrl.edgeDistNull];    

AP_vec = [nucleus_struct_ctrl.ap_vector];

frame_vec = [nucleus_struct_ctrl.frames]; 
time_vec = [nucleus_struct_ctrl.time];    

fluo_int_vec = [nucleus_struct_ctrl.fluo];
fluo_spot_vec = [nucleus_struct_ctrl.fluo_spot];
fluo_null_vec = [nucleus_struct_ctrl.fluo_null];

pt_spot_vec = [nucleus_struct_ctrl.pt_spot];
pt_null_vec = [nucleus_struct_ctrl.pt_null];
pt_mf_vec = [nucleus_struct_ctrl.protein];

for i = 1:numel(nucleus_struct_ctrl)
    NucleusID = nucleus_struct_ctrl(i).ncID;
    ParticleID = nucleus_struct_ctrl(i).ParticleID;    
    fluo = nucleus_struct_ctrl(i).fluo;
    % make on field
    on = NaN(size(fluo));
    on(find(~isnan(fluo),1):end) = 0;
    on(find(~isnan(fluo),1):find(~isnan(fluo),1,'last')) = 1;
    on_vec = [on_vec on];        
    % get mf pt field    
    ParticleID_vec = [ParticleID_vec repelem(ParticleID,numel(fluo))];        
    NucleusID_vec = [NucleusID_vec repelem(NucleusID,numel(fluo))];                    
end
setID_vec = floor(ParticleID_vec);
% vector of particle distances from nucleus center
CenterDist_vec = sqrt((xNC_vec-xParticle_vec).^2+(yNC_vec-yParticle_vec).^2);

% concatenate fields
data_array = double(NucleusID_vec');
for f = 2:numel(var_name_cell)
    var_vec = eval([var_name_cell{f} '_vec']);
    data_array = [data_array double(var_vec')];
end

out_table = array2table(data_array,'VariableNames', var_name_cell);        

%%% Quick QC checks
close all
% fluo fields
fluo_fluo_fig = figure;
hold on
title('Fluo Field Consistency')
scatter(out_table.fluo_int,out_table.fluo_spot,20)
scatter(out_table.fluo_int,out_table.fluo_null,20)
xlabel('pipeline (AU)')
ylabel('raw (AU)')
legend('spot','control')
grid on
saveas(fluo_fluo_fig,[FigPath 'fluo_checks.png'])

% protein fields
pt_pt_fig = figure;
hold on
title('Protein Field Consistency')
scatter(out_table.pt_mf,out_table.pt_spot,20)
scatter(out_table.pt_mf,out_table.pt_null,20)
xlabel('pipeline (AU)')
ylabel('raw (AU)')
legend('spot','control')
grid on
saveas(pt_pt_fig,[FigPath 'protein_checks.png'])
% position variables
x_x_fig = figure;
hold on
title('X Field Consistency')
scatter(out_table.xNC,out_table.xParticle,20)
scatter(out_table.xNC,out_table.xNull,20)
xlabel('nucleus (AU)')
ylabel('locus (AU)')
legend('spot','control')
grid on
saveas(x_x_fig,[FigPath 'x_checks.png'])

y_y_fig = figure;
hold on
title('Y Field Consistency')
scatter(out_table.yNC,out_table.yParticle,20)
scatter(out_table.yNC,out_table.yNull,20)
xlabel('nucleus (AU)')
ylabel('locus (AU)')
legend('spot','control')
grid on
saveas(y_y_fig,[FigPath 'y_checks.png'])
% save
writetable(out_table,[DataPath 'analysis_data.csv'])    