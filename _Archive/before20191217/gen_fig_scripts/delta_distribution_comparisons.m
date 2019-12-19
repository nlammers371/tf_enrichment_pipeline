% Script to compare pt distirbutions across different conditions
clear
close all
% Set write path
FigPath = ['../../../fig/various_analyses/'];
mkdir(FigPath)
% specify projects to compare
project_cell = {'Bcd_GFP_hb_mCherry_Zoom2x','Bcd_GFP_snail_mCherry_Zoom2x'};
gene_name_cell = {};
protein_name_cell = {};
gene_fluor_cell = {};
protein_fluor_cell = {};
% set basic analysis parameters
dist_lim = .6;
% initailze structure to store results
master_struct = struct;

for i = 1:numel(project_cell)
    project = project_cell{i};
    % extract protein, gene, fluorophore info
    ReadPath = ['../../../dat/' project '/'];
    underscores = strfind(project,'_');
    protein_name_cell = [protein_name_cell{:} {project(1:underscores(1)-1)}];
    protein_fluor_cell = [protein_fluor_cell{:} {project(underscores(1)+1:underscores(2)-1)}];
    gene_name_cell = [gene_name_cell{:} {project(underscores(2)+1:underscores(3)-1)}];
    gene_fluor_cell = [gene_fluor_cell{:} {project(underscores(3)+1:underscores(4)-1)}];
    % parameters for plots
    load([ReadPath 'nucleus_struct_ctrl.mat']);
    master_struct(i).nucleus_struct_ctrl = nucleus_struct_ctrl;
end

%%% compare delta distributions
pt_sp_cell = cell(1,numel(master_struct));
pt_nn_cell = cell(1,numel(master_struct));
pt_dt_cell = cell(1,numel(master_struct));
pt_time_cell = cell(1,numel(master_struct));
for i = 1:numel(master_struct)
    ctrl_vec = [master_struct(i).nucleus_struct_ctrl.ctrl_flags_final];
    dist_vec = [master_struct(i).nucleus_struct_ctrl.edgeDistSpot]*master_struct(i).nucleus_struct_ctrl(1).PixelSize;
    pt_vec_spot = [master_struct(i).nucleus_struct_ctrl.pt_spot];
    pt_sp_cell{i} = pt_vec_spot(ctrl_vec==1&dist_vec>dist_lim);
    pt_vec_null = [master_struct(i).nucleus_struct_ctrl.pt_null];
    pt_nn_cell{i} = pt_vec_null(ctrl_vec==1&dist_vec>dist_lim);        
end

delta_bins = linspace(-2,2,100);
% make comparison histogram
delta_fig = figure;
hold on
for i = 1:numel(pt_sp_cell)
    pt_ct = histc((pt_sp_cell{i}-pt_nn_cell{i}) ./ mean(pt_nn_cell{1}),delta_bins);
    pt_ct = pt_ct / sum(pt_ct);
    bar(delta_bins,pt_ct,1,'FaceAlpha',.4,'EdgeAlpha',0)
end
xlabel('normalized \DeltaF','Fontsize',14)
ylabel('share','Fontsize',14)
grid on
legend([protein_name_cell{1} '-' gene_name_cell{1}],[protein_name_cell{2} '-' gene_name_cell{2}])
saveas(delta_fig,[FigPath 'delta_comparison_' protein_name_cell{1} '-' gene_name_cell{1} '_' protein_name_cell{2} '-' gene_name_cell{2} '.png'])

for i = 1:2
    dt_fig = figure;
    hold on
    delta_vec_norm = (pt_sp_cell{i}-pt_nn_cell{i}) ./ mean(pt_nn_cell{i});
    pd = fitdist(delta_vec_norm','Normal');    
    pt_ct = histc(delta_vec_norm,delta_bins);    
    pt_ct = pt_ct / sum(pt_ct);
    bar(delta_bins,pt_ct,1,'FaceAlpha',.8,'EdgeAlpha',0)
    norm_vec = exp(-.5*((delta_bins-pd.mu)/pd.sigma).^2);
    plot(delta_bins,norm_vec/sum(norm_vec),'LineWidth',2) 
    xlabel('normalized \DeltaF','Fontsize',14)
    ylabel('share','Fontsize',14)
    grid on
    legend([protein_name_cell{i} '-' gene_name_cell{i}],['gaussian fit (\mu=' num2str(round(pd.mu,2)) ' \sigma=' num2str(round(pd.sigma,2))])
    saveas(dt_fig,[FigPath 'delta_' protein_name_cell{i} '-' gene_name_cell{i}  '.png'])
end


%% Make exmple unimodal and bimodal plots
n1 = normrnd(1,.7,1,10000);
n2 = normrnd(3.5,.7,1,2000);

n1_fig = figure;
histogram(n1,'Normalization','probability')
saveas(n1_fig,[FigPath 'n1_ex_norm.pdf'])

n2_fig = figure;
histogram([n1 n2],'Normalization','probability')
saveas(n2_fig,[FigPath 'n2_ex_norm.pdf'])
