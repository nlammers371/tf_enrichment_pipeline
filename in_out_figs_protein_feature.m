clear
close all
% script to compare Dl enrichment correlations hbP2P and snaBAC
project_cell = {'Dl-Ven x snaBAC','Dl-Ven x hbP2P'};
dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
figPath = [dropboxFolder '\LocalEnrichmentFigures\Dl_sna_hbP2P\'];
mkdir(figPath);
master_struct = struct;
for i = 1:numel(project_cell)
    dataPath = [dropboxFolder '\ProcessedEnrichmentData\' project_cell{i} '\'];
    load([dataPath 'input_output_results_protein.mat'])
    master_struct(i).results = results_struct;
end
% Define some colors
yw = [234 194 100]/256; % yellow
bl = [115 143 193]/256; % blue
rd = [213 108 85]/256; % red
gr = [191 213 151]/256; % green
br = [207 178 147]/256; % brown
% useful vectors
n_lags = floor(numel(master_struct(1).results(1).spot_protein_mean)/2);
x_axis = (-n_lags:n_lags)*20 / 60;
% make comparison plots
for i = 1:numel(master_struct(1).results)
    ID = master_struct(1).results(i).ID;
    target_protein = master_struct(1).results(i).spot_protein_mean;
    
    target_fluo = master_struct(1).results(i).response_mean;
    target_fluo_se = master_struct(1).results(i).response_ste;
    
    bio_control_fluo = master_struct(2).results(i).response_mean;
    bio_control_fluo_se = master_struct(2).results(i).response_ste;
        
    % perform simple linear fits to remove offset
    X = [ones(numel(x_axis),1) x_axis'];    
    bc = X\bio_control_fluo';
    bt = X\target_fluo';    
    bio_control_fluo = bio_control_fluo - bc(1) - bc(2)*x_axis;
    target_fluo = target_fluo - bt(1) - bt(2)*x_axis;
    % generate upper and lower bound vectors 
    ubt = target_fluo + target_fluo_se;
    lbt = target_fluo - target_fluo_se;
      
    % make figure
    in_out_fig = figure;
    hold on
    % error regions    
%     fill([x_axis fliplr(x_axis)],[ubcb fliplr(lbcb)],bl,'FaceAlpha',.1,'EdgeAlpha',0);
    fill([x_axis fliplr(x_axis)],[ubt fliplr(lbt)],rd,'FaceAlpha',.3,'EdgeAlpha',0);
    % trend plots        
    p2 = plot(x_axis,bio_control_fluo,'Color',bl,'LineWidth',1.5);    
    p1 = plot(x_axis,target_fluo,'Color',rd,'LineWidth',1.5);
    ax = gca;
    ax.YColor = 'black';
    ylabel('sna activity (au)')
    grid on
    yyaxis right
    p3 = plot(x_axis,target_protein,'-','Color','black','LineWidth',1.3);
    ax = gca;
    ax.YColor = 'black';
    set(gca,'ytick',[])
    ylabel('response (au)')
    legend([p1 p2 p3],'target (sna)','control (hbP2P)','Dorsal','Location','southeast')
    title(ID)
    % save
    saveas(in_out_fig,[figPath ID '.png'])
end

