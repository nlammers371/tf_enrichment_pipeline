clear 
close all

dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
project_cell = {'Dl-Ven_snaBAC-mCh','Dl-Ven_hbP2P-mCh'};
figPath = [dropboxFolder 'LocalEnrichmentFigures\input_output_dl_sna\'];
mkdir(figPath)

master_struct = struct;
for i = 1:2
    project = project_cell{i};
    dataPath = [dropboxFolder 'ProcessedEnrichmentData\' project '/'];
    load([dataPath 'input_output_results.mat'])
    master_struct(i).results_struct = results_struct;
end

%5% Make figures
fluo_activity_array = master_struct(1).results_struct.fluo_activity_array;
spot_protein_array = master_struct(1).results_struct.spot_protein_array;
swap_protein_array = master_struct(1).results_struct.swap_protein_array;
virtual_protein_array = master_struct(1).results_struct.virtual_protein_array;
off_spot_protein_array = master_struct(2).results_struct.spot_protein_array;

window_vec = (1:size(fluo_activity_array,2)) - floor(size(fluo_activity_array,2)/2)-1;

% make rise series figs
yw = [234 194 100]/256; % yellow
bl = [115 143 193]/256; % blue
rd = [213 108 85]/256; % red
gr = [191 213 151]/256; % green
br = [207 178 147]/256; % brown

% plot_range = 1:11;
for pi = 1:size(fluo_activity_array,1)
%     pi = plot_range(p);
    fluo_fig = figure;
    hold on
    yyaxis left
    plot(window_vec*20/60,fluo_activity_array(pi,:,6),'-','Color','black','LineWidth',1.3)
    ylim([-50 50])
    ylabel('transcriptional activity (au)')
    yyaxis right
    plot(window_vec*20/60,spot_protein_array(pi,:,6),'-','Color',bl,'LineWidth',1.3)
    plot(window_vec*20/60,swap_protein_array(pi,:,6),'-','Color',rd,'LineWidth',1)
    plot(window_vec*20/60,virtual_protein_array(pi,:,6),'-','Color',gr,'LineWidth',1)
    plot(window_vec*20/60,off_spot_protein_array(pi,:,6),'-','Color',yw,'LineWidth',1)
    ylabel('Dl concentration (au)')
    ylim([-12 12])
    xlabel('offset (minutes)')
    legend('transcription','target locus (sna)','control (swap)','control (virtual spot)','control (hbP2P)','Location','northeast')    
    saveas(fluo_fig,[figPath 'fluo_feature_d' num2str(round(100*pi/size(fluo_activity_array,1))) '.png'])
end    