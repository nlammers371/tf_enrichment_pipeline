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
    load([dataPath 'fluo_input_output_results.mat'])
    master_struct(i).results_struct = results_struct;
end

%5% Make figures
fluo_activity_array = master_struct(1).results_struct(1).fluo_activity_array;
spot_protein_array = master_struct(1).results_struct(1).spot_protein_array;
swap_protein_array = master_struct(1).results_struct(1).swap_protein_array;
virtual_protein_array = master_struct(1).results_struct(1).virtual_protein_array;
off_spot_protein_array = master_struct(2).results_struct(1).spot_protein_array;

window_vec = (1:size(fluo_activity_array,2)) - floor(size(fluo_activity_array,2)/2)-1;
quantile_plot_vec = [3 6 9];
% make rise series figs
yw = [234 194 100]/256; % yellow
bl = [115 143 193]/256; % blue
rd = [213 108 85]/256; % red
gr = [191 213 151]/256; % green
br = [207 178 147]/256; % brown
y_range_cell = {[-10 10],[-12 12],[-12 12]};
% plot_range = 1:11;
iter = 1;
for qi = quantile_plot_vec
    qi_string = ['q' sprintf('%02d',round(100*qi/11))];
    qi_path = [figPath qi_string filesep];
    mkdir(qi_path)
    for pi = 1:-1%size(fluo_activity_array,1)
    %     pi = plot_range(p);
        fluo_fig = figure('Visible','off');
        hold on
        yyaxis left
        plot(window_vec*20/60,fluo_activity_array(pi,:,qi),'-','Color','black','LineWidth',1.3)
        ylim([-50 50])
        ylabel('transcriptional activity (au)')
        spot_vec = spot_protein_array(pi,:,qi) - nanmean(spot_protein_array(pi,:,qi));
        swap_vec = swap_protein_array(pi,:,qi) - nanmean(swap_protein_array(pi,:,qi));
        virtual_vec = virtual_protein_array(pi,:,qi) - nanmean(virtual_protein_array(pi,:,qi));
        off_vec = off_spot_protein_array(pi,:,qi) - nanmean(off_spot_protein_array(pi,:,qi));
        yyaxis right
        plot(window_vec*20/60,spot_vec,'-','Color',bl,'LineWidth',2)
        plot(window_vec*20/60,swap_vec,'-','Color',[rd .5],'LineWidth',1)
        plot(window_vec*20/60,virtual_vec,'-','Color',[gr .5],'LineWidth',1)
        plot(window_vec*20/60,off_vec,'-','Color',[yw .5],'LineWidth',1)
        ylabel('Dl concentration (au)')
        ylim([y_range_cell{iter}]);
        xlabel('offset (minutes)')
%         legend('transcription','target locus (sna)','control (swap)','control (virtual spot)','control (hbP2P)','Location','northeast')    
        saveas(fluo_fig,[qi_path 'fluo_feature_' qi_string '_d' num2str(round(100*pi/size(fluo_activity_array,1))) '.tif'])
    end       
    
    %%% make summary figures
    % rises
    cm = flipud(jet(128));    
    fluo_summary_fig = figure('Visible','off');
    hold on
    yyaxis left    
    fluo_vec_one = fluo_activity_array(end,:,qi);
    fluo_vec_one = fluo_vec_one-min(fluo_vec_one);
    fluo_vec_one = fluo_vec_one/max(fluo_vec_one);
    plot(window_vec*20/60,fluo_vec_one,'Color',cm(128,:));%'-','Color','black','LineWidth',1.3)
    
    fluo_vec_two = fluo_activity_array(1,:,qi);
    fluo_vec_two = fluo_vec_two-min(fluo_vec_two);
    fluo_vec_two = fluo_vec_two/max(fluo_vec_two);
    
    area(window_vec*20/60,fluo_vec_one,'FaceColor',cm(128,:),'FaceAlpha',.1);%'-','Color','black','LineWidth',1.3)
    area(window_vec*20/60,fluo_vec_two,'FaceColor',cm(1,:),'FaceAlpha',.1);
    
    ylim([0 1])
    ax = gca;
    ax.YColor = 'black';        
    ylabel('transcriptional activity (au)')
    
    inc = floor(128/size(fluo_activity_array,1));
    yyaxis right
    spot_vec = spot_protein_array(1,:,qi) - nanmean(spot_protein_array(1,:,qi));
    plot(window_vec*20/60,spot_vec,'-','Color',cm(1,:),'LineWidth',2)
%     for j = round(size(fluo_activity_array,1)/2)
%         spot_vec = spot_protein_array(j,:,qi) - nanmean(spot_protein_array(j,:,qi));
%         plot(window_vec*20/60,spot_vec,'-','Color',[cm(1+inc*j,:) .5],'LineWidth',1)
%     end
    spot_vec = spot_protein_array(end,:,qi) - nanmean(spot_protein_array(end,:,qi));
    plot(window_vec*20/60,spot_vec,'-','Color',cm(end,:),'LineWidth',2)
    
    ylabel('Dl concentration (au)')
    ax = gca;
    ax.YColor = 'black';
    ylim([-10 15]);
    xlabel('offset (minutes)')
    ax.FontSize = 12;
    saveas(fluo_summary_fig,[figPath 'fluo_summary' qi_string '.tif'])
    iter = iter + 1;
end