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
    load([dataPath 'hmm_input_output_results.mat'])
    master_struct(i).results_struct = results_struct;
end

%5% Make figures
fluo_activity_array = master_struct(1).results_struct(1).fluo_quantile_array;
hmm_activity_array = master_struct(1).results_struct(1).hmm_quantile_array;
spot_protein_array = master_struct(1).results_struct(1).spot_quantile_array;
swap_protein_array = master_struct(1).results_struct(1).swap_quantile_array;
virtual_protein_array = master_struct(1).results_struct(1).virtual_quantile_array;
off_spot_protein_array = master_struct(2).results_struct(1).spot_quantile_array;
% get lag and lead keys
z_dur_key = master_struct(1).results_struct(1).z_dur_key;
z_sep_key = master_struct(1).results_struct(1).z_sep_key;

window_vec = (1:size(fluo_activity_array,2)) - floor(size(fluo_activity_array,2)/2)-1;
quantile_plot_vec = [6];
% make rise series figs
yw = [234 194 100]/256; % yellow
bl = [115 143 193]/256; % blue
rd = [213 108 85]/256; % red
gr = [191 213 151]/256; % green
br = [207 178 147]/256; % brown
% y_range_cell = {[-10 10],[-12 12],[-12 12]};
% plot_range = 1:11;
summary_sep = 6;
for qi = quantile_plot_vec
    qi_string = ['hmm_q' sprintf('%02d',round(100*qi/11))];
    qi_path = [figPath qi_string filesep];
    mkdir(qi_path)
    for pi = 1:size(fluo_activity_array,1)   
        z_size = z_dur_key(pi);
        min_sep = z_sep_key(pi);
        feature_string = 'rise';
        y_bounds = [-25 25];
        if sign(z_size)== -1
            feature_string = 'fall';
            y_bounds = [-20 20];
        end
        fluo_fig = figure('Visible','off');
        hold on
        yyaxis left
        fluo_vec = fluo_activity_array(pi,:,qi);
        fluo_vec = fluo_vec-min(fluo_vec);
        fluo_vec = fluo_vec/max(fluo_vec);
        plot(window_vec*20/60,hmm_activity_array(pi,:,qi),'-','Color','black','LineWidth',1.3)
        plot(window_vec*20/60,fluo_vec,'--','Color',[.3 .3 .3],'LineWidth',1.3)
        ylim([0 1])
        ax = gca;
        ax.YColor = 'black';        
        ylabel('transcriptional activity (au)')
        
        % subtract mean
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
        ax = gca;
        ax.YColor = 'black';
        ylim(y_bounds);
        xlabel('offset (minutes)')
        ax.FontSize = 12;
%         legend('transcription','target locus (sna)','control (swap)','control (virtual spot)','control (hbP2P)','Location','northeast')    
       saveas(fluo_fig,[qi_path 'hmm_' feature_string '_' qi_string '_sz' sprintf('%02d',abs(z_size)) '_sep' sprintf('%02d',abs(min_sep)) '.tif'])
    end    
    %%% make summary figures
    % rises
    cm = flipud(jet(128));
    rise_ids = find(abs(z_dur_key)<15&z_dur_key>0&z_sep_key==summary_sep);
    hmm_rise_fig = figure('Visible','off');
    hold on
    yyaxis left    
    area(window_vec*20/60,hmm_activity_array(rise_ids(end),:,qi),'FaceColor',cm(128,:),'FaceAlpha',.15);%'-','Color','black','LineWidth',1.3)
    area(window_vec*20/60,hmm_activity_array(rise_ids(1),:,qi),'FaceColor',cm(1,:),'FaceAlpha',.25);
    ylim([0 1])
    ax = gca;
    ax.YColor = 'black';        
    ylabel('transcriptional activity (au)')
    
    inc = floor(128/numel(rise_ids));
    yyaxis right
    spot_vec = spot_protein_array(rise_ids(1),:,qi) - nanmean(spot_protein_array(rise_ids(1),:,qi));
    plot(window_vec*20/60,spot_vec,'-','Color',cm(1,:),'LineWidth',2)
    for j = 2:numel(rise_ids)-1
        spot_vec = spot_protein_array(rise_ids(j),:,qi) - nanmean(spot_protein_array(rise_ids(j),:,qi));
        plot(window_vec*20/60,spot_vec,'-','Color',[cm(1+inc*j,:) .5],'LineWidth',1)
    end
    spot_vec = spot_protein_array(rise_ids(end),:,qi) - nanmean(spot_protein_array(rise_ids(end),:,qi));
    plot(window_vec*20/60,spot_vec,'-','Color',cm(end,:),'LineWidth',2)
    
    ylabel('Dl concentration (au)')
    ax = gca;
    ax.YColor = 'black';
    ylim([-33 33]);
    xlabel('offset (minutes)')
    ax.FontSize = 12;
    saveas(hmm_rise_fig,[figPath 'hmm_rise_summary' qi_string '_sep' sprintf('%02d',abs(summary_sep)) '.tif'])
    
    % falls
    cm = jet(128);
    fall_ids = find(abs(z_dur_key)<15&z_dur_key<0&z_sep_key==summary_sep);
%     fall_ids = fall_ids(2:end)
    hmm_fall_fig = figure('Visible','off');
    hold on
    yyaxis left    
    area(window_vec*20/60,hmm_activity_array(fall_ids(1),:,qi),'FaceColor',cm(1,:),'FaceAlpha',.15);
    area(window_vec*20/60,hmm_activity_array(fall_ids(end),:,qi),'FaceColor',cm(128,:),'FaceAlpha',.25);%'-','Color','black','LineWidth',1.3)    
    ylim([0 1])
    ax = gca;
    ax.YColor = 'black';        
    ylabel('transcriptional activity (au)')
    
    inc = floor(128/numel(fall_ids));
    yyaxis right
    spot_vec = spot_protein_array(fall_ids(1),:,qi) - nanmean(spot_protein_array(fall_ids(1),:,qi));
    plot(window_vec*20/60,spot_vec,'-','Color',cm(1,:),'LineWidth',2)
    for j = 2:numel(fall_ids)-1
        spot_vec = spot_protein_array(fall_ids(j),:,qi) - nanmean(spot_protein_array(fall_ids(j),:,qi));
        plot(window_vec*20/60,spot_vec,'-','Color',[cm(1+inc*j,:) .5],'LineWidth',1)
    end
    spot_vec = spot_protein_array(fall_ids(end),:,qi) - nanmean(spot_protein_array(fall_ids(end),:,qi));
    plot(window_vec*20/60,spot_vec,'-','Color',cm(end,:),'LineWidth',2)
    
    ylabel('Dl concentration (au)')
    ax = gca;
    ax.YColor = 'black';
    ylim([-20 20]);
    xlabel('offset (minutes)')
    ax.FontSize = 12;
    saveas(hmm_fall_fig,[figPath 'hmm_fall_summary' qi_string '_sep' sprintf('%02d',abs(summary_sep)) '.tif'])
end