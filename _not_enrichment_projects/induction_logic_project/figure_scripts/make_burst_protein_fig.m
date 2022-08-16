% Script to generate figures establishing presence of enrichment bursts at
% start of transcription bursts
clear 
close all
addpath('../utilities')
% set ID variables
project = 'Dl-Ven_snaBAC-mCh_v3';
DropboxFolder = 'S:\Nick\Dropbox (Personal)\';
[~, DataPath, FigureRoot] =   header_function(DropboxFolder, project); 
% define HMM parameters
K = 3;
w = 7;
% load data structure
load([DataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '_dt.mat'],'hmm_input_output')
% FigPath = [FigureRoot '\' project '\burst_alginment_fig\'];
FigPath = [FigureRoot  project '\RO1_grant_figs\'];
mkdir(FigPath)

%% pull illustrative trace, HMM trajectory, and protein
close all
for i = randsample(1:numel(hmm_input_output),numel(hmm_input_output),false)
    plot(hmm_input_output(i).time,hmm_input_output(i).fluo_check)
    title(num2str(i))
    pause(1.5)
end
%%
Tres = 20;
mSize = 45;
cmap1 = brewermap([],'Set2');
blue = [115 143 193]/256;
plot_id_vec = [900 941 4208];
% plot_id = 941;
PBoC_flag_vec = [false true];

% define colors 
red = [213 108 85]/256;
green = [122 169 116]/256;


for p = 1:length(plot_id_vec)
    close all
    plot_id = plot_id_vec(p);
    
    % generate vector of burst-specific loading rates (this is for illustrative
    % purposes only)
    r_vec = hmm_input_output(plot_id).r_vec;  
    z_vec = (hmm_input_output(plot_id).z_vec'-1)>0;  
    z_diff_vec = [0 diff(z_vec)];
    z_chpts = find(z_diff_vec ~=0);
    if z_vec(1) == 1
        z_chpts = [1 z_chpts];
    end
    if z_vec(end) == 1
         z_chpts = [z_chpts numel(z_vec)];
    end
    init_vec = zeros(size(z_vec));
    for c = 1:2:numel(z_chpts)
        init_vec(z_chpts(c):z_chpts(c+1)-1) = nanmean(r_vec(z_chpts(c):z_chpts(c+1)-1))*Tres;
    end    
    
    for i = 1:2
        PBoC_flag = PBoC_flag_vec(i);
        % raw fluorescence
        fluo = hmm_input_output(plot_id).fluo;
        time = hmm_input_output(plot_id).time/60;
        % Protein channel
        spot_protein = hmm_input_output(plot_id).mf_protein;%*PixelSize;
        nn_ids = find(~isnan(spot_protein)&~isnan(fluo));        

        input_output_fig = figure;
        hold on
        p = plot(0,0);
        p1 = plot(time,fluo,'-','Color',red,'LineWidth',1.5');
        scatter(time,fluo,mSize,'MarkerFaceColor',red,'MarkerEdgeAlpha',0);
        ylabel('snail fluorescence (AU)')
        
        yyaxis right
                
%         area(time,spot_protein,'FaceColor',green);
        p2 = plot(time,spot_protein,'-','Color',green,'LineWidth',1.5');
        scatter(time,spot_protein,mSize,'MarkerFaceColor',green,'MarkerEdgeAlpha',0);
    
        xlabel('time (minutes)')
        ylabel('Dorsal-venus fluorescence (AU)')                
        
%         legend([p1 p2],'snail transcription','Dorsal concentration')
        box on
        xlim([time(nn_ids(1)) time(nn_ids(end))])
        if PBoC_flag
            suffix = 'PBoC';
            StandardFigurePBoC(p,gca);
            input_output_fig.InvertHardcopy = 'off';
        else
            suffix = 'standard';
            StandardFigure(p,gca);
        end
        
        ax = gca;
        ax.YAxis(1).Color = 'k';
        ax.YAxis(2).Color = 'k';
        ax.XAxis(1).Color = 'k';
%         ax.XAxis(2).Color = 'k';
        
        saveas(input_output_fig,[FigPath 'fluo_trend_' sprintf('%03d',plot_id) suffix '.pdf'])
        saveas(input_output_fig,[FigPath 'fluo_trend_' sprintf('%03d',plot_id) suffix '.png'])
        
        % hmm-decoded
      
        hmm_fig = figure('Position',[100 100 512 256]);
        hold on
        p = plot(0,0);
        s = stairs(time,init_vec,'-','Color',blue,'LineWidth',1.5);    
        xlabel('time (minutes)')
        ylabel('initiation rate (AU)')
        xlim([time(nn_ids(1)) time(nn_ids(end))])
        ylim([0 1.05*max(init_vec(nn_ids))])
        set(gca,'Ytick',[0:5:25])
        box on
        if PBoC_flag        
            StandardFigurePBoC(p,gca);
            hmm_fig.InvertHardcopy = 'off';
        else        
            StandardFigure(p,gca);
        end
        ax = gca;
        ax.YAxis(1).Color = 'k';        
        ax.XAxis(1).Color = 'k';
        
        saveas(hmm_fig,[FigPath 'hmm_trend_' sprintf('%03d',plot_id) suffix '.pdf'])
        saveas(hmm_fig,[FigPath 'hmm_trend_' sprintf('%03d',plot_id) suffix '.png']) 
                          
    end
end   
