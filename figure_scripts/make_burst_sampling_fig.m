% Script to generate figures establishing presence of enrichment bursts at
% start of transcription bursts
clear 
% close all
addpath('utilities')
% set ID variables
project = 'Dl-Ven_snaBAC-mCh';
DropboxFolder = 'E:\Nick\LivemRNA\Dropbox (Personal)\';
[~, DataPath, FigureRoot] =   header_function(DropboxFolder, project); 
% define HMM parameters
K = 3;
w = 7;
% load data structure
load([DataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '_dt.mat'],'hmm_input_output')
load([DataPath 'nucleus_struct.mat'])
FigPath = [FigureRoot '\' project '\burst_alginment_fig\'];
mkdir(FigPath)

%% pull illustrative trace, HMM trajectory, and protein
close all
% for i = randsample(1:numel(hmm_input_output),numel(hmm_input_output),false)
%     plot(hmm_input_output(i).time,hmm_input_output(i).fluo_check)
%     title(num2str(i))
%     pause(1.5)
% end
PixelSize = nucleus_struct(1).PixelSize;
zStep = nucleus_struct(1).zStep;
VoxelSize = PixelSize^2 * zStep;
cmap1 = brewermap([],'Set2');
blue = [115 143 193]/256;
plot_id = 1068;
PBoC_flag_vec = [false true];
for i = 1:2
    PBoC_flag = PBoC_flag_vec(i);
    % raw fluorescence
    fluo = hmm_input_output(plot_id).fluo_check;
    time = hmm_input_output(plot_id).time/60;
    
    fluo_fig = figure;
    hold on
    p = plot(0,0);
    plot(time,fluo,'--','Color','black','LineWidth',1.5');
%     scatter(time,fluo,'MarkerFaceColor','black','MarkerEdgeColor','black')
    xlabel('time (minutes)')
    ylabel('fluorescence (AU)')
    box on
    xlim([10,max(time)])
    if PBoC_flag
        suffix = '_PBoC';
        StandardFigurePBoC(p,gca);
        fluo_fig.InvertHardcopy = 'off';
    else
        suffix = '_standard';
        StandardFigure(p,gca);
    end
    saveas(fluo_fig,[FigPath 'fluo_trend_' suffix '.pdf'])
    saveas(fluo_fig,[FigPath 'fluo_trend_' suffix '.png'])
    
    % raw fluorescence
    z_vec = (hmm_input_output(plot_id).z_vec-1)>0;    
    hmm_fig = figure;
    hold on
    p = plot(0,0);
    s = stairs(time,z_vec,'-','Color',blue,'LineWidth',1.5');    
    xlabel('time (minutes)')
    ylabel('promoter state')
    xlim([10,max(time)])
    set(gca,'Ytick',[0 1])
    box on
    if PBoC_flag        
        StandardFigurePBoC(p,gca);
        hmm_fig.InvertHardcopy = 'off';
    else        
        StandardFigure(p,gca);
    end
    saveas(hmm_fig,[FigPath 'hmm_trend_' suffix '.pdf'])
    saveas(hmm_fig,[FigPath 'hmm_trend_' suffix '.png'])
    
    % Protein channel
    spot_protein = hmm_input_output(plot_id).spot_protein_dt*PixelSize;
    nn_ids = find(~isnan(spot_protein));
    protein_fig = figure;
    hold on
    p = plot(0,0);
    plot(time(nn_ids),spot_protein(nn_ids),'Color',cmap1(2,:),'LineWidth',1.5);
    xlabel('time (minutes)')
    ylabel('Dl concentration (AU)')    
    xlim([10,max(time)])
    ylim([-30 30])
    box on
    if PBoC_flag        
        StandardFigurePBoC(p,gca);
        protein_fig.InvertHardcopy = 'off';
    else        
        StandardFigure(p,gca);
    end
    saveas(protein_fig,[FigPath 'protein_trend_' suffix '.pdf'])
    saveas(protein_fig,[FigPath 'protein_trend_' suffix '.png'])    
end
    
