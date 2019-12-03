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
FigPath = [FigureRoot '\' project '\burst_decoding_fig\'];
mkdir(FigPath)

%% pull illustrative trace, HMM trajectory, and protein
close all
for i = randsample(1:numel(hmm_input_output),numel(hmm_input_output),false)
    plot(hmm_input_output(i).time,hmm_input_output(i).fluo_check)
    title(num2str(i))
    pause(1.5)
end

%%
nc_pt_index = [nucleus_struct.ParticleID];
PixelSize = nucleus_struct(1).PixelSize;
zStep = nucleus_struct(1).zStep;
VoxelSize = PixelSize^2 * zStep;
cmap1 = brewermap([],'Set2');
blue = [115 143 193]/256;
red = [213 108 85]/256;
% plot_id = 1068;
plot_id = 420;
PBoC_flag_vec = [false true];
Tres = hmm_input_output(1).Tres;

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


% extract other data vectors
fluo = hmm_input_output(plot_id).fluo_check;
fluo_hmm = hmm_input_output(plot_id).fluo_hmm;
fluo_err = repelem(5,numel(fluo)); %NL: use fake errror temporarily
time = hmm_input_output(plot_id).time/60;

for i = 1%1:2
    PBoC_flag = PBoC_flag_vec(i);
   
    
    fluo_fig = figure('Position',[0 0 512 256]);
    hold on
    p = plot(0,0);
    e = errorbar(time,fluo,fluo_err,'-','Color','black','LineWidth',1.5');    
    e.CapSize = 0;
    scatter(time,fluo,'MarkerFaceColor','black','MarkerEdgeAlpha',0)
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
    
    % now include fit
    fluo_hmm_fig = figure('Position',[0 0 512 256]);
    hold on
    p = plot(0,0);
    e = errorbar(time,fluo,fluo_err,'-','Color','black','LineWidth',1.5);    
    e.CapSize = 0;    
    scatter(time,fluo,'MarkerFaceColor','black','MarkerEdgeAlpha',0)
    plot(time,fluo_hmm,'Color',red,'LineWidth',1.5)
%     scatter(time,fluo,'MarkerFaceColor','black','MarkerEdgeColor','black')
    xlabel('time (minutes)')
    ylabel('fluorescence (AU)')
    box on
    xlim([10,max(time)])
    if PBoC_flag
        suffix = '_PBoC';
        StandardFigurePBoC(p,gca);
        fluo_hmm_fig.InvertHardcopy = 'off';
    else
        suffix = '_standard';
        StandardFigure(p,gca);
    end
    saveas(fluo_hmm_fig,[FigPath 'fluo_trend_fit_' suffix '.pdf'])
    saveas(fluo_hmm_fig,[FigPath 'fluo_trend_fit_' suffix '.png'])
    
    % hmm-decoded
      
    hmm_fig = figure('Position',[0 0 512 256]);
    hold on
    p = plot(0,0);
    s = stairs(time,init_vec,'-','Color',blue,'LineWidth',1.5);    
    xlabel('time (minutes)')
    ylabel('initiation rate (AU)')
    xlim([10,max(time)])
    set(gca,'Ytick',[0:5:25])
    box on
    if PBoC_flag        
        StandardFigurePBoC(p,gca);
        hmm_fig.InvertHardcopy = 'off';
    else        
        StandardFigure(p,gca);
    end
    saveas(hmm_fig,[FigPath 'hmm_trend_' suffix '.pdf'])
    saveas(hmm_fig,[FigPath 'hmm_trend_' suffix '.png'])      
end
    
