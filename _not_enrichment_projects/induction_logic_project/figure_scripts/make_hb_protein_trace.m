clear 
close all
addpath('utilities')

% set ID variables
project = '2xDl-Ven_hbP2P-mCh';

% Params
fluo_dim = 2;
protein_dim = 2;
K = 3;
w = 7;

% paths
DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\';
[~, DataPath, FigureRoot] =   header_function(DropboxFolder, project); 
FigPath = [FigureRoot '\' project '\input_output01\'];
mkdir(FigPath)

% load data 
load([DataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '_f' num2str(fluo_dim) 'D_p' num2str(protein_dim) 'D.mat'])

%% see if we can find trace that stretches 40 min
long_trace_ids = [];
for i = 1:length(hmm_input_output)
  if hmm_input_output(i).time(end) - hmm_input_output(i).time(1) >= 25*60
    long_trace_ids(end+1) = i;
  end
end

plot_id = long_trace_ids(2);

pt_fig = figure;
hold on
cmap1 = brewermap([],'Set2');

plot(hmm_input_output(plot_id).time/60,hmm_input_output(plot_id).spot_protein,'-','Color',cmap1(3,:),'LineWidth',2);
ylabel('local [Dorsal] (au)')

p = plot(0,0);

ax = gca;
ax.YColor = 'black';

% grid on
xlabel('time (min)')
legend('Dl at {\it hbP2P}', 'Location','northwest');

set(gca,'Fontsize',14);%,'xtick',-4:2:4)
chH = get(gca,'Children');
set(gca,'Children',flipud(chH));

ylim([0 5])
xlim([0 40])

set(gca,    'Box','off',...
            'Color',[228,221,209]/255,...            
            'TickLength',[0.02,0.05])    
pt_fig.Color = 'white';        
pt_fig.InvertHardcopy = 'off';

% save
saveas(pt_fig,[FigPath 'Dorsal_at_hbP2P_control.tif'])
saveas(pt_fig,[FigPath 'Dorsal_at_hbP2P_control.pdf'])
