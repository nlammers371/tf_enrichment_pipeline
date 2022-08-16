% script to make basic performance figures
clear
close all

addpath(genpath('../utilities'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%% Load experimental data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

resultsRoot = 'S:\Nick\Dropbox\ProcessedEnrichmentData\parameterSweeps\';
if ~isfolder(resultsRoot)
  resultsRoot = 'C:\Users\nlamm\Dropbox (Personal)\ProcessedEnrichmentData\parameterSweeps\';
end

% mak figure directory  
date_string = '22-Sep-2021';
resultsPath = [resultsRoot date_string filesep];

try
    FigurePath = [liveProject.figurePath 'io_sim_results' filesep date_string filesep];
catch
    FigurePath = 'C:\Users\nlamm\Dropbox (Personal)\Postdoc\UChicago\Figures\';
end
mkdir(FigurePath)

% load data
load([resultsRoot 'io_ref_wt.mat'],'io_ref_wt')
%%
plot_id = 11;
gr = [122 169 116]/255;

% de-mean kni trace for illustrative purposes
kni_trace = io_ref_wt.knirps_array(:,plot_id);
kni_dm = detrend(kni_trace,3)+1;

kni_fig = figure('Position',[100 100 512 256]);

plot(io_ref_wt.time_axis/60,kni_dm,'LineWidth',3,'Color',gr)

xlabel('time');
ylabel('Knirps concentration (au)')

set(gca,'FontSize',14)

box on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
%     xlim(ap_lims)
still_on_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
saveas(kni_fig,[FigurePath 'kni_trace.pdf'])