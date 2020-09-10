clear 
close all
addpath('../utilities')

% set ID variables and paths

project = 'Rbp1-GFP_snaBAC-mCh';
% DropboxFolder = 'S:\Nick\Dropbox\';
DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\';
[~, DataPath, FigureRoot] =   header_function(DropboxFolder, project); 
FigPath = [FigureRoot '\' project '\exploratory_analyses\'];
mkdir(FigPath)

% load data
load([DataPath 'nucleus_struct_protein.mat'])

%% plot MS2 fluorescence vs Pol II
close all
ms2_fluo_vec = [nucleus_struct_protein.fluo];
polII_fluo_vec = [nucleus_struct_protein.spot_protein_vec];

figure(1);
cmap = brewermap([],'Set2');
scatter(ms2_fluo_vec,polII_fluo_vec,'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)

set(gca,'Fontsize',14)
xlabel('MS2-MCP-mCherry spot intensity (au)')
ylabel('local Rbp1-GFP intensity (au)')
grid on
xlim([0 2e5])
saveas(gcf,[FigPath 'xcal_scatter.png'])


% try to estiamte lower bound
% find subset of points that are in lowest 3 percent for their flourescence
% bin
close all
fluo_index = linspace(prctile(ms2_fluo_vec,5),prctile(ms2_fluo_vec,97.5),50);
low_filter = false(size(polII_fluo_vec));
for f = 1:length(fluo_index)-1
  bin_filter = ms2_fluo_vec>=fluo_index(f)&ms2_fluo_vec<fluo_index(f+1);
  pct = prctile(polII_fluo_vec(bin_filter),3);
  low_filter(bin_filter & polII_fluo_vec <= pct) = true;
end

% fit line to these low points
mdl = fitlm(ms2_fluo_vec(low_filter),polII_fluo_vec(low_filter));
ms2_pd_vec = [0 sort(ms2_fluo_vec(low_filter)) max(ms2_fluo_vec)]';
polII_pred = predict(mdl,ms2_pd_vec);

figure(2);
hold on
cmap = brewermap([],'Set2');
scatter(ms2_fluo_vec,polII_fluo_vec,'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
scatter(ms2_fluo_vec(low_filter),polII_fluo_vec(low_filter),'MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
plot(ms2_pd_vec,polII_pred,'Color','k','LineWidth',1.5)

set(gca,'Fontsize',14)
xlabel('MS2-MCP-mCherry spot intensity (au)')
ylabel('local Rbp1-GFP intensity (au)')
grid on
ylim([0 1.5e4])
xlim([0 2e5])
saveas(gcf,[FigPath 'lower_bound_fit.png'])

% %%
% close all
% i = 2;
% figure(2);
% hold on
% yyaxis left
% plot(nucleus_struct_protein(i).time,nucleus_struct_protein(i).fluo)
% 
% yyaxis right
% plot(nucleus_struct_protein(i).time,nucleus_struct_protein(i).spot_protein_vec)
