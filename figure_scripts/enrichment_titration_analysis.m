clear
close all
addpath('utilities')

% set paths
DropboxFolder =  'E:\Nick\LivemRNA\Dropbox (Personal)\';
project1x = 'Dl-Ven_snaBAC-mCh';
project2x = '2xDl-Ven_snaBAC-mCh';
[~, DataPath1x, FigureRoot] =   header_function(DropboxFolder, project1x);
[~, DataPath2x, ~] =   header_function(DropboxFolder, project2x);
FigPath = [FigureRoot 'dorsal_sna_enrichment_titration/'];
mkdir(FigPath)
% load data
load([DataPath1x 'nucleus_struct_protein.mat'])
nc_data_1x = nucleus_struct_protein;
load([DataPath2x 'nucleus_struct_protein.mat'])
nc_data_2x = nucleus_struct_protein;
clear nucleus_struct_protein;

% basic plot and data qc params 
DistLim = 0.8; % min distance from edge permitted (um)
nBoots = 100;%00; % number of bootstrap samples to use for estimating SE

PixelSize = nc_data_1x(1).PixelSize;
zStep = nc_data_1x(1).zStep;
VoxelSize = PixelSize^2 * zStep;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pull useful vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% distance from nucleus edge
dist1x_vec = [nc_data_1x.spot_edge_dist_vec];
dist1x_ft_vec = dist1x_vec >= DistLim;
dist2x_vec = [nc_data_2x.spot_edge_dist_vec];
dist2x_ft_vec = dist2x_vec >= DistLim;

% absolute spot enrichment
delta_protein1x_vec = ([nc_data_1x.spot_protein_vec] - [nc_data_1x.edge_null_protein_vec])*VoxelSize;
delta_protein1x_dist_vec = delta_protein1x_vec(dist1x_ft_vec);
delta_protein2x_vec = [nc_data_2x.spot_protein_vec] - [nc_data_2x.edge_null_protein_vec];
delta_protein2x_dist_vec = delta_protein2x_vec(dist2x_ft_vec);

% average nuclear concentrations
mf_protein1x_vec = [nc_data_1x.mf_null_protein_vec]*VoxelSize;
mf_protein1x_dist_vec = mf_protein1x_vec(dist1x_ft_vec);
mf_protein2x_vec = [nc_data_2x.mf_null_protein_vec];
mf_protein2x_dist_vec = mf_protein2x_vec(dist2x_ft_vec);

%% Simple histograms to compare nucleus concentrations for 1 and 2x
n_bins = 100;
% get bounds
mf_lb = min([mf_protein1x_vec mf_protein2x_vec]);
mf_ub = max([mf_protein1x_vec mf_protein2x_vec]);
dt_lb = prctile([delta_protein1x_vec mf_protein2x_vec],1);
dt_ub = prctile([delta_protein1x_vec mf_protein2x_vec],99);
% generate bins
delta_bins = linspace(dt_lb,dt_ub,n_bins);
mf_bins = linspace(mf_lb,mf_ub,n_bins);

close all
% make figure2
mf_hist = figure;
hold on
histogram(mf_protein1x_dist_vec,mf_bins,'Normalization','probability')
histogram(mf_protein2x_dist_vec,mf_bins,'Normalization','probability')
xlabel('nuclear Dl concentration (au)')
ylabel('share')
legend('1x','2x')
grid on 
box on
set(gca,'Fontsize',14)
saveas(mf_hist,[FigPath 'nuclear_dorsal_1x_vs_2x.png'])
saveas(mf_hist,[FigPath 'nuclear_dorsal_1x_vs_2x.pdf'])


delta_hist = figure;
hold on
histogram(delta_protein1x_dist_vec,mf_bins,'Normalization','probability')
histogram(delta_protein2x_dist_vec,mf_bins,'Normalization','probability')
xlabel('Dl enrichment at locus (au)')
ylabel('share')
legend('1x','2x')
grid on 
box on
set(gca,'Fontsize',14)
saveas(mf_hist,[FigPath 'enrichment_1x_vs2x.png'])
saveas(mf_hist,[FigPath 'enrichment_1x_vs2x.pdf'])


%%% Look at enrichment as a function of nucleus concentration

% calculate average enrichment as a function of nuclear concentration
% combine sets for now
n_points = 25;
mf_lb = prctile([mf_protein1x_vec mf_protein2x_vec],.5);
mf_ub = prctile([mf_protein1x_vec mf_protein2x_vec],99.5);
mf_plot_index = linspace(mf_lb,mf_ub,n_points);

% concatenate 1x and 2x sets
mf_vec_full = [mf_protein1x_vec mf_protein2x_vec];
delta_vec_full = [delta_protein1x_vec delta_protein2x_vec];
id_vec_full = [repelem(1,numel(mf_protein1x_vec)) repelem(2,numel(mf_protein2x_vec))];

% initialize boot arrays
mf_sigma = median(diff(mf_plot_index))/2;
delta_boot_array = NaN(nBoots,n_points);
index_vec = 1:numel(mf_vec_full);
for n = 1:n_points
    mf_pt = mf_plot_index(n);
    mf_weights = exp(-.5*((mf_pt-mf_vec_full)/mf_sigma).^2); % gaussian weights
    for b = 1:nBoots
        boot_indices = randsample(index_vec,numel(index_vec),true);
        boot_weights = mf_weights(boot_indices);
        delta_boot_array(b,n) = nansum(delta_vec_full(boot_indices).*boot_weights) / nansum(boot_weights);
    end
end

%% plot
% calculate mean and se
mean_enrichment_mean = nanmean(delta_boot_array);
mean_enrichment_se = nanstd(delta_boot_array);


mean_titration = figure;
hold on
e = errorbar(mf_plot_index,mean_enrichment_mean,mean_enrichment_se,'--','Color','black','LineWidth',1);
e.CapSize = 0;
scatter(mf_plot_index,mean_enrichment_mean,'MarkerFaceColor','black','MarkerEdgeAlpha',0);
ylabel('Dl enrichment at locus (au)')
xlabel('nuclear Dl concentration (au)')
grid on 
box on
set(gca,'Fontsize',14)
saveas(mean_titration,[FigPath 'avg_titration_trend.png'])
saveas(mean_titration,[FigPath 'avg_titration_trend.pdf'])


%% Fit 3rd order polynomial
close all
nan_filter = ~isnan(delta_vec_full) & ~isnan(mf_vec_full);
p3 = polyfit(mf_vec_full(nan_filter), delta_vec_full(nan_filter), 3);
p2 = polyfit(mf_vec_full(nan_filter), delta_vec_full(nan_filter), 3);

p2_trend = polyval(p2,mf_plot_index);
p3_trend = polyval(p3,mf_plot_index);

fit_titration = figure;
hold on
pl2 = plot(mf_plot_index,p2_trend,'-','LineWidth',1.5);
e = errorbar(mf_plot_index,mean_enrichment_mean,mean_enrichment_se,'o','Color','black','LineWidth',1);
e.CapSize = 0;
s = scatter(mf_plot_index,mean_enrichment_mean,'MarkerFaceColor','black','MarkerEdgeAlpha',0);

% pl3 = plot(mf_plot_index,p3_trend,'-d','LineWidth',1.5);
legend([s pl2], 'mean','2nd order polynomial','Location','northwest')
ylabel('Dl enrichment at locus (au)')
xlabel('nuclear Dl concentration (au)')
grid on 
box on
set(gca,'Fontsize',14)
saveas(fit_titration,[FigPath 'fit_titration_trend.png'])
saveas(fit_titration,[FigPath 'fit_titration_trend.pdf'])


% fit_fig = figure;
% fit_ax = gca;
% e = errorbar(mf_index,delta_v_mf_c_tt_mean(:,end),delta_v_mf_c_tt_ste(:,end),'Color',[.6 .6 .6],'LineWidth',2);
% hold on
% e.CapSize = 0;
% p1 = plot(mf_index,enrichment_pd1,'Color',cm2(120,:),'LineWidth',1.5);
% p2 = plot(mf_index,enrichment_pd3,'Color',cm2(10,:),'LineWidth',1.5);
% p3 = plot(0,0);
% legend([e p1 p2],'data','linear','sigmoid','Location','northwest')
% xlabel(['average ' protein_name ' concentration (au)'])
% ylabel(['absolute ' protein_name ' enrichment (au)'])
% grid on
% StandardFigure([p3],gca)
% xlim([min(mf_index),max(mf_index)])
% saveas(fit_fig,[FigPath 'mf_enrichment_prediction.png'])
% 
% % make paper fig in PBoC style
% fit_ax = gca;
% StandardFigurePBoC(e,fit_ax);
% saveas(fit_fig, [paperFigPath write_string '_mf_enrichment_prediction.pdf'])
% 
% 
% pd3_delta_protein_tt = param(1)+(param(2)-param(1))./(1+10.^((param(3)-mf_protein_vec_dist)*param(4)));
% pd3_delta_vec = NaN(size(tt_index));
% for t = 1:numel(tt_index)
%     if tt_index(t) < min(time_vec_dist) || tt_index(t) > max(time_vec_dist)
%         continue
%     end
%     t_weights = exp(-.5*((time_vec_dist-tt_index(t))/tt_sigma).^2);
%     pd3_delta_vec(t) = nansum(pd3_delta_protein_tt.*t_weights)/nansum(t_weights);        
% end
% nan_filter = ~isnan(mf_protein_vec_dist) & ~isnan(delta_protein_vec_dist);
% X = [ones(sum(nan_filter),1) mf_protein_vec_dist(nan_filter)'];
% % linear model
% beta = X \ delta_protein_vec_dist(nan_filter)';
% % third order polynomial
% p = polyfit(mf_protein_vec_dist(nan_filter),delta_protein_vec_dist(nan_filter),3);
% [param]=sigm_fit(mf_protein_vec_dist(nan_filter),delta_protein_vec_dist(nan_filter));
% 
% enrichment_pd1 =  beta(1) + beta(2)*mf_index;
% enrichment_pd2 = polyval(p,mf_index);
% enrichment_pd3 = param(1)+(param(2)-param(1))./(1+10.^((param(3)-mf_index)*param(4)));