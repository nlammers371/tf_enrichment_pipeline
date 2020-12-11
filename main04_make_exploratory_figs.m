% Script conduct locus enrichment analyses
function main04_make_exploratory_figs(projectName, varargin)

%% %%%%%%%%%%%%%%%%%% Handle paths and set default options %%%%%%%%%%%%%%%% 
close all force
addpath(genpath('utilities'))

[liveProject, ~, dataName, hasAPInfo, has3DSpotInfo, hasProteinInfo, hasNucleusProbFiles] = headerFunction(projectName);

if ~hasProteinInfo
  warning('No input protein info associated with this project. Aborting protein sampling')
  return
end

DistLim = 0.8; % min distance from edge permitted (um)
NBoots = 100; % number of bootstrap samples to use for estimating SE
cMapHeat = 'viridis';%flipud(brewermap([],'Spectral'));
lbRaw = [];
ubRaw = [];
lbFold = [];
ubFold = [];
lbDiff = [];
ubDiff = [];


% make figure path
FigurePath = [liveProject.figurePath '\basic_figs\'];
mkdir(FigurePath);

% % check for optional inputs
% for i = 1:(numel(varargin)-1)  
%     if i ~= numel(varargin)
%         if ~ischar(varargin{i+1})
%             eval([varargin{i} '=varargin{i+1};']);        
%         end
%     end    
% end

%% %%%%%%%% Load  data sets to analyze and extract experiment info %%%%%%%%

% Load analysis data
load([liveProject.dataPath 'spot_struct_protein.mat'], 'spot_struct_protein');
load([liveProject.dataPath 'snip_data.mat'], 'snip_data');

% extract protein and fluorophore info
% NL: would be nice to have a way to automatically identify target
% gene/reporter
inputStringRaw = liveProject.includedExperiments{1}.Channels{liveProject.includedExperiments{1}.inputChannels(1)};
inputString = inputStringRaw(1:strfind(inputStringRaw, ':')-1);

outputStringRaw = liveProject.includedExperiments{1}.Channels{liveProject.includedExperiments{1}.spotChannels(1)};
outputString = outputStringRaw(1:strfind(outputStringRaw, ':')-1);

% Generate distance vector for filtering snip stacks
FrameInfo = getFrameInfo(liveProject.includedExperiments{1});
PixelSize = FrameInfo(1).PixelSize;

%% %%%%%%%%%%%%%%%%% extract basic sampling info %%%%%%%%%%%%%%%%%%%%%%%%%%

% generate protein vectors 
spot_protein_vec = [spot_struct_protein.spot_protein_vec];
null_protein_vec = [spot_struct_protein.edge_null_protein_vec];
dist_vec = [spot_struct_protein.spot_edge_dist_vec]*PixelSize;

% make figure to look for enrichment edge artifact 
% makeEdgeArtifactFigure()

%% %%%%%%%%%%%%%%%%%%%%% Compare protein snippets %%%%%%%%%%%%%%%%%%%%%%%%%

% Apply distance filter and make histogram figures
null_protein_vec_dist = null_protein_vec(dist_vec>=DistLim);
spot_protein_vec_dist = spot_protein_vec(dist_vec>=DistLim);
delta_protein_vec_dist = spot_protein_vec_dist-null_protein_vec_dist;

% extract time & fluo vector
time_vec_dist = [spot_struct_protein.time];
time_vec_dist = time_vec_dist(dist_vec>=DistLim);
fluo_vec_dist = [spot_struct_protein.fluo];
fluo_vec_dist = fluo_vec_dist(dist_vec>=DistLim);
nucleus_protein_vec_dist = [spot_struct_protein.nuclear_protein_vec];
nucleus_protein_vec_dist = nucleus_protein_vec_dist(dist_vec>=DistLim);


%% %%%%%%%%%%%%%%%%%%%%%%%%% Generate snip stacks %%%%%%%%%%%%%%%%%%%%%%%%%

[spot_protein_snips, spot_mcp_snips, edge_control_protein_snips, edge_control_mcp_snips] = ...
                                      generateProteinSnips(snip_data, ~isnan(null_protein_vec)&dist_vec>=DistLim);

% take average
spot_protein_snip_mean = nanmean(spot_protein_snips,3);
edge_control_protein_snip_mean = nanmean(edge_control_protein_snips,3);
edge_control_mcp_snip_mean = nanmean(edge_control_mcp_snips,3);
spot_mcp_snip_mean = nanmean(spot_mcp_snips,3);


% Plot heatmaps
visibleOn = true;
spot_protein_snip_mean_title = [inputString ' at ' 'Target Locus'];
spot_protein_snip_mean_ylabel = [inputString ' concentration (au)'];
[spot_protein_snip_heatmap, ubRaw, lbRaw] = makeHeatmapPlots(spot_protein_snip_mean, ...
    visibleOn, spot_protein_snip_mean_title, spot_protein_snip_mean_ylabel,...
    cMapHeat,PixelSize,lbRaw,ubRaw);
  
saveas(spot_protein_snip_heatmap,[FigurePath 'mean_pt_snippet_spot' '.png']);
saveas(spot_protein_snip_heatmap, [FigurePath 'mean_pt_snippet_spot' '.pdf']);

null_protein_snip_mean_title = [inputString ' at Control Locus'];
null_protein_snip_mean_ylabel = [inputString ' concentration (au)'];
null_protein_snip_heatmap = makeHeatmapPlots(edge_control_protein_snip_mean, ...
    visibleOn, null_protein_snip_mean_title, null_protein_snip_mean_ylabel,...
    cMapHeat,PixelSize,lbRaw,ubRaw);
saveas(null_protein_snip_heatmap,[FigurePath 'mean_pt_snippet_null' '.png']);
saveas(null_protein_snip_heatmap, [FigurePath 'mean_pt_snippet_null' '.pdf']);


% Make fold diff image
rel_protein_snip_mean = (spot_protein_snip_mean) ./ edge_control_protein_snip_mean;
% caxis([lb ub])
rel_protein_snip_mean_title = ['Relative ' inputString ' Enrichment at Target Locus'];
rel_protein_snip_mean_clabel = [inputString ' fold enrichment'];
rel_protein_snip_heatmap = makeHeatmapPlots(rel_protein_snip_mean, ...
    visibleOn, rel_protein_snip_mean_title, rel_protein_snip_mean_clabel, ...
    cMapHeat,PixelSize,lbFold,ubFold);
saveas(rel_protein_snip_heatmap,[FigurePath 'mean_pt_snippet_rel' '.png']);
saveas(rel_protein_snip_heatmap, [FigurePath 'mean_pt_snippet_rel' '.pdf']);

% Make absolute diff image
absDiff_protein_snip_mean = (spot_protein_snip_mean) - edge_control_protein_snip_mean;
% sumEnrichedProtein = sum(sum(absDiff_protein_snip_mean));
% disp(['Total additional protein (au) at locus (sum of all pixels of absolute different between spot and null snips)' num2str(sumEnrichedProtein)])
% caxis([lb ub])
absDiff_protein_snip_mean_title = ['Absolute Difference ' inputString ' Enrichment at Target Locus'];
absDiff_protein_snip_mean_clabel = [inputString ' absolute enrichment (au)'];
absDiff_protein_snip_heatmap = makeHeatmapPlots(absDiff_protein_snip_mean,...
    visibleOn, absDiff_protein_snip_mean_title, absDiff_protein_snip_mean_clabel,...
    cMapHeat,PixelSize,lbDiff,ubDiff);

saveas(absDiff_protein_snip_heatmap,[FigurePath 'mean_pt_snippet_absDiff' '.png']);
saveas(absDiff_protein_snip_heatmap, [FigurePath 'mean_pt_snippet_absDiff' '.pdf']);

%% %%%%%%%%%%%%%%%%%%%%%%%%%MS2 spot fluorescence %%%%%%%%%%%%%%%%%%%%%%%%%

% Plot MS2 spot heatmaps
spot_mean_fluo_title = [outputString ' Intensity at Active Locus'];
null_mean_fluo_title = [outputString ' Intensity at Control Locus'];
mean_fluo_clabel = ['MS2 intensity (au)'];
spot_mean_fluo_heatmap = makeHeatmapPlots(spot_mcp_snip_mean, visibleOn, ...
    spot_mean_fluo_title, mean_fluo_clabel,cMapHeat,PixelSize,lbRaw,ubRaw);
saveas(spot_mean_fluo_heatmap,[FigurePath 'mean_fluo_snippet_spot' '.png']);
saveas(spot_mean_fluo_heatmap, [FigurePath 'mean_fluo_snippet_spot' '.pdf']);

edge_control_mean_fluo_heatmap = makeHeatmapPlots(edge_control_mcp_snip_mean, visibleOn, ...
    null_mean_fluo_title, mean_fluo_clabel, cMapHeat,PixelSize,lbRaw,ubRaw);
saveas(edge_control_mean_fluo_heatmap,[FigurePath 'mean_fluo_snippet_null' '.png']);
saveas(edge_control_mean_fluo_heatmap, [FigurePath 'mean_fluo_snippet_null' '.pdf']);


%% %%%%%%%%%%%%%%%%%% Make Radial Profile Figure %%%%%%%%%%%%%%%%%%%%%%%%%%
% Make r reference array where each element is its distance, in um, from 
% the center pixel in the array
snip_size = size(spot_protein_snips,1);
[y_ref, x_ref] = meshgrid(1:snip_size,1:snip_size);
r_ref = sqrt((x_ref - ceil(snip_size/2)).^2 + (y_ref - ceil(snip_size/2)).^2)*PixelSize;
dist_index = unique(r_ref);
% generate smoothly interpolated axis for plot
dist_plot_axis = linspace(dist_index(1),dist_index(end),50);

% indexing vector for sampling
index_vec = 1:size(edge_control_protein_snips,3);

% define scale for moving average
r_sigma = PixelSize;
maxBootSize = 5e3;
% initialize profile arrays
r_spot_mat = NaN(NBoots,numel(dist_plot_axis));
r_control_mat = NaN(NBoots,numel(dist_plot_axis));

for n = 1:NBoots
    s_ids = randsample(index_vec,min([maxBootSize length(index_vec)]),true);
    pt_spot = nanmean(spot_protein_snips(:,:,s_ids),3);
    pt_null = nanmean(edge_control_protein_snips(:,:,s_ids),3);
    for r = 1:length(dist_plot_axis)
        r_weights = exp(-.5*((r_ref-dist_plot_axis(r))/r_sigma).^2);
        r_spot_mat(n,r) = nansum(pt_spot(:).*r_weights(:)) ./ nansum(r_weights(:));
        r_control_mat(n,r) = nansum(pt_null(:).*r_weights(:)) ./ nansum(r_weights(:));
    end
end

r_spot_mean = nanmean(r_spot_mat);
r_spot_ste = nanstd(r_spot_mat);

r_control_mean = nanmean(r_control_mat);
r_control_ste = nanstd(r_control_mat);

% % save table with radial profile results
% radial_table = array2table([dist_index',r_spot_mean',r_spot_ste',r_control_mean',r_control_ste'],...
%             'VariableNames',{'radial_distance','mean_locus_protein','ste_locus_protein','mean_control_protein','ste_control_protein'});
% writetable(radial_table,[DataPath 'radial_profile_data.csv'])

% make figure
cm2 = brewermap([],'Set2');

r_fig = figure;
r_ax = gca;
hold on
e = errorbar(dist_plot_axis,r_control_mean / r_control_mean(1),r_control_ste / r_control_mean(1),'Color','black','LineWidth',1.75);
e.CapSize = 0;
e = errorbar(dist_plot_axis,r_spot_mean / r_control_mean(1),r_spot_ste / r_control_mean(1),'Color',cm2(3,:),'LineWidth',1.75);
e.CapSize = 0;
grid off
r_ax.XLabel.String = 'radius (\mu m)';
r_ax.YLabel.String ='relative enrichment';
legend('control','active locus')
title(['Radial Concentration Profile (' inputString ')'])
r_ax.XLim = [0 1.2];
StandardFigure(e,r_ax);
% r_ax.YLim = [relEnrich_lb relEnrich_ub];
saveas(r_fig, [FigurePath '_radial_enrichment.png'])


% %%%%%%% plot relative enrichment at locus as a function of time %%%%%%%%%%%
% 
% % take bootstrap samples
% tt_index = 0:60:round(max(time_vec_dist));
% delta_protein_time_mat = NaN(numel(tt_index),NBoots);
% fluo_time_mat = NaN(numel(tt_index),NBoots);
% spot_protein_time_mat = NaN(numel(tt_index),NBoots);
% null_protein_time_mat = NaN(numel(tt_index),NBoots);
% tt_sigma = 60; % seconds
% for n = 1:NBoots
%     s_ids = randsample(1:numel(delta_protein_vec_dist),numel(delta_protein_vec_dist),true);
%     dv_samp1 = delta_protein_vec_dist(s_ids);    
%     nn_samp = null_protein_vec_dist(s_ids);
%     ff_samp = fluo_vec_dist(s_ids);
%     sp_samp = spot_protein_vec_dist(s_ids);
%     time_samp = time_vec_dist(s_ids);
%     for t = 1:numel(tt_index)
%         if tt_index(t) < min(time_vec_dist) || tt_index(t) > max(time_vec_dist)
%             continue
%         end
%         t_weights = exp(-.5*((time_samp-tt_index(t))/tt_sigma).^2);
%         delta_protein_time_mat(t,n) = nansum(dv_samp1.*t_weights)/nansum(t_weights);        
%         null_protein_time_mat(t,n) =  nansum(nn_samp.*t_weights)/nansum(t_weights);
%         spot_protein_time_mat(t,n) = nansum(sp_samp.*t_weights)/nansum(t_weights);
%         fluo_time_mat(t,n) = nansum(ff_samp.*t_weights)/nansum(t_weights);
%     end
% end
% % calculate average and standard error
% time_delta_mean = nanmean(delta_protein_time_mat,2);
% time_delta_ste = nanstd(delta_protein_time_mat,[],2);
% time_null_mean = nanmean(null_protein_time_mat,2);
% time_null_ste = nanstd(null_protein_time_mat,[],2);
% time_spot_mean = nanmean(spot_protein_time_mat,2);
% time_spot_ste = nanstd(spot_protein_time_mat,[],2);
% time_fluo_mean = nanmean(fluo_time_mat,2);
% time_fluo_ste = nanstd(fluo_time_mat,[],2);
% cm1 = brewermap([],'Set3');
% 
% delta_time_fig = figure;
% yyaxis left
% e = errorbar(tt_index/60,time_delta_mean,time_delta_ste,'Color',cm1(5,:),'LineWidth',2);
% e.CapSize = 0;
% ylabel([protein_name ' enrichment at locus (au)'])
% grid on
% y_max = nanmax(time_delta_mean);
% y_min = nanmin(time_delta_mean);
% ylim([y_min-.1*y_min*sign(y_min) 1.1*y_max])
% 
% yyaxis right
% plot(tt_index/60, time_fluo_mean,'LineWidth',2,'Color',cm1(9,:));
% e.CapSize = 0;
% hold on
% p = plot(0,0);
% ax = gca;
% ax.YColor = [.2 .2 .2];
% xlim([tt_index(1) tt_index(end)]/60)
% ylabel([gene_name ' activity (au)'])
% legend('enrichment trend','activity trend', 'Location','northeast')
% xlabel('minutes')
% % title(['Enrichment of ' id_string])
% StandardFigure(p,gca);
% saveas(delta_time_fig, [FigurePath '_temporal_enrichment_w_mcp.png'])
% 
% % make individual plots
% delta_time_fig = figure;
% e = errorbar(tt_index/60,time_delta_mean,time_delta_ste,'Color',cm1(5,:),'LineWidth',1.5);
% e.CapSize = 0;
% ylabel([protein_name ' enrichment at locus (au)'])
% grid on
% StandardFigure(e,gca);
% xlabel('minutes')
% xlim([tt_index(1) tt_index(end)]/60)
% saveas(delta_time_fig, [FigurePath '_temporal_enrichment.png'])
% 
% 
% null_time_fig = figure;
% yyaxis left
% e = errorbar(tt_index/60,time_null_mean,time_null_ste,'Color',cm1(4,:),'LineWidth',2);
% e.CapSize = 0;
% ylabel(['average ' protein_name ' conc.(au)'])
% ax = gca;
% ax.YColor = cm1(4,:);
% 
% yyaxis right
% plot(tt_index/60, time_fluo_mean,'LineWidth',2,'Color',cm1(9,:));
% hold on
% p = plot(0,0);
% ax = gca;
% ax.YColor = cm1(9,:);
% xlim([tt_index(1) tt_index(end)]/60)
% ylabel([gene_name ' activity (au)'])
% 
% legend('background trend','activity trend', 'Location','northeast')
% StandardFigure(p,gca);
% xlabel('minutes')
% xlim([6 60])
% saveas(null_time_fig, [FigurePath '_temporal_background_w_fluo.png'])
% 
% %%
% %%%%%%%%%%% Is trend a function of time or protein concentration? %%%%%%%%%
% prctile_vec = [0 20 40 60 80 100];
% mf_prctile_vec = NaN(size(prctile_vec));
% tt_prctile_vec =  60*(0:10:50);
% for p = 1:numel(prctile_vec)
%     mf_prctile_vec(p) = prctile(nucleus_protein_vec_dist,prctile_vec(p));
% end
% mf_id_vec = NaN(size(nucleus_protein_vec_dist));
% tt_id_vec = NaN(size(time_vec_dist));
% for p = 1:numel(mf_prctile_vec)-1
%     mf_ft = nucleus_protein_vec_dist >= mf_prctile_vec(p) & nucleus_protein_vec_dist < mf_prctile_vec(p+1);
%     mf_id_vec(mf_ft) = p;
%     tt_ft = time_vec_dist >= tt_prctile_vec(p) & time_vec_dist < tt_prctile_vec(p+1);
%     tt_id_vec(tt_ft) = p;
% end
% mf_vec_dist = nucleus_protein_vec_dist;
% tt_vec_dist = time_vec_dist;
% mf_index = linspace(prctile(nucleus_protein_vec_dist,1),prctile(nucleus_protein_vec_dist,99),50);
% % track enrichemnt trends as function of space and protein concentration
% delta_v_tt_c_mf_array = NaN(numel(tt_index),numel(prctile_vec)-1,NBoots);
% delta_v_mf_c_tt_array = NaN(numel(mf_index),numel(prctile_vec)-1,NBoots);
% 
% mf_sigma = median(diff(mf_index));
% 
% dep_var_cell = {'mf','tt'};
% dynamic_ids = [];
% boot_sigma = [];
% dynamic_var_vec = [];
% index = [];
% for i = 1:numel(dep_var_cell)
%     dynamic_var = dep_var_cell{i};
%     static_var = dep_var_cell{[1,2]~=i};
%     eval(['var_array = delta_v_' dynamic_var '_c_' static_var '_array;'])
%     for j = 1:numel(prctile_vec)
%         if j < numel(prctile_vec)
%             eval(['dynamic_ids = find(' static_var '_id_vec == j);'])
%         else
%             eval(['dynamic_ids = 1:numel(' static_var '_id_vec);'])
%         end
%         eval(['boot_sigma = ' dynamic_var '_sigma;'])
%         eval(['dynamic_var_vec = ' dynamic_var '_vec_dist;'])
%         eval(['index = ' dynamic_var '_index;'])
%         if numel(dynamic_ids) < 100
%             continue
%         end
%         for n = 1:NBoots
%             s_ids = randsample(dynamic_ids,numel(dynamic_ids),true);
%             delta_boot = delta_protein_vec_dist(s_ids);
%             dynamic_boot = dynamic_var_vec(s_ids);
%             for k = 1:numel(index)
%                 weights = exp(-.5*((dynamic_boot-index(k))/boot_sigma).^2);
%                 var_array(k,j,n) = nansum(delta_boot.*weights)/nansum(weights); 
%             end
%         end
%     end
%     eval(['delta_v_' dynamic_var '_c_' static_var 'mean = nanmean(var_array,3);'])
%     eval(['delta_v_' dynamic_var '_c_' static_var '_ste = nanstd(var_array,[],3);'])
%     eval(['delta_v_' dynamic_var '_c_' static_var '_array = var_array;'])
%     clear var_array
% end    
% 
% % make protein cohort figures
% cm2 = brewermap(128,'Spectral');
% color_array = cm2(1:30:128,:);
% constant_mf_fig = figure;
% hold on
% lgd_str = {};
% e = [];
% for i = 1:numel(prctile_vec)-1
%     e = [e errorbar(tt_index/60, delta_v_tt_c_mf_mean(:,i), delta_v_tt_c_mf_ste(:,i),'Color',color_array(i,:),'CapSize',0)];
%     lgd_str = [lgd_str{:} {['protein cohort ' num2str(i)]}];
% end
% legend(e,lgd_str{:})
% xlabel('time (minutes)')
% ylabel(['absolute ' protein_name 'enrichment (au)'])
% saveas(constant_mf_fig,[FigurePath 'protein_cohort_plot.png'])
% 
% constant_tt_fig = figure;
% hold on
% lgd_str = {};
% e = [];
% for i = 1:numel(prctile_vec)-1
%     e = [e errorbar(mf_index, delta_v_mf_c_tt_mean(:,i), delta_v_mf_c_tt_ste(:,i),'Color',color_array(i,:),'CapSize',0)];
%     lgd_str = [lgd_str{:} {['time cohort ' num2str(i)]}];
% end
% legend(e,lgd_str{:})
% xlabel(['average ' protein_name ' concentration (au)'])
% ylabel(['absolute ' protein_name 'enrichment (au)'])
% saveas(constant_tt_fig,[FigurePath 'time_cohort_plot.png'])
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Test hunch that there's a simple linear relationship btw mf and
% % enrichment
% nan_filter = ~isnan(nucleus_protein_vec_dist) & ~isnan(delta_protein_vec_dist);
% X = [ones(sum(nan_filter),1) nucleus_protein_vec_dist(nan_filter)'];
% % linear model
% beta = X \ delta_protein_vec_dist(nan_filter)';
% % third order polynomial
% p = polyfit(nucleus_protein_vec_dist(nan_filter),delta_protein_vec_dist(nan_filter),3);
% [param]=sigm_fit(nucleus_protein_vec_dist(nan_filter),delta_protein_vec_dist(nan_filter));
% 
% enrichment_pd1 =  beta(1) + beta(2)*mf_index;
% enrichment_pd2 = polyval(p,mf_index);
% enrichment_pd3 = param(1)+(param(2)-param(1))./(1+10.^((param(3)-mf_index)*param(4)));
% 
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
% saveas(fit_fig,[FigurePath 'mf_enrichment_prediction.png'])
% 
% % make paper fig in PBoC style
% fit_ax = gca;
% StandardFigurePBoC(e,fit_ax);
% saveas(fit_fig, [FigurePath '_mf_enrichment_prediction.pdf'])
% 
% 
% pd3_delta_protein_tt = param(1)+(param(2)-param(1))./(1+10.^((param(3)-nucleus_protein_vec_dist)*param(4)));
% pd3_delta_vec = NaN(size(tt_index));
% for t = 1:numel(tt_index)
%     if tt_index(t) < min(time_vec_dist) || tt_index(t) > max(time_vec_dist)
%         continue
%     end
%     t_weights = exp(-.5*((time_vec_dist-tt_index(t))/tt_sigma).^2);
%     pd3_delta_vec(t) = nansum(pd3_delta_protein_tt.*t_weights)/nansum(t_weights);        
% end
% 
% % remake figure
% delta_time_fig = figure;
% hold on
% e = errorbar(tt_index/60,time_delta_mean,time_delta_ste,'LineWidth',1.5,'Color',[.6 .6 .6]);
% e.CapSize = 0;
% s = plot(tt_index/60,pd3_delta_vec,'-','LineWidth',1.5,'Color',cm2(10,:));
% ylabel([protein_name ' enrichment at locus (au)'])
% grid on
% xlim([6 60])
% xlabel('minutes')
% % title(['Enrichment of ' id_string])
% box on
% p = plot(0,0);
% legend([e s],'enrichment trend (actual)','predicted enrichment (sigmoid)',...
%     'Location','northeast')
% StandardFigure(p,gca)
% saveas(delta_time_fig, [FigurePath '_temporal_enrichment_prediction.png'])
% 
% 

% 
% 
% 
% % % take bootstrap samples
% % ap_index = unique(round(snip_ap_vec_dist));
% % delta_protein_ap_mat = NaN(numel(ap_index),NBoots);
% % spot_protein_ap_mat = NaN(numel(ap_index),NBoots);
% % null_protein_ap_mat = NaN(numel(ap_index),NBoots);
% % ap_sigma= 2; % seconds
% % for n = 1:NBoots
% %     s_ids = randsample(1:numel(delta_protein_vec_dist),numel(delta_protein_vec_dist),true);
% %     dv_samp = delta_protein_vec_dist(s_ids);    
% %     nn_samp = null_protein_vec_dist(s_ids);
% %     sp_samp = spot_protein_vec_dist(s_ids);
% %     ap_samp = snip_ap_vec_dist(s_ids);
% %     for a = 1:numel(ap_index)        
% %         ap_weights = exp(-.5*((ap_samp-ap_index(a))/ap_sigma).^2);
% %         delta_protein_ap_mat(a,n) = nansum(dv_samp.*ap_weights)/nansum(ap_weights);        
% %         null_protein_ap_mat(a,n) =  nansum(nn_samp.*ap_weights)/nansum(ap_weights);
% %         spot_protein_ap_mat(a,n) = nansum(sp_samp.*ap_weights)/nansum(ap_weights);
% %     end
% % end
% % % calculate average and standard error
% % ap_delta_mean = nanmean(delta_protein_ap_mat./null_protein_ap_mat,2);
% % ap_delta_ste = nanstd(delta_protein_ap_mat./null_protein_ap_mat,[],2);
% % ap_null_mean = nanmean(null_protein_ap_mat,2);
% % ap_null_ste = nanstd(null_protein_ap_mat,[],2);
% % ap_spot_mean = nanmean(spot_protein_ap_mat,2);
% % ap_spot_ste = nanstd(spot_protein_ap_mat,[],2);
% % 
% % cm = Colormap_plot;
% % delta_ap_fig = figure;
% % 
% % yyaxis left
% % e = errorbar(ap_index,100*ap_delta_mean,100 * ap_delta_ste,'-','Color',cm(60,:)/1.2,'LineWidth',1.5);
% % e.CapSize = 0;
% % ylabel('% enrichment at locus')
% % y_max = nanmax(100 * ap_delta_mean);
% % y_min = nanmin(100 * ap_delta_mean);
% % ylim([y_min-.1*y_min*sign(y_min) 1.1*y_max])
% % 
% % yyaxis right
% % plot(ap_index, ap_spot_mean,'-','Color',cm(30,:));
% % hold on
% % plot(ap_index, ap_null_mean,'-','Color','black');
% % ylabel('concentration (au)')
% % 
% % legend('percent enrichment','locus concentration', 'control concentration','Location','northwest')
% % 
% % title(['Percent Enrichment of ' id_string])
% % xlabel('% AP')
% % grid on
% % saveas(delta_ap_fig, [FigurePath '_ap_percent_enrichment.png'])
% end