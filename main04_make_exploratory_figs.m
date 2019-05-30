% Script conduct locus enrichment analyses
function main04_make_exploratory_figs(project,protein_string,gene_string,varargin)

close all

dataPath = ['../../dat/' project '/'];
figPath = ['../../fig/' project '/'];
DistLim = .8; % min distance from edge permitted (um)
NBoots = 20;%00; % number of bootstrap samples to use for estimating SE
ManualDistThreshold = 0;
Colormap_plot = jet(128); %specifies the colomap used to make plots/graphs
Colormap_heat = viridis(128); %specifies the colormap used to make heatmaps
relEnrich_ub = 1.3; %upper bound of relative enrichment for consistency
relEnrich_lb = 0.85; %lower bound of relative enrichment for consistency
dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
dataPath = [dropboxFolder 'ProcessedEnrichmentData\' project '/'];

for i = 1:numel(varargin)
    if strcmpi(varargin{i}, 'dropboxFolder')        
        dataPath = [varargin{i+1} '/ProcessedEnrichmentData/' project '/'];
        figPath = [varargin{i+1} '/LocalEnrichmentFigures/' project '/'];        
    elseif ischar(varargin{i})
        if ismember(varargin{i}, {'ControlType','ROIRadiusSpot','DistLim','NBoots','ManualDistThreshold'}) 
            eval([varargin{i} '=varargin{i+1};']);
        end
    end
end

paperFigPath = [figPath 'paperFigs/'];
mkdir(paperFigPath);

figPath = [figPath 'basicFigs/'];
mkdir(figPath)



% Load analysis data
load([dataPath 'nucleus_struct_protein.mat'], 'nucleus_struct_protein');

% extract protein, gene, fluorophore info
pt_dash = strfind(protein_string,'-');
protein_name = protein_string(1:pt_dash-1);
protein_fluor = protein_string(pt_dash+1:end);
gene_dash = strfind(gene_string,'-');
gene_name = gene_string(1:gene_dash-1);
gene_fluor = gene_string(gene_dash+1:end);

id_string = [protein_name '-' protein_fluor ' : ' gene_name '-' gene_fluor]; 
write_string = [protein_name '-' protein_fluor '__' gene_name '-' gene_fluor]; 
% Generate distance vector for filtering snip stacks
PixelSize = nucleus_struct_protein(1).PixelSize;

% Snip stacks
spot_protein_snips = cat(3,nucleus_struct_protein.spot_protein_snips);
null_protein_snips = cat(3,nucleus_struct_protein.edge_null_protein_snips);
spot_mcp_snips = cat(3,nucleus_struct_protein.spot_mcp_snips);
null_mcp_snips = cat(3,nucleus_struct_protein.edge_null_mcp_snips);

% Make r reference array
snip_size = size(spot_protein_snips,1);
[y_ref, x_ref] = meshgrid(1:snip_size,1:snip_size);
r_ref = sqrt((x_ref-ceil(snip_size/2)).^2 + (y_ref-ceil(snip_size/2)).^2)*PixelSize;

spot_protein_vec = [nucleus_struct_protein.spot_protein_vec];
null_protein_vec = [nucleus_struct_protein.edge_null_protein_vec];
dist_vec = [nucleus_struct_protein.spot_edge_dist_vec]*PixelSize;

% First look for presence of edge artifact
dist_sigma = .1; %(um)
protein_delta_vec =  spot_protein_vec - null_protein_vec;
dist_index = 0:.1:floor(prctile(dist_vec,99)*10)/10;

delta_dist_mat = NaN(numel(dist_index),NBoots);
null_dist_mat = NaN(numel(dist_index),NBoots);
spot_dist_mat = NaN(numel(dist_index),NBoots);


for n = 1:NBoots
    s_ids = randsample(1:numel(protein_delta_vec),numel(protein_delta_vec),true);
    dv_samp1 = protein_delta_vec(s_ids);    
    nn_samp = spot_protein_vec(s_ids);
    sp_samp = null_protein_vec(s_ids);
    dist_samp = dist_vec(s_ids);        
    for t = 1:numel(dist_index)
        d_weights = exp(-.5*((dist_samp-dist_index(t))/dist_sigma).^2);
        delta_dist_mat(t,n) = nansum(dv_samp1.*d_weights) / nansum(d_weights);
        null_dist_mat(t,n) = (nansum(nn_samp.*d_weights) / nansum(d_weights));
        spot_dist_mat(t,n) = (nansum(sp_samp.*d_weights) / nansum(d_weights));
    end
end
% calcualte mean and standard error
delta_dist_mean = nanmean(delta_dist_mat,2);
delta_dist_ste = nanstd(delta_dist_mat,[],2);

null_dist_mean = nanmean(null_dist_mat,2);

% make dist-dependent fold enrichment figure and select distance threshold
pass = 0;
while ~pass
    delta_dist_fig = figure;
    e = errorbar(dist_index,100* delta_dist_mean ./ null_dist_mean,100* delta_dist_ste ./ null_dist_mean);
    e.CapSize = 0;
    hold on
    xlabel('distance from edge (\mu m)')
    ylabel('apparent % enrichment')
    
    grid on
    y_max = 100 * nanmax(delta_dist_mean ./ null_dist_mean);
    y_min = 100 * nanmin(delta_dist_mean ./ null_dist_mean);
    p = plot([DistLim DistLim],[y_min y_max],'Color', 'red');
    legend(p, 'current limit')
    ylim([min([0 , y_min-.1*y_min*sign(y_min)]) 1.1*y_max])  
    if ManualDistThreshold
        title(strvcat('Setting Edge Distance Threshold:',...
                'If current value is staisfactory press "Enter"',...
                'Else click a desired cutoff and click "Enter"'))
        [x,~] = ginput;
        if isempty(x)
            pass = 1;
        end
    else
        pass = 1;
    end
end
title(['Enrichment vs. Distance from Nucleus Edge (' id_string ')'])
saveas(delta_dist_fig, [figPath write_string '_edge_artifact_plot.png'])

% Apply distance filter and make histogram figures
null_protein_vec_dist = null_protein_vec(dist_vec>=DistLim);
spot_protein_vec_dist = spot_protein_vec(dist_vec>=DistLim);
delta_protein_vec_dist = spot_protein_vec_dist-null_protein_vec_dist;
% extract time & fluo vector
time_vec_dist = [nucleus_struct_protein.time];
time_vec_dist = time_vec_dist(dist_vec>=DistLim);
fluo_vec_dist = [nucleus_struct_protein.fluo];
fluo_vec_dist = fluo_vec_dist(dist_vec>=DistLim);
mf_protein_vec_dist = [nucleus_struct_protein.mf_null_protein_vec];
mf_protein_vec_dist = mf_protein_vec_dist(dist_vec>=DistLim);


% delta_norm_protein_vec_dist = (spot_protein_vec_dist-null_protein_vec_dist) / nanmean(null_protein_vec_dist);
% 
% n_unique = numel(unique([null_protein_vec_dist spot_protein_vec_dist]));
% % plot raw distributions
% lb = .9*prctile([null_protein_vec_dist spot_protein_vec_dist],1);
% ub = 1.1*prctile([null_protein_vec_dist spot_protein_vec_dist],99);
% bins = linspace(lb,ub,min([ceil(n_unique/20),100]));
% 
% spot_ct = histc(spot_protein_vec_dist,bins);
% null_ct = histc(null_protein_vec_dist,bins);
% spot_ct = spot_ct / sum(spot_ct);
% null_ct = null_ct / sum(null_ct);
% 
% raw_hist_fig = figure;
% cm = Colormap_plot;
% hold on
% % hist plots
% yyaxis left
% b1 = bar(bins,null_ct,1,'FaceAlpha',.5,'FaceColor',cm(30,:));
% b2 = bar(bins,spot_ct,1,'FaceAlpha',.5,'FaceColor',cm(100,:));
% ylabel('share')
% ax = gca;
% ax.YColor = 'black';
% % cumulative plots
% yyaxis right
% plot(bins,cumsum(null_ct),'-','LineWidth',1.5,'Color',cm(30,:));
% plot(bins,cumsum(spot_ct),'-','LineWidth',1.5,'Color',cm(100,:));
% ylabel('cumulative share')
% ax = gca;
% ax.YColor = 'black';
% 
% legend('control',['active locus (' gene_name ')'])
% xlabel('concentration (au)')
% title([protein_name '-' protein_fluor ' Concentrations'])
% grid on
% saveas(raw_hist_fig,[figPath 'hist_plots_' write_string '_dist_' num2str(DistLim) '.png']); 
% 
% % plot delta distribution
% pd = fitdist(delta_norm_protein_vec_dist','Normal');
% lb = .9*prctile(delta_norm_protein_vec_dist,1);
% ub = 1.1*prctile(delta_norm_protein_vec_dist,99);
% bins = linspace(lb,ub,min([n_unique/20,100]));
% 
% delta_ct = histc(delta_norm_protein_vec_dist,bins);
% delta_ct = delta_ct / sum(delta_ct);
%    
% 
% delta_hist_fig = figure;
% cm = Colormap_plot;
% hold on
% % hist plots
% b = bar(bins,delta_ct,1,'FaceAlpha',.5,'FaceColor',cm(60,:)/1.2);
% norm_vec = exp(-.5*((bins-pd.mu)/pd.sigma).^2);
% p = plot(bins,norm_vec/sum(norm_vec),'LineWidth',2);
% ylabel('share')
% legend([b,p],['\Delta F ('  gene_name ')'], ...
%     ['Gaussian Fit (\mu=' num2str(pd.mu) ' \sigma=' num2str(pd.sigma)],'Location','best')
% xlabel('concentration (au)')
% title(['Enrichment at Active Locus Relative to Control (' id_string ')'])
% grid on
% saveas(delta_hist_fig,[figPath 'delta_hist_' write_string '_dist_' num2str(DistLim) '.png']); 


%%%%%%%%%%%%%%%%%%%%%%%% Compare protein snippets %%%%%%%%%%%%%%%%%%%%%%%%%
% apply distance filter
spot_protein_snips_dist = spot_protein_snips(:,:,dist_vec>=DistLim);
null_protein_snips_dist = null_protein_snips(:,:,dist_vec>=DistLim);
spot_mcp_snips_dist = spot_mcp_snips(:,:,dist_vec>=DistLim);
null_mcp_snips_dist = null_mcp_snips(:,:,dist_vec>=DistLim);

% randomize snip orientation....
inv_mat = [fliplr(1:snip_size); 1:snip_size]' ;
spot_protein_snips_mixed = NaN(size(spot_protein_snips_dist));
null_protein_snips_mixed = NaN(size(spot_protein_snips_dist));
spot_mcp_snips_mixed = NaN(size(spot_protein_snips_dist));
null_mcp_snips_mixed = NaN(size(spot_protein_snips_dist));
for i = 1:size(spot_protein_snips_dist,3)
    h = ceil(rand()*2);
    v = ceil(rand()*2);
    
    spot_protein_snips_mixed(:,:,i) = spot_protein_snips_dist(inv_mat(:,v),inv_mat(:,h),i);
    null_protein_snips_mixed(:,:,i) = null_protein_snips_dist(inv_mat(:,v),inv_mat(:,h),i);
    spot_mcp_snips_mixed(:,:,i) = spot_mcp_snips_dist(inv_mat(:,v),inv_mat(:,h),i);
    null_mcp_snips_mixed(:,:,i) = null_mcp_snips_dist(inv_mat(:,v),inv_mat(:,h),i);
end

spot_protein_snip_mean = nanmean(spot_protein_snips_mixed,3);
null_protein_snip_mean = nanmean(null_protein_snips_mixed,3);

lb = prctile([spot_protein_snip_mean(:)'  null_protein_snip_mean(:)'],1);
ub = prctile([spot_protein_snip_mean(:)'  null_protein_snip_mean(:)'],99);

xtick_string = "set(gca,'xtick',1:5:snip_size,'xticklabel',round(((1:5:snip_size)-round(snip_size/2))*PixelSize,2))";
ytick_string = "set(gca,'ytick',1:5:snip_size,'yticklabel',round(((1:5:snip_size)-round(snip_size/2))*PixelSize,2))";

% Plot heatmaps
spot_protein_snip_mean_title = [protein_name '-' protein_fluor ' at ' 'Active ' gene_name ' Locus'];
spot_protein_snip_mean_ylabel = [protein_name '-' protein_fluor ' concentration (au)'];
makeHeatmapPlots(spot_protein_snip_mean, spot_protein_snip_mean_title, spot_protein_snip_mean_ylabel, '_mean_pt_snippet_spot')

null_protein_snip_mean_title = [protein_name '-' protein_fluor ' at Control Locus'];
null_protein_snip_mean_ylabel = [protein_name '-' protein_fluor ' concentration (au)'];
makeHeatmapPlots(null_protein_snip_mean, null_protein_snip_mean_title, null_protein_snip_mean_ylabel, '_mean_pt_snippet_null')


% Make diff image
rel_protein_snip_mean = (spot_protein_snip_mean) ./ null_protein_snip_mean;
% caxis([lb ub])
rel_protein_snip_mean_title = ['Relative ' protein_name '-' protein_fluor ' Enrichment at Active ' gene_name ' Locus'];
rel_protein_snip_mean_ylabel = [protein_name '-' protein_fluor ' fold enrichment'];
makeHeatmapPlots(rel_protein_snip_mean, rel_protein_snip_mean_title, rel_protein_snip_mean_ylabel, '_mean_pt_snippet_rel')


% Compare mcp snippets

null_mean_fluo = nanmean(null_mcp_snips_mixed,3);
spot_mean_fluo = nanmean(spot_mcp_snips_mixed,3);
ub = prctile([reshape(null_mean_fluo,1,[]) reshape(spot_mean_fluo,1,[])],99);
lb = prctile([reshape(null_mean_fluo,1,[]) reshape(spot_mean_fluo,1,[])],1);
% Plot heatmaps
spot_mean_fluo_title = [gene_fluor ' Intensity at Active Locus (' gene_name ')'];
null_mean_fluo_title = [gene_fluor ' Intensity at Control Locus'];
mean_fluo_ylabel = [gene_name ' expression (au)'];
makeHeatmapPlots(spot_mean_fluo, spot_mean_fluo_title, mean_fluo_ylabel, '_mean_fluo_snippet_spot')
makeHeatmapPlots(null_mean_fluo, null_mean_fluo_title, mean_fluo_ylabel, '_mean_fluo_snippet_null')


%%%%%%% Nested function to make basic and PBoC-style heatmap plots %%%%%%%
function makeHeatmapPlots(Data, Title, YLabel, FileName)
    snippet_fig = figure;
    colormap(Colormap_heat)
    snippet_im = imagesc(Data);
    title(Title)
    if contains(FileName, 'rel')
        caxis([relEnrich_lb relEnrich_ub])
    else
        caxis([lb ub])
    end
    h = colorbar;
    ylabel('\mum','FontSize',12)
    xlabel('\mum','FontSize',12)
    ylabel(h,YLabel,'FontSize',12)
    eval(xtick_string)
    eval(ytick_string)
    saveas(snippet_fig,[figPath write_string FileName '.png']);

    % make paper fig in PBoC style
    snippet_ax = gca;
    StandardFigurePBoC(snippet_im,snippet_ax);
    saveas(snippet_fig, [paperFigPath write_string FileName '.pdf'])
end


%%%%%%% plot relative enrichment at locus as a function of time %%%%%%%%%%%

% take bootstrap samples
tt_index = 0:120:round(max(time_vec_dist));
delta_protein_time_mat = NaN(numel(tt_index),NBoots);
fluo_time_mat = NaN(numel(tt_index),NBoots);
spot_protein_time_mat = NaN(numel(tt_index),NBoots);
null_protein_time_mat = NaN(numel(tt_index),NBoots);
tt_sigma = 120; % seconds
for n = 1:NBoots
    s_ids = randsample(1:numel(delta_protein_vec_dist),numel(delta_protein_vec_dist),true);
    dv_samp1 = delta_protein_vec_dist(s_ids);    
    nn_samp = null_protein_vec_dist(s_ids);
    ff_samp = fluo_vec_dist(s_ids);
    sp_samp = spot_protein_vec_dist(s_ids);
    time_samp = time_vec_dist(s_ids);
    for t = 1:numel(tt_index)
        if tt_index(t) < min(time_vec_dist) || tt_index(t) > max(time_vec_dist)
            continue
        end
        t_weights = exp(-.5*((time_samp-tt_index(t))/tt_sigma).^2);
        delta_protein_time_mat(t,n) = nansum(dv_samp1.*t_weights)/nansum(t_weights);        
        null_protein_time_mat(t,n) =  nansum(nn_samp.*t_weights)/nansum(t_weights);
        spot_protein_time_mat(t,n) = nansum(sp_samp.*t_weights)/nansum(t_weights);
        fluo_time_mat(t,n) = nansum(ff_samp.*t_weights)/nansum(t_weights);
    end
end
% calculate average and standard error
time_delta_mean = nanmean(delta_protein_time_mat,2);
time_delta_ste = nanstd(delta_protein_time_mat,[],2);
% time_null_mean = nanmean(null_protein_time_mat,2);
% time_null_ste = nanstd(null_protein_time_mat,[],2);
% time_spot_mean = nanmean(spot_protein_time_mat,2);
% time_spot_ste = nanstd(spot_protein_time_mat,[],2);
time_fluo_mean = nanmean(fluo_time_mat,2);
time_fluo_ste = nanstd(fluo_time_mat,[],2);

delta_time_fig = figure;
yyaxis left
e = errorbar(tt_index/60,time_delta_mean,time_delta_ste,'LIneWidth',1.5);
e.CapSize = 0;

ylabel([protein_name ' enrichment at locus (au)'])
grid on
y_max = nanmax(time_delta_mean);
y_min = nanmin(time_delta_mean);
ylim([y_min-.1*y_min*sign(y_min) 1.1*y_max])

yyaxis right
e = errorbar(tt_index/60, time_fluo_mean,time_fluo_ste,'LineWidth',1,'Color',[.2 .2 .2]);
e.CapSize = 0;
ax = gca;
ax.YColor = [.2 .2 .2];
ylim([100 240])
ylabel([gene_name ' activity (au)'])

legend('enrichment trend','activity trend', 'Location','northeast')
xlabel('minutes')
title(['Enrichment of ' id_string])

saveas(delta_time_fig, [figPath write_string '_temporal_enrichment.png'])

%%%%%%%%%%% Is trend a function of time or protein concentration? %%%%%%%%%
prctile_vec = [0 20 40 60 80 100];
mf_prctile_vec = NaN(size(prctile_vec));
tt_prctile_vec =  60*(0:10:50);
for p = 1:numel(prctile_vec)
    mf_prctile_vec(p) = prctile(mf_protein_vec_dist,prctile_vec(p));
end
mf_id_vec = NaN(size(mf_protein_vec_dist));
tt_id_vec = NaN(size(time_vec_dist));
for p = 1:numel(mf_prctile_vec)-1
    mf_ft = mf_protein_vec_dist >= mf_prctile_vec(p) & mf_protein_vec_dist < mf_prctile_vec(p+1);
    mf_id_vec(mf_ft) = p;
    tt_ft = time_vec_dist >= tt_prctile_vec(p) & time_vec_dist < tt_prctile_vec(p+1);
    tt_id_vec(tt_ft) = p;
end
mf_vec_dist = mf_protein_vec_dist;
tt_vec_dist = time_vec_dist;
mf_index = linspace(min(mf_protein_vec_dist),prctile(mf_protein_vec_dist,99),50);
% track enrichemnt trends as function of space and protein concentration
delta_v_tt_c_mf_array = NaN(numel(tt_index),numel(prctile_vec)-1,NBoots);
delta_v_mf_c_tt_array = NaN(numel(mf_index),numel(prctile_vec)-1,NBoots);

mf_sigma = 2*median(diff(mf_index));

dep_var_cell = {'mf','tt'};
dynamic_ids = [];
boot_sigma = [];
dynamic_var_vec = [];
index = [];
for i = 1:numel(dep_var_cell)
    dynamic_var = dep_var_cell{i};
    static_var = dep_var_cell{[1,2]~=i};
    eval(['var_array = delta_v_' dynamic_var '_c_' static_var '_array;'])
    for j = 1:numel(prctile_vec)
        if j < numel(prctile_vec)
            eval(['dynamic_ids = find(' static_var '_id_vec == j);'])
        else
            eval(['dynamic_ids = 1:numel(' static_var '_id_vec);'])
        end
        eval(['boot_sigma = ' dynamic_var '_sigma;'])
        eval(['dynamic_var_vec = ' dynamic_var '_vec_dist;'])
        eval(['index = ' dynamic_var '_index;'])
        if numel(dynamic_ids) < 100
            continue
        end
        for n = 1:NBoots
            s_ids = randsample(dynamic_ids,numel(dynamic_ids),true);
            delta_boot = delta_protein_vec_dist(s_ids);
            dynamic_boot = dynamic_var_vec(s_ids);
            for k = 1:numel(index)
                weights = exp(-.5*((dynamic_boot-index(k))/boot_sigma).^2);
                var_array(k,j,n) = nansum(delta_boot.*weights)/nansum(weights); 
            end
        end
    end
    eval(['delta_v_' dynamic_var '_c_' static_var '_mean = nanmean(var_array,3);'])
    eval(['delta_v_' dynamic_var '_c_' static_var '_ste = nanstd(var_array,[],3);'])
    eval(['delta_v_' dynamic_var '_c_' static_var '_array = var_array;'])
    clear var_array
end    

% make protein cohort figures
cm = jet(128);
color_array = cm(1:30:128,:);
constant_mf_fig = figure;
hold on
lgd_str = {};
e = [];
for i = 1:numel(prctile_vec)-1
    e = [e errorbar(tt_index/60, delta_v_tt_c_mf_mean(:,i), delta_v_tt_c_mf_ste(:,i),'Color',color_array(i,:),'CapSize',0)];
    lgd_str = [lgd_str{:} {['protein cohort ' num2str(i)]}];
end
legend(e,lgd_str{:})
xlabel('time (minutes)')
ylabel(['absolute ' protein_name 'enrichment (au)'])
saveas(constant_mf_fig,[figPath 'protein_cohort_plot.png'])

constant_tt_fig = figure;
hold on
lgd_str = {};
e = [];
for i = 1:numel(prctile_vec)-1
    e = [e errorbar(mf_index, delta_v_mf_c_tt_mean(:,i), delta_v_mf_c_tt_ste(:,i),'Color',color_array(i,:),'CapSize',0)];
    lgd_str = [lgd_str{:} {['time cohort ' num2str(i)]}];
end
legend(e,lgd_str{:})
xlabel(['average ' protein_name ' concentration (au)'])
ylabel(['absolute ' protein_name 'enrichment (au)'])
saveas(constant_tt_fig,[figPath 'time_cohort_plot.png'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test hunch that there's a simple linear relationship btw mf and
% enrichment
nan_filter = ~isnan(mf_protein_vec_dist) & ~isnan(delta_protein_vec_dist);
X = [ones(sum(nan_filter),1) mf_protein_vec_dist(nan_filter)'];
% linear model
beta = X \ delta_protein_vec_dist(nan_filter)';
% second order polynomial
p = polyfit(mf_protein_vec_dist(nan_filter),delta_protein_vec_dist(nan_filter),3);

enrichment_pd1 =  beta(1) + beta(2)*mf_index;
enrichment_pd2 = polyval(p,mf_index);

cm = jet(128);
lin_fit_fig = figure;
hold on
% scatter(mf_protein_vec_dist,delta_protein_vec_dist,'MarkerFaceColor',[.3 .3 .3],'MarkerFaceAlpha',.05)
e = errorbar(mf_index,delta_v_mf_c_tt_mean(:,end),delta_v_mf_c_tt_ste(:,end),'Color','black');
e.CapSize = 0;
plot(mf_index,enrichment_pd1,'Color',cm(10,:))
plot(mf_index,enrichment_pd2,'Color',cm(120,:))
legend('data','poly 1','poly 3','Location','northwest')
xlabel(['average ' protein_name ' concentration (au)'])
ylabel(['absolute ' protein_name ' enrichment (au)'])
grid on
saveas(lin_fit_fig,[figPath 'mf_enrichment_linear.png'])

% Not truly linear, but close enough for now. Let's see if we can re-create
% temporal profile using this model
% pd1_delta_protein_time_mat = NaN(numel(tt_index),NBoots);
pd2_delta_protein_time_mat = NaN(numel(tt_index),NBoots);
% pd1_delta_protein_vec = beta(1) + beta(2)*mf_protein_vec_dist;
pd2_delta_protein_vec = polyval(p,mf_protein_vec_dist);
for n = 1:NBoots
    s_ids = randsample(1:numel(delta_protein_vec_dist),numel(delta_protein_vec_dist),true);
%     dv_samp1 = pd1_delta_protein_vec(s_ids);     
    dv_samp2 = pd2_delta_protein_vec(s_ids);     
    time_samp = time_vec_dist(s_ids);
    for t = 1:numel(tt_index)
        if tt_index(t) < min(time_vec_dist) || tt_index(t) > max(time_vec_dist)
            continue
        end
        t_weights = exp(-.5*((time_samp-tt_index(t))/tt_sigma).^2);
%         pd1_delta_protein_time_mat(t,n) = nansum(dv_samp1.*t_weights)/nansum(t_weights);     
        pd2_delta_protein_time_mat(t,n) = nansum(dv_samp2.*t_weights)/nansum(t_weights); 
    end
end

% pd1_time_delta_mean = nanmean(pd1_delta_protein_time_mat,2);
% pd1_time_delta_ste = nanstd(pd1_delta_protein_time_mat,[],2);

pd2_time_delta_mean = nanmean(pd2_delta_protein_time_mat,2);
pd2_time_delta_ste = nanstd(pd2_delta_protein_time_mat,[],2);

% remake figure
delta_time_fig = figure;
yyaxis left
hold on
plot(tt_index/60,time_delta_mean,'LineWidth',1.5);
% plot(tt_index/60,pd1_time_delta_mean,'LineWidth',1.5);
plot(tt_index/60,pd2_time_delta_mean,'-','LineWidth',1.5,'Color',cm(120,:));

ylabel([protein_name ' enrichment at locus (au)'])
grid on
y_max = nanmax(time_delta_mean);
y_min = nanmin(time_delta_mean);
ylim([y_min-.1*y_min*sign(y_min) 1.1*y_max])
ax = gca;
ax.YColor = 'black';

yyaxis right
plot(tt_index/60, time_fluo_mean,'--','LineWidth',1,'Color',[.2 .2 .2 .4]);
ax = gca;
ax.YColor = [.2 .2 .2];
ylabel([gene_name ' activity (au)'])
ylim([100 240])
legend('enrichment trend (actual)','predicted enrichment (poly 2)','activity trend',...
    'Location','northeast')
xlabel('minutes')
title(['Enrichment of ' id_string])

saveas(delta_time_fig, [figPath write_string '_temporal_enrichment_pd.png'])


%%%%%%%%%%%%%%%%%%%%% Make Radial Profile Figure %%%%%%%%%%%%%%%%%%%%%%%%%%
% generate position reference matrices
r_sigma = .1;
% indexing vector for sampling
index_vec = 1:size(null_protein_snips_mixed,3);
% initialize profile arrays
r_spot_mat = NaN(NBoots,numel(dist_index));
r_null_mat = NaN(NBoots,numel(dist_index));

for n = 1:NBoots
    s_ids = randsample(index_vec,numel(index_vec),true);
    pt_spot = nanmean(spot_protein_snips_mixed(:,:,s_ids),3);
    pt_null = nanmean(null_protein_snips_mixed(:,:,s_ids),3);
    for r = 1:numel(dist_index)
        r_weights = exp(-.5*((r_ref-dist_index(r))/r_sigma).^2);
        r_spot_mat(n,r) = nansum(pt_spot(:).*r_weights(:)) ./ nansum(r_weights(:));
        r_null_mat(n,r) = nansum(pt_null(:).*r_weights(:)) ./ nansum(r_weights(:));
    end
end

r_spot_mean = nanmean(r_spot_mat);
r_spot_ste = nanstd(r_spot_mat);

r_null_mean = nanmean(r_null_mat);
r_null_ste = nanstd(r_null_mat);

% save table with radial profile results
radial_table = array2table([dist_index',r_spot_mean',r_spot_ste',r_null_mean',r_null_ste'],...
            'VariableNames',{'radial_distance','mean_locus_protein','ste_locus_protein','mean_control_protein','ste_control_protein'});
writetable(radial_table,[dataPath 'radial_profile_data.csv'])

% make figure
cm = Colormap_plot;
r_fig = figure;
r_ax = gca;
hold on
e = errorbar(dist_index,r_null_mean / r_null_mean(1),r_null_ste / r_null_mean(1),'Color','black','LineWidth',1.75);
e.CapSize = 0;
e = errorbar(dist_index,r_spot_mean / r_null_mean(1),r_spot_ste / r_null_mean(1),'Color',cm(35,:),'LineWidth',1.75);
e.CapSize = 0;
grid on
r_ax.XLabel.String = 'radius (\mu m)';
r_ax.YLabel.String ='relative enrichment';
legend('control','active locus')
title(['Radial Concentration Profile (' id_string ')'])
r_ax.XLim = [0 1.2];
r_ax.YLim = [relEnrich_lb relEnrich_ub];
saveas(r_fig, [figPath write_string '_radial_enrichment.png'])

% make paper fig in PBoC style
r_ax = gca;
StandardFigurePBoC(e,r_ax);
saveas(r_fig, [paperFigPath write_string '_radial_enrichment.pdf'])



% % take bootstrap samples
% ap_index = unique(round(snip_ap_vec_dist));
% delta_protein_ap_mat = NaN(numel(ap_index),NBoots);
% spot_protein_ap_mat = NaN(numel(ap_index),NBoots);
% null_protein_ap_mat = NaN(numel(ap_index),NBoots);
% ap_sigma= 2; % seconds
% for n = 1:NBoots
%     s_ids = randsample(1:numel(delta_protein_vec_dist),numel(delta_protein_vec_dist),true);
%     dv_samp = delta_protein_vec_dist(s_ids);    
%     nn_samp = null_protein_vec_dist(s_ids);
%     sp_samp = spot_protein_vec_dist(s_ids);
%     ap_samp = snip_ap_vec_dist(s_ids);
%     for a = 1:numel(ap_index)        
%         ap_weights = exp(-.5*((ap_samp-ap_index(a))/ap_sigma).^2);
%         delta_protein_ap_mat(a,n) = nansum(dv_samp.*ap_weights)/nansum(ap_weights);        
%         null_protein_ap_mat(a,n) =  nansum(nn_samp.*ap_weights)/nansum(ap_weights);
%         spot_protein_ap_mat(a,n) = nansum(sp_samp.*ap_weights)/nansum(ap_weights);
%     end
% end
% % calculate average and standard error
% ap_delta_mean = nanmean(delta_protein_ap_mat./null_protein_ap_mat,2);
% ap_delta_ste = nanstd(delta_protein_ap_mat./null_protein_ap_mat,[],2);
% ap_null_mean = nanmean(null_protein_ap_mat,2);
% ap_null_ste = nanstd(null_protein_ap_mat,[],2);
% ap_spot_mean = nanmean(spot_protein_ap_mat,2);
% ap_spot_ste = nanstd(spot_protein_ap_mat,[],2);
% 
% cm = Colormap_plot;
% delta_ap_fig = figure;
% 
% yyaxis left
% e = errorbar(ap_index,100*ap_delta_mean,100 * ap_delta_ste,'-','Color',cm(60,:)/1.2,'LineWidth',1.5);
% e.CapSize = 0;
% ylabel('% enrichment at locus')
% y_max = nanmax(100 * ap_delta_mean);
% y_min = nanmin(100 * ap_delta_mean);
% ylim([y_min-.1*y_min*sign(y_min) 1.1*y_max])
% 
% yyaxis right
% plot(ap_index, ap_spot_mean,'-','Color',cm(30,:));
% hold on
% plot(ap_index, ap_null_mean,'-','Color','black');
% ylabel('concentration (au)')
% 
% legend('percent enrichment','locus concentration', 'control concentration','Location','northwest')
% 
% title(['Percent Enrichment of ' id_string])
% xlabel('% AP')
% grid on
% saveas(delta_ap_fig, [FigPath write_string '_ap_percent_enrichment.png'])
end