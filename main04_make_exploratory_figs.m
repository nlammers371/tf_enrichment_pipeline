% Script conduct locus enrichment analyses
function main04_make_exploratory_figs(project, varargin)

close all

DataPath = ['../../dat/' project '/'];
FigPath = ['../../fig/' project '/'];
ControlType = 'edge'; % specify type of control to use
ROIRadius = .3; % radus (um) of region used to query and compare TF concentrations
DistLim = .6; % min distance from edge permitted (um)
NBoots = 100; % number of bootstrap samples to use for estimating SE
ManualDistThreshold = 0;
Colormap_plot = jet(128); %specifies the colomap used to make plots/graphs
Colormap_heat = viridis(128); %specifies the colormap used to make heatmaps

for i = 1:numel(varargin)
    if strcmpi(varargin{i}, 'dropboxFolder')        
        DataPath = [varargin{i+1} '/ProcessedEnrichmentData/' project '/'];
        FigPath = [varargin{i+1} '/LocalEnrichmentFigures/' project '/'];        
    elseif ischar(varargin{i})
        if ismember(varargin{i}, {'ControlType','ROIRadius','DistLim','NBoots','ManualDistThreshold'}) 
            eval([varargin{i} '=varargin{i+1};']);
        end
    end
end
FigPath = [FigPath ControlType '/'];
mkdir(FigPath)

% Load analysis data
load([DataPath 'nucleus_struct_protein.mat']);

% extract protein, gene, fluorophore info
underscores = strfind(project,'_');
protein_name = project(1:underscores(1)-1);
protein_fluor = project(underscores(1)+1:underscores(2)-1);
gene_name = project(underscores(2)+1:underscores(3)-1);
if numel(underscores) == 3
    ind = numel(project);
else
    ind = underscores(4)-1;
end
gene_fluor = project(underscores(3)+1:end);

id_string = [protein_name '-' protein_fluor ' : ' gene_name '-' gene_fluor]; 
write_string = [protein_name '-' protein_fluor '__' gene_name '-' gene_fluor]; 
% Generate distance vector for filtering snip stacks
PixelSize = nucleus_struct_protein(1).PixelSize;
snip_dist_vec = [];
snip_time_vec = [];
snip_ap_vec = [];
for i = 1:numel(nucleus_struct_protein)
    snip_frame_vec = nucleus_struct_protein(i).snip_frame_vec;
    nc_dist_vec = nucleus_struct_protein(i).(['spot_' ControlType '_dist_vec']);
    time_vec = nucleus_struct_protein(i).time;
    ap_vec = nucleus_struct_protein(i).ap_vector;
    nc_frames = nucleus_struct_protein(i).frames;
    snip_dist_vec = [snip_dist_vec PixelSize*nc_dist_vec(ismember(nc_frames,snip_frame_vec))];
    snip_time_vec = [snip_time_vec time_vec(ismember(nc_frames,snip_frame_vec))];
    snip_ap_vec = [snip_ap_vec ap_vec(ismember(nc_frames,snip_frame_vec))];
end    
% fixing error in a data set
if strcmpi(project,'Hb_NbGFP_hbBAC_mCherry')
    ap_ft = snip_ap_vec > .5;
    snip_ap_vec(ap_ft) = 1 - snip_ap_vec(ap_ft);
end
% invert distance vector if we're using centroid metric
if strcmpi(ControlType, 'centroid')
    snip_dist_vec = max(snip_dist_vec) - snip_dist_vec;
end
% Snip stacks
spot_protein_snips = cat(3,nucleus_struct_protein.spot_protein_snips);
null_protein_snips = cat(3,nucleus_struct_protein.([ControlType '_null_protein_snips']));
spot_mcp_snips = cat(3,nucleus_struct_protein.spot_mcp_snips);
null_mcp_snips = cat(3,nucleus_struct_protein.([ControlType '_null_mcp_snips']));

% Make r reference array
snip_size = size(spot_protein_snips,1);
[y_ref, x_ref] = meshgrid(1:snip_size,1:snip_size);
r_ref = sqrt((x_ref-ceil(snip_size/2)).^2 + (y_ref-ceil(snip_size/2)).^2)*PixelSize;

r_ft = 1*repmat(r_ref <= ROIRadius,1,1,size(spot_protein_snips,3));
spot_protein_vec = reshape(nanmean(nanmean(r_ft.*spot_protein_snips,1),2),1,[]);
null_protein_vec = reshape(nanmean(nanmean(r_ft.*null_protein_snips,1),2),1,[]);

% First look for presence of edge artifact
dist_sigma = .1; %(um)
protein_delta_vec =  spot_protein_vec - null_protein_vec;
dist_index = 0:.1:floor(prctile(snip_dist_vec,99)*10)/10;

delta_dist_mat = NaN(numel(dist_index),NBoots);
null_dist_mat = NaN(numel(dist_index),NBoots);
spot_dist_mat = NaN(numel(dist_index),NBoots);

for n = 1:NBoots
    s_ids = randsample(1:numel(protein_delta_vec),numel(protein_delta_vec),true);
    dv_samp = protein_delta_vec(s_ids);    
    nn_samp = spot_protein_vec(s_ids);
    sp_samp = null_protein_vec(s_ids);
    dist_samp = snip_dist_vec(s_ids);        
    for t = 1:numel(dist_index)
        d_weights = exp(-.5*((dist_samp-dist_index(t))/dist_sigma).^2);
        delta_dist_mat(t,n) = nansum(dv_samp.*d_weights) / nansum(d_weights);
        null_dist_mat(t,n) = (nansum(nn_samp.*d_weights) / nansum(d_weights));
        spot_dist_mat(t,n) = (nansum(sp_samp.*d_weights) / nansum(d_weights));
    end
end
% calcualte mean and standard error
delta_dist_mean = nanmean(delta_dist_mat,2);
delta_dist_ste = nanstd(delta_dist_mat,[],2);

null_dist_mean = nanmean(null_dist_mat,2);
null_dist_ste = nanstd(null_dist_mat,[],2);

spot_dist_mean = nanmean(spot_dist_mat,2);
spot_dist_ste = nanstd(spot_dist_mat,[],2);


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
saveas(delta_dist_fig, [FigPath write_string '_edge_artifact_plot.png'])

% Apply distance filter and make histogram figures
null_protein_vec_dist = null_protein_vec(snip_dist_vec>=DistLim);
spot_protein_vec_dist = spot_protein_vec(snip_dist_vec>=DistLim);
delta_protein_vec_dist = spot_protein_vec_dist-null_protein_vec_dist;
delta_norm_protein_vec_dist = (spot_protein_vec_dist-null_protein_vec_dist) / nanmean(null_protein_vec_dist);
snip_time_vec_dist = snip_time_vec(snip_dist_vec>=DistLim);
snip_ap_vec_dist = 100*snip_ap_vec(snip_dist_vec>=DistLim);
% plot raw distributions
lb = .9*prctile([null_protein_vec_dist spot_protein_vec_dist],1);
ub = 1.1*prctile([null_protein_vec_dist spot_protein_vec_dist],99);
bins = linspace(lb,ub,100);

spot_ct = histc(spot_protein_vec_dist,bins);
null_ct = histc(null_protein_vec_dist,bins);
spot_ct = spot_ct / sum(spot_ct);
null_ct = null_ct / sum(null_ct);

raw_hist_fig = figure;
cm = Colormap_plot;
hold on
% hist plots
yyaxis left
b1 = bar(bins,null_ct,1,'FaceAlpha',.5,'FaceColor',cm(30,:));
b2 = bar(bins,spot_ct,1,'FaceAlpha',.5,'FaceColor',cm(100,:));
ylabel('share')
ax = gca;
ax.YColor = 'black';
% cumulative plots
yyaxis right
plot(bins,cumsum(null_ct),'-','LineWidth',1.5,'Color',cm(30,:));
plot(bins,cumsum(spot_ct),'-','LineWidth',1.5,'Color',cm(100,:));
ylabel('cumulative share')
ax = gca;
ax.YColor = 'black';

legend('control',['active locus (' gene_name ')'])
xlabel('concentration (au)')
title([protein_name '-' protein_fluor ' Concentrations'])
grid on
saveas(raw_hist_fig,[FigPath 'hist_plots_' write_string '_dist_' num2str(DistLim) '.png']); 

% plot delta distribution
pd = fitdist(delta_norm_protein_vec_dist','Normal');
lb = .9*prctile(delta_norm_protein_vec_dist,1);
ub = 1.1*prctile(delta_norm_protein_vec_dist,99);
bins = linspace(lb,ub,100);

delta_ct = histc(delta_norm_protein_vec_dist,bins);
delta_ct = delta_ct / sum(delta_ct);
   

delta_hist_fig = figure;
cm = Colormap_plot;
hold on
% hist plots
b = bar(bins,delta_ct,1,'FaceAlpha',.5,'FaceColor',cm(60,:)/1.2);
norm_vec = exp(-.5*((bins-pd.mu)/pd.sigma).^2);
p = plot(bins,norm_vec/sum(norm_vec),'LineWidth',2);
ylabel('share')
legend([b,p],['\Delta F ('  gene_name ')'], ...
    ['Gaussian Fit (\mu=' num2str(pd.mu) ' \sigma=' num2str(pd.sigma)],'Location','best')
xlabel('concentration (au)')
title(['Enrichment at Active Locus Relative to Control (' id_string ')'])
grid on
saveas(delta_hist_fig,[FigPath 'delta_hist_' write_string '_dist_' num2str(DistLim) '.png']); 


%%%%%%%%%%%%%%%%%%%%%%%% Compare protein snippets %%%%%%%%%%%%%%%%%%%%%%%%%
% apply distance filter
spot_protein_snips_dist = spot_protein_snips(:,:,snip_dist_vec>=DistLim);
null_protein_snips_dist = null_protein_snips(:,:,snip_dist_vec>=DistLim);
spot_mcp_snips_dist = spot_mcp_snips(:,:,snip_dist_vec>=DistLim);
null_mcp_snips_dist = null_mcp_snips(:,:,snip_dist_vec>=DistLim);

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

spot_pt_snip_fig = figure;
colormap(Colormap_heat)
imagesc(spot_protein_snip_mean)
title([protein_name '-' protein_fluor ' at ' 'Active ' gene_name ' Locus'])
caxis([lb ub])
h = colorbar;
ylabel('\mum','FontSize',12)
xlabel('\mum','FontSize',12)
ylabel(h,[protein_name '-' protein_fluor ' concentration (au)'],'FontSize',12)
eval(xtick_string)
eval(ytick_string)
saveas(spot_pt_snip_fig,[FigPath write_string '_mean_pt_snippet_spot.png']);    

null_pt_snip_fig = figure;
colormap(Colormap_heat)
imagesc(null_protein_snip_mean)
title([protein_name '-' protein_fluor ' at Control Locus'])
caxis([lb ub])
h = colorbar;
ylabel('\mum','FontSize',12)
xlabel('\mum','FontSize',12)
ylabel(h,[protein_name '-' protein_fluor ' concentration (au)'],'FontSize',12)
eval(xtick_string)
eval(ytick_string)
saveas(null_pt_snip_fig,[FigPath write_string '_mean_pt_snippet_null.png']);  

% Make diff image
rel_protein_snip_mean = (spot_protein_snip_mean) ./ null_protein_snip_mean;
rel_pt_snip_fig = figure;
colormap(Colormap_heat)
imagesc(rel_protein_snip_mean)
title(['Relative ' protein_name '-' protein_fluor ' Enrichment at Active ' gene_name ' Locus'])
% caxis([lb ub])
h = colorbar;
ylabel('\mum','FontSize',12)
xlabel('\mum','FontSize',12)
ylabel(h,[protein_name '-' protein_fluor ' fold enrichment'],'FontSize',12)
eval(xtick_string)
eval(ytick_string)
saveas(rel_pt_snip_fig,[FigPath write_string '_mean_pt_snippet_rel.png']);    

% compare mcp snippets

null_mean_fluo = nanmean(null_mcp_snips_mixed,3);
spot_mean_fluo = nanmean(spot_mcp_snips_mixed,3);
ub = prctile([reshape(null_mean_fluo,1,[]) reshape(spot_mean_fluo,1,[])],99);
lb = prctile([reshape(null_mean_fluo,1,[]) reshape(spot_mean_fluo,1,[])],1);

fluo_snippet_spot_fig = figure;
colormap(Colormap_heat)
imagesc(spot_mean_fluo)
title([gene_fluor ' Intensity at Active Locus (' gene_name ')'])
caxis([lb ub])
h = colorbar;
ylabel('\mum','FontSize',12)
xlabel('\mum','FontSize',12)
ylabel(h,[gene_name ' expression (au)'],'FontSize',12)
eval(xtick_string)
eval(ytick_string)
saveas(fluo_snippet_spot_fig,[FigPath write_string '_mean_fluo_snippet_spot.png']);    

fluo_snippet_null_fig = figure;
colormap(Colormap_heat)
imagesc(null_mean_fluo)
title([gene_fluor ' Intensity at Control Locus'])
caxis([lb ub])
h = colorbar;
ylabel('\mum','FontSize',12)
xlabel('\mum','FontSize',12)
ylabel(h,[gene_name ' expression (au)'],'FontSize',12)
eval(xtick_string)
eval(ytick_string)
saveas(fluo_snippet_null_fig,[FigPath write_string '_mean_fluo_snippet_null.png']);    


%%%%%%% plot relative enrichment at locus as a function of time %%%%%%%%%%%

% take bootstrap samples
time_index = 0:20:round(max(snip_time_vec_dist));
delta_protein_time_mat = NaN(numel(time_index),NBoots);
spot_protein_time_mat = NaN(numel(time_index),NBoots);
null_protein_time_mat = NaN(numel(time_index),NBoots);
t_window = 120; % seconds
for n = 1:NBoots
    s_ids = randsample(1:numel(delta_protein_vec_dist),numel(delta_protein_vec_dist),true);
    dv_samp = delta_protein_vec_dist(s_ids);    
    nn_samp = null_protein_vec_dist(s_ids);
    sp_samp = spot_protein_vec_dist(s_ids);
    time_samp = snip_time_vec_dist(s_ids);
    for t = 1:numel(time_index)
        if time_index(t) < min(snip_time_vec_dist) || time_index(t) > max(snip_time_vec_dist)
            continue
        end
        t_weights = exp(-.5*((time_samp-time_index(t))/t_window).^2);
        delta_protein_time_mat(t,n) = nansum(dv_samp.*t_weights)/nansum(t_weights);        
        null_protein_time_mat(t,n) =  nansum(nn_samp.*t_weights)/nansum(t_weights);
        spot_protein_time_mat(t,n) = nansum(sp_samp.*t_weights)/nansum(t_weights);
    end
end
% calculate average and standard error
time_delta_mean = nanmean(delta_protein_time_mat./null_protein_time_mat,2);
time_delta_ste = nanstd(delta_protein_time_mat./null_protein_time_mat,[],2);
time_null_mean = nanmean(null_protein_time_mat,2);
time_null_ste = nanstd(null_protein_time_mat,[],2);
time_spot_mean = nanmean(spot_protein_time_mat,2);
time_spot_ste = nanstd(spot_protein_time_mat,[],2);


delta_time_fig = figure;
yyaxis left
e = errorbar(time_index/60,100 * time_delta_mean,100 * time_delta_ste,'-','Color',cm(60,:)/1.2,'LIneWidth',1.5);
e.CapSize = 0;

ylabel('% enrichment at locus')
grid on
y_max = nanmax(100 * time_delta_mean);
y_min = nanmin(100 * time_delta_mean);
ylim([y_min-.1*y_min*sign(y_min) 1.1*y_max])

yyaxis right
plot(time_index/60, time_spot_mean,'-','Color',cm(30,:));
hold on
plot(time_index/60, time_null_mean,'-','Color','black');
ylabel('concentration (au)')

legend('percent enrichment','locus concentration', 'control concentration','Location','northwest')
xlabel('minutes')
title(['Percent Enrichment of ' id_string])

saveas(delta_time_fig, [FigPath write_string '_time_percent_enrichment.png'])


%%%%%%% plot relative enrichment at locus as a function of AP %%%%%%%%%%%

% take bootstrap samples
ap_index = unique(round(snip_ap_vec_dist));
delta_protein_ap_mat = NaN(numel(ap_index),NBoots);
spot_protein_ap_mat = NaN(numel(ap_index),NBoots);
null_protein_ap_mat = NaN(numel(ap_index),NBoots);
ap_sigma= 2; % seconds
for n = 1:NBoots
    s_ids = randsample(1:numel(delta_protein_vec_dist),numel(delta_protein_vec_dist),true);
    dv_samp = delta_protein_vec_dist(s_ids);    
    nn_samp = null_protein_vec_dist(s_ids);
    sp_samp = spot_protein_vec_dist(s_ids);
    ap_samp = snip_ap_vec_dist(s_ids);
    for a = 1:numel(ap_index)        
        ap_weights = exp(-.5*((ap_samp-ap_index(a))/ap_sigma).^2);
        delta_protein_ap_mat(a,n) = nansum(dv_samp.*ap_weights)/nansum(ap_weights);        
        null_protein_ap_mat(a,n) =  nansum(nn_samp.*ap_weights)/nansum(ap_weights);
        spot_protein_ap_mat(a,n) = nansum(sp_samp.*ap_weights)/nansum(ap_weights);
    end
end
% calculate average and standard error
ap_delta_mean = nanmean(delta_protein_ap_mat./null_protein_ap_mat,2);
ap_delta_ste = nanstd(delta_protein_ap_mat./null_protein_ap_mat,[],2);
ap_null_mean = nanmean(null_protein_ap_mat,2);
ap_null_ste = nanstd(null_protein_ap_mat,[],2);
ap_spot_mean = nanmean(spot_protein_ap_mat,2);
ap_spot_ste = nanstd(spot_protein_ap_mat,[],2);

cm = Colormap_plot;
delta_ap_fig = figure;

yyaxis left
e = errorbar(ap_index,100*ap_delta_mean,100 * ap_delta_ste,'-','Color',cm(60,:)/1.2,'LineWidth',1.5);
e.CapSize = 0;
ylabel('% enrichment at locus')
y_max = nanmax(100 * ap_delta_mean);
y_min = nanmin(100 * ap_delta_mean);
ylim([y_min-.1*y_min*sign(y_min) 1.1*y_max])

yyaxis right
plot(ap_index, ap_spot_mean,'-','Color',cm(30,:));
hold on
plot(ap_index, ap_null_mean,'-','Color','black');
ylabel('concentration (au)')

legend('percent enrichment','locus concentration', 'control concentration','Location','northwest')

title(['Percent Enrichment of ' id_string])
xlabel('% AP')
grid on
saveas(delta_ap_fig, [FigPath write_string '_ap_percent_enrichment.png'])

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

% make figure
cm = Colormap_plot;
r_fig = figure;
hold on
e = errorbar(dist_index,r_null_mean / r_null_mean(1),r_null_ste / r_null_mean(1),'Color','black','LineWidth',1.75);
e.CapSize = 0;
e = errorbar(dist_index,r_spot_mean / r_null_mean(1),r_spot_ste / r_null_mean(1),'Color',cm(35,:),'LineWidth',1.75);
e.CapSize = 0;
grid on
xlabel('radius (\mu m)')
ylabel('relative enrichment')
legend('control','active locus')
title(['Radial Concentration Profile (' id_string ')'])

saveas(r_fig, [FigPath write_string '_radial_enrichment.png'])
