clear
close all

projectName = 'Bcd-GFP_hbMS2-mCh_Airy_fast';

liveProject = LiveEnrichmentProject(projectName);
resultsRoot = [liveProject.dataPath filesep];
FigurePath = liveProject.figurePath;

projectNameCell = {'Bcd-GFP_hbMS2-mCh_Airy_fast','Bcd-GFP_hbMS2-mCh_Airy_fast_int','Bcd-GFP_hbMS2-mCh_NoAiry_02','Bcd-GFP_hbP2P-mCh'};
master_struct = struct;

for p = 1:length(projectNameCell)
    projectName = projectNameCell{p};
    liveProject = LiveEnrichmentProject(projectName);
    resultsRoot = [liveProject.dataPath filesep];
    
    % load data
    load([resultsRoot 'spot_struct.mat'])
    load([resultsRoot 'spot_struct_protein.mat'])
    load([resultsRoot 'proteinSamplingInfo.mat'])
    load([resultsRoot 'snip_data.mat'])
    master_struct(p).spot_struct = spot_struct;
    master_struct(p).spot_struct_protein = spot_struct_protein;
    master_struct(p).pixelSize = liveProject.includedExperiments{1}.pixelSize_um;
    master_struct(p).snip_data = snip_data;
    master_struct(p).projectName = projectName;
    clear spot_struct
end
%% Extract the snip date
DistLim = 0.8;
for p = 1:length(projectNameCell)
    spot_struct_protein = master_struct(p).spot_struct_protein;
    snip_data = master_struct(p).snip_data;

    null_protein_vec = [spot_struct_protein.edge_null_protein_vec];
    dist_vec = [spot_struct_protein.spot_edge_dist_vec]*master_struct(p).pixelSize;
    
    % extract snips
    [master_struct(p).spot_protein_snips, master_struct(p).spot_mcp_snips, ...
     master_struct(p).edge_control_protein_snips, master_struct(p).edge_control_mcp_snips] = ...
                                          generateProteinSnips(snip_data, ~isnan(null_protein_vec)&dist_vec>=DistLim);

    % average
    master_struct(p).spot_protein_snip_mean = nanmean(master_struct(p).spot_protein_snips,3);
    master_struct(p).edge_control_protein_snip_mean = nanmean(master_struct(p).edge_control_protein_snips,3);
    master_struct(p).edge_control_mcp_snip_mean = nanmean(master_struct(p).edge_control_mcp_snips,3);
    master_struct(p).spot_mcp_snip_mean = nanmean(master_struct(p).spot_mcp_snips,3);
end
%% Calculate radial concentration profiles
NBoots = 10;

for p = 1:length(projectNameCell)

    snip_size = size(master_struct(p).spot_protein_snips,1);
    [y_ref, x_ref] = meshgrid(1:snip_size,1:snip_size);
    r_ref = sqrt((x_ref - ceil(snip_size/2)).^2 + (y_ref - ceil(snip_size/2)).^2)*master_struct(p).pixelSize;
    dist_index = unique(r_ref);
    % generate smoothly interpolated axis for plot
    dist_plot_axis = linspace(dist_index(1),dist_index(end),50);

    % indexing vector for sampling
    index_vec = 1:size(master_struct(p).edge_control_protein_snips,3);
    % define scale for moving average
    r_sigma = master_struct(p).pixelSize;
    maxBootSize = 5e3;
    % initialize profile arrays
    r_spot_mat = NaN(NBoots,numel(dist_plot_axis));
    r_control_mat = NaN(NBoots,numel(dist_plot_axis));

    for n = 1:NBoots
        s_ids = randsample(index_vec,min([maxBootSize length(index_vec)]),true);
        pt_spot = nanmean(master_struct(p).spot_protein_snips(:,:,s_ids),3);
        pt_null = nanmean(master_struct(p).edge_control_protein_snips(:,:,s_ids),3);
        for r = 1:length(dist_plot_axis)
            r_weights = exp(-.5*((r_ref-dist_plot_axis(r))/r_sigma).^2);
            r_spot_mat(n,r) = nansum(pt_spot(:).*r_weights(:)) ./ nansum(r_weights(:));
            r_control_mat(n,r) = nansum(pt_null(:).*r_weights(:)) ./ nansum(r_weights(:));
        end
    end

    master_struct(p).r_spot_mean = nanmean(r_spot_mat);
    master_struct(p).r_spot_ste = nanstd(r_spot_mat);

    master_struct(p).r_control_mean = nanmean(r_control_mat);
    master_struct(p).r_control_ste = nanstd(r_control_mat);
    master_struct(p).dist_plot_axis = dist_plot_axis;
end

%% make plots

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First, interpolated vs raw airyscan

dist_plot_axis = master_struct(1).dist_plot_axis;
r_control_mean = master_struct(1).r_control_mean;
r_control_ste = master_struct(1).r_control_ste;

r_spot_mean = master_struct(1).r_spot_mean;
r_spot_ste = master_struct(1).r_spot_ste;

r_spot_mean_int = master_struct(2).r_spot_mean;
r_spot_ste_int = master_struct(2).r_spot_ste;

[~,norm_ind] = min(abs(dist_plot_axis-1));
norm_shift = r_control_mean(1) / r_control_mean(norm_ind) - 1;

r_fig = figure;
cm2 = brewermap([],'Set2');
r_ax = gca;
hold on
e = errorbar(dist_plot_axis,r_control_mean / r_control_mean(norm_ind)-norm_shift,r_control_ste / r_control_mean(norm_ind),'Color','black','LineWidth',1.75);
e.CapSize = 0;
e = errorbar(dist_plot_axis,r_spot_mean / r_spot_mean(norm_ind)-norm_shift,r_spot_ste / r_spot_mean(norm_ind),'Color',cm2(3,:),'LineWidth',1.75);
e.CapSize = 0;
e = errorbar(dist_plot_axis,r_spot_mean_int / r_spot_mean_int(norm_ind)-norm_shift,r_spot_ste_int / r_spot_mean_int(norm_ind),'Color',cm2(2,:),'LineWidth',1.75);
e.CapSize = 0;
grid off
r_ax.XLabel.String = 'radius (\mu m)';
r_ax.YLabel.String ='relative enrichment';
legend('control','active locus','active locus (interpolated)')
% title(['Radial Concentration Profile (' inputString ')'])
r_ax.XLim = [0 1];
% r_ax.YLim = [.98 1.08];
StandardFigure([],r_ax);
% set(gca,'ytick',0.95:0.05:1.25)
% r_ax.YLim = [relEnrich_lb relEnrich_ub];
saveas(r_fig, [FigurePath '_radial_enrichment_int_comparison.png'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now compare old and new confocal sets


dist_plot_axis = master_struct(3).dist_plot_axis;
r_control_mean = master_struct(3).r_control_mean;
r_control_ste = master_struct(3).r_control_ste;

r_spot_mean = master_struct(3).r_spot_mean;
r_spot_ste = master_struct(3).r_spot_ste;

r_spot_mean_old = master_struct(4).r_spot_mean;
r_spot_ste_old = master_struct(4).r_spot_ste;

[~,norm_ind] = min(abs(dist_plot_axis-1));
norm_shift = r_control_mean(1) / r_control_mean(norm_ind) - 1;

r_fig = figure;
cm2 = brewermap([],'Set2');
r_ax = gca;
hold on
e = errorbar(dist_plot_axis,r_control_mean / r_control_mean(norm_ind)-norm_shift,r_control_ste / r_control_mean(norm_ind),'Color','black','LineWidth',1.75);
e.CapSize = 0;
e = errorbar(dist_plot_axis,r_spot_mean / r_spot_mean(norm_ind)-norm_shift,r_spot_ste / r_spot_mean(norm_ind),'Color',cm2(3,:),'LineWidth',1.75);
e.CapSize = 0;
e = errorbar(dist_plot_axis,r_spot_mean_old / r_spot_mean_old(norm_ind)-norm_shift,r_spot_ste_old / r_spot_mean_int(norm_ind),'Color',cm2(4,:),'LineWidth',1.75);
e.CapSize = 0;
grid off
r_ax.XLabel.String = 'radius (\mu m)';
r_ax.YLabel.String ='relative enrichment';
legend('control','active locus (Zeiss 980)','active locus (Leica SP8)')
% title(['Radial Concentration Profile (' inputString ')'])
r_ax.XLim = [0 1];
% r_ax.YLim = [.98 1.08];
StandardFigure([],r_ax);
% set(gca,'ytick',0.95:0.05:1.25)
% r_ax.YLim = [relEnrich_lb relEnrich_ub];
saveas(r_fig, [FigurePath '_radial_enrichment_old_new_comparison.png'])