% Script to conduct exploratory analyses relating to spatiotemporal
% dynamics of enrichment

function make_comparison_plots(ProjectCell, varargin)
close all

DataRoot = '../dat/';
ControlType = 'edge'; % specify type of control to use
ROIRadius = .3; % radus (um) of region used to query and compare TF concentrations
DistLim = .6; % min distance from edge permitted (um)
NBoots = 100; % number of bootstrap samples to use for estimating SE
for i = 1:numel(varargin)
    if strcmpi(varargin{i}, 'DropboxFolder')        
        DataRoot = [varargin{i+1} '/ProcessedEnrichmentData/'];
        FigPath = [varargin{i+1} '/LocalEnrichmentFigures/' ProjectCell{1} '/'];
    elseif strcmpi(varargin{i}, 'ControlType')  
        ControlType = varargin{i+1};
    elseif strcmpi(varargin{i}, 'ROIRadius')  
        ROIRadius = varargin{i+1};
    elseif strcmpi(varargin{i}, 'DistLim')  
        DistLim = varargin{i+1};
    elseif strcmpi(varargin{i}, 'NBoots')  
        NBoots = varargin{i+1};        
    elseif strcmpi(varargin{i}, 'ManualDistThreshold')  
        ManualDistThreshold = 1;    
    end    
end
% core variables
% primary_project = 'Hb_NbGFP_hbBAC_mCherry';
% control_project = 'Hb_NbGFP_snaBAC_mCherry';

% ProjectCell = [{'Hb_NbGFP_hbBAC_mCherry'} {'Hb_NbGFP_snaBAC_mCherry'}];
% DropboxFolder = 'E:\Nick\Dropbox (Garcia Lab)\';

% DataRoot = [DropboxFolder '/ProcessedEnrichmentData/'];
if numel(ProjectCell) == 2
    FigPath = [FigPath '/comparison_' ProjectCell{2} '/'];
end

FigPath = [FigPath ControlType '/'];
mkdir(FigPath)

% DataPathVec = [{[DataRoot primary_project '/']} {[DataRoot control_project '/']}];

master_struct = struct;

% iterate through projects 
for i = 1:numel(ProjectCell)
    % Load analysis data
    load([DataRoot '/' ProjectCell{i} '/nucleus_struct_protein.mat']);

    % Extract protein, gene, fluorophore info
    project = ProjectCell{i};
    underscores = strfind(project,'_');
    protein_name =project(1:underscores(1)-1);
    protein_fluor = project(underscores(1)+1:underscores(2)-1);
    gene_name = project(underscores(2)+1:underscores(3)-1);
    if numel(underscores) == 3
        ind = numel(project);
    else
        ind = underscores(4)-1;
    end
    gene_fluor = project(underscores(3)+1:end);
    
    master_struct(i).gene_name = gene_name;
    master_struct(i).protein_name = protein_name;
    master_struct(i).id_string = [protein_name '-' protein_fluor ' : ' gene_name '-' gene_fluor]; 
    master_struct(i).write_string = [protein_name '-' protein_fluor '__' gene_name '-' gene_fluor]; 
    % Generate distance vector for filtering snip stacks
    PixelSize = nucleus_struct_protein(1).PixelSize;
    snip_dist_vec = [];
    snip_time_vec = [];
    snip_ap_vec = [];
    snip_mf_protein_vec = [];
    for j = 1:numel(nucleus_struct_protein)
        snip_frame_vec = nucleus_struct_protein(j).snip_frame_vec;
        nc_dist_vec = nucleus_struct_protein(j).(['spot_' ControlType '_dist_vec']);
        time_vec = nucleus_struct_protein(j).time;
        ap_vec = nucleus_struct_protein(j).ap_vector;
        nc_protein = nucleus_struct_protein(j).protein;

        nc_frames = nucleus_struct_protein(j).frames;
        snip_dist_vec = [snip_dist_vec PixelSize*nc_dist_vec(ismember(nc_frames,snip_frame_vec))];
        snip_time_vec = [snip_time_vec time_vec(ismember(nc_frames,snip_frame_vec))];
        snip_ap_vec = [snip_ap_vec ap_vec(ismember(nc_frames,snip_frame_vec))];
        snip_mf_protein_vec = [snip_mf_protein_vec nc_protein(ismember(nc_frames,snip_frame_vec))];
    end  
    % invert distance vector if we're using centroid metric
    if strcmpi(ControlType, 'centroid')
        snip_dist_vec = max(snip_dist_vec) - snip_dist_vec;
    end
    % fixing error in a data set
    if strcmpi(ProjectCell{i},'Hb_NbGFP_hbBAC_mCherry')
        ap_ft = snip_ap_vec > .5;
        snip_ap_vec(ap_ft) = 1 - snip_ap_vec(ap_ft);
    end
    % Snip stacks
    spot_protein_snips = cat(3,nucleus_struct_protein.spot_protein_snips);
    spot_protein_snips = spot_protein_snips(:,:,snip_dist_vec>=DistLim);
    null_protein_snips = cat(3,nucleus_struct_protein.([ControlType '_null_protein_snips']));
    null_protein_snips = null_protein_snips(:,:,snip_dist_vec>=DistLim);
    
    % Make r reference array
    snip_size = size(spot_protein_snips,1);
    [y_ref, x_ref] = meshgrid(1:snip_size,1:snip_size);
    r_ref = sqrt((x_ref-ceil(snip_size/2)).^2 + (y_ref-ceil(snip_size/2)).^2)*PixelSize;

    % Make protein cencentration vectors
    r_ft = 1*repmat(r_ref <= ROIRadius,1,1,size(spot_protein_snips,3));
    master_struct(i).spot_protein_vec = reshape(nanmean(nanmean(r_ft.*spot_protein_snips,1),2),1,[]);
    master_struct(i).null_protein_vec = reshape(nanmean(nanmean(r_ft.*null_protein_snips,1),2),1,[]);
    master_struct(i).delta_protein_vec = master_struct(i).spot_protein_vec - master_struct(i).null_protein_vec;
    master_struct(i).project = ProjectCell{i};
    master_struct(i).time_vec = snip_time_vec(snip_dist_vec>=DistLim)/60;
    master_struct(i).ap_vec = snip_ap_vec(snip_dist_vec>=DistLim)*100;
    master_struct(i).mf_protein_vec = snip_mf_protein_vec(snip_dist_vec>=DistLim);
    master_struct(i).max_pt = nanmax(master_struct(i).mf_protein_vec);
    master_struct(i).mf_protein_vec = master_struct(i).mf_protein_vec / master_struct(1).max_pt;
    master_struct(i).spot_protein_snips = spot_protein_snips;
    master_struct(i).null_protein_snips = null_protein_snips;
end



%% Examine how local TF levels scale with nucleus average
PlotVariableCell = [{'mf_protein_vec'}, {'ap_vec'}, {'time_vec'}];
AxisNameCell = [{'mean nucleus concentration (au)'},{'AP (%)'}, {'minutes into nc14'}];
SigmaVec = [1/30 1 2];
for pv = 2:numel(PlotVariableCell)
    plot_var = PlotVariableCell{pv};
    axis_name = AxisNameCell{pv};
    var_sigma = SigmaVec(pv);
    pv_bins = linspace(min(master_struct(1).(plot_var)),max(master_struct(1).(plot_var)),30);

    for k = 1:numel(master_struct)

        spot_pt = NaN(NBoots,numel(pv_bins));
        null_pt = NaN(NBoots,numel(pv_bins));
        delta_pt = NaN(NBoots,numel(pv_bins));
        counts_pt = NaN(NBoots,numel(pv_bins));
        sample_index = 1:numel(master_struct(1).spot_protein_vec);

        for n = 1:NBoots
            samp_ids = randsample(sample_index,numel(sample_index),true);      
            var_boot = master_struct(k).(plot_var)(samp_ids);

            for i = 1:size(spot_pt,2)
                diff_vec = var_boot-pv_bins(i);
                if max(diff_vec) < 0 || min(diff_vec) > 0
                    continue
                end
                var_weights = exp(-.5*((var_boot-pv_bins(i))/var_sigma).^2);

                counts_pt(n,i) = nansum(var_weights);
                spot_pt(n,i) = nansum(master_struct(k).spot_protein_vec(samp_ids).*var_weights)/counts_pt(n,i);
                null_pt(n,i) = nansum(master_struct(k).null_protein_vec(samp_ids).*var_weights)/counts_pt(n,i);
                delta_pt(n,i) = nansum(master_struct(k).delta_protein_vec(samp_ids).*var_weights)/counts_pt(n,i);
            end
        end

        master_struct(k).counts = nanmean(counts_pt);

        master_struct(k).spot_mean = nanmean(spot_pt);
        master_struct(k).spot_ste = nanstd(spot_pt);
        master_struct(k).null_mean = nanmean(null_pt);
        master_struct(k).null_ste = nanstd(null_pt);
        master_struct(k).delta_mean = nanmean(delta_pt);
        master_struct(k).delta_ste = nanstd(delta_pt);
    end

    cm = jet(128); 
    prim_counts = master_struct(1).counts;    
    xlims = [pv_bins(find(~isnan(prim_counts),1)) pv_bins(find(~isnan(prim_counts),1,'last'))];
    pt_scale_fig = figure;
    hold on
    e = errorbar(pv_bins, master_struct(1).spot_mean, master_struct(1).spot_ste,'Color',cm(30,:),'LineWidth',1.5);
    e.CapSize = 0;
    e = errorbar(pv_bins, master_struct(1).null_mean, master_struct(1).null_ste,'Color','black','LineWidth',1.5);
    e.CapSize = 0;
    lgd_str = [{['locus (' master_struct(1).gene_name ')']}, {'control'}];
    if numel(master_struct) == 2
        e = errorbar(pv_bins, master_struct(2).spot_mean, master_struct(2).spot_ste,'Color',cm(120,:),'LineWidth',1.5);
        e.CapSize = 0;
        lgd_str = [lgd_str{:} {['locus (' master_struct(2).gene_name ')']}];
    end
    grid on
    legend(lgd_str,'Location','northwest')
    ylabel('local concentration (au)')
    xlabel(axis_name)
    xlim(xlims)
    saveas(pt_scale_fig,[FigPath 'pt_scaling' plot_var '.png'])

    pt_diff_fig = figure;
    hold on
    e = errorbar(pv_bins, master_struct(1).delta_mean, master_struct(1).delta_ste,'Color',cm(30,:),'LineWidth',1.5);
    e.CapSize = 0;
    lgd_str = [{master_struct(1).gene_name}];
    if numel(master_struct) == 2
        e = errorbar(pv_bins, master_struct(2).delta_mean, master_struct(2).delta_ste,'Color',cm(120,:),'LineWidth',1.5);
        e.CapSize = 0;
        lgd_str = [lgd_str{:} {master_struct(2).gene_name}];
    end
    xlim(xlims)
    xlabel(axis_name)
    ylabel('difference between locus and control (au)')
    legend(lgd_str,'Location','northwest')
    grid on
    saveas(pt_diff_fig,[FigPath 'delta_scaling' plot_var '.png'])
end

% if strcmpi(ProjectCell{i},'Hb_NbGFP_hbBAC_mCherry') && numel(ProjectCell) == 2
%     time_filters = {15:30,1:15};
%     for i = 1:2
%         spot_protein_snip_mean = nanmean(spot_protein_snips_mixed,3);
%         null_protein_snip_mean = nanmean(null_protein_snips_mixed,3);
% 
%         lb = prctile([spot_protein_snip_mean(:)'  null_protein_snip_mean(:)'],1);
%         ub = prctile([spot_protein_snip_mean(:)'  null_protein_snip_mean(:)'],99);
% 
%         xtick_string = "set(gca,'xtick',1:5:snip_size,'xticklabel',round(((1:5:snip_size)-round(snip_size/2))*PixelSize,2))";
%         ytick_string = "set(gca,'ytick',1:5:snip_size,'yticklabel',round(((1:5:snip_size)-round(snip_size/2))*PixelSize,2))";
% 
%         spot_pt_snip_fig = figure;
%         colormap(jet(128))
%         imagesc(spot_protein_snip_mean)
%         title([protein_name '-' protein_fluor ' at ' 'Active ' gene_name ' Locus'])
%         caxis([lb ub])
%         h = colorbar;
%         ylabel('\mum','FontSize',12)
%         xlabel('\mum','FontSize',12)
%         ylabel(h,[protein_name '-' protein_fluor ' concentration (au)'],'FontSize',12)
%         eval(xtick_string)
%         eval(ytick_string)
%         saveas(spot_pt_snip_fig,[FigPath write_string '_mean_pt_snippet_spot.png']);    
% 
%         null_pt_snip_fig = figure;
%         colormap(jet(128))
%         imagesc(null_protein_snip_mean)
%         title([protein_name '-' protein_fluor ' at Control Locus'])
%         caxis([lb ub])
%         h = colorbar;
%         ylabel('\mum','FontSize',12)
%         xlabel('\mum','FontSize',12)
%         ylabel(h,[protein_name '-' protein_fluor ' concentration (au)'],'FontSize',12)
%         eval(xtick_string)
%         eval(ytick_string)
%         saveas(null_pt_snip_fig,[FigPath write_string '_mean_pt_snippet_null.png']);  
% 
%         % Make diff image
%         rel_protein_snip_mean = (spot_protein_snip_mean) ./ null_protein_snip_mean;
%         rel_pt_snip_fig = figure;
%         colormap(jet(128))
%         imagesc(rel_protein_snip_mean)
%         title(['Relative ' protein_name '-' protein_fluor ' Enrichment at Active ' gene_name ' Locus'])
%         % caxis([lb ub])
%         h = colorbar;
%         ylabel('\mum','FontSize',12)
%         xlabel('\mum','FontSize',12)
%         ylabel(h,[protein_name '-' protein_fluor ' fold enrichment'],'FontSize',12)
%         eval(xtick_string)
%         eval(ytick_string)
%         saveas(null_pt_snip_fig,[FigPath write_string '_mean_pt_snippet_rel.png']);  
%     end
% end
