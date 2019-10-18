% script to delve into details of spatio-temporal enrichment dynamics
clear
close all
addpath('utilities')
% define core ID variables
project = 'Dl-Ven_snaBAC-mCh';
% project = 'Dl-Ven_hbP2P-mCh';
dropboxFolder =  'E:\Nick\LivemRNA\Dropbox\';
dataPath = [dropboxFolder 'ProcessedEnrichmentData\' project '\'];
figPath = [dropboxFolder 'LocalEnrichmentFigures\' project '\enrichment_disaggregation\'];
mkdir(figPath)
% load data
load([dataPath 'nucleus_struct_protein.mat'])
%%
PixelSize = nucleus_struct_protein(1).PixelSize;
DistLim = .8;
dist_vec = [nucleus_struct_protein.spot_edge_dist_vec]*PixelSize;
dist_filter = dist_vec >= DistLim;
% Snip stacks
spot_protein_snips = cat(3,nucleus_struct_protein.spot_protein_snips);
null_protein_snips = cat(3,nucleus_struct_protein.edge_null_protein_snips);
% spot_mcp_snips = cat(3,nucleus_struct_protein.spot_mcp_snips);
% null_mcp_snips = cat(3,nucleus_struct_protein.edge_null_mcp_snips);
delta_protein_snips = spot_protein_snips - null_protein_snips;
delta_bins = linspace(prctile(delta_protein_snips(:),5),prctile(delta_protein_snips(:),95),101);
% Make r reference array
snip_size = size(spot_protein_snips,1);
[y_ref, x_ref] = meshgrid(1:snip_size,1:snip_size);
r_ref = sqrt((x_ref-ceil(snip_size/2)).^2 + (y_ref-ceil(snip_size/2)).^2)*PixelSize;
r_ref_rnd = round(r_ref,1);
r_index = unique(r_ref_rnd(:));
% generate radial array same size as snip stacks
delta_protein_snips = imgaussfilt(delta_protein_snips(:,:,dist_filter),.75);
r_ref_stack = repmat(r_ref_rnd,1,1,size(delta_protein_snips,3));
r_ref_vec = r_ref_stack(:);
delta_protein_vec = delta_protein_snips(:); 
% initialize profile arrays
r_delta_mat = NaN(numel(delta_bins)-1,numel(r_index));
nan_ft = ~isnan(delta_protein_vec);
for r = 1:numel(r_index)
    delta_sub_vec = delta_protein_vec(r_ref_vec==r_index(r)&nan_ft);
    ct_vec = histcounts(delta_sub_vec,delta_bins,'Normalization','probability');
    r_delta_mat(:,r) = ct_vec;
end

%% Make figure

radial_dist_fig = figure;
hm_cm = flipud(brewermap([],'RdYlBu'));
colormap(hm_cm)
imagesc(flipud(r_delta_mat));
xlabel('distance from spot center (\mu m)')
set(gca,'xtick',1:3:22,'xticklabels',r_index(1:3:22))
ylabel('Dl enrichment (au)')
set(gca,'ytick',(1:10:101),'yticklabels',fliplr(delta_bins(1:10:101)))
h = colorbar;
ylabel(h, 'share')
set(gca,'FontSize', 12);
saveas(radial_dist_fig,[figPath 'radial_enrichment_distributions.tif'])

%% find distribution over position of brightest pixels
max_ind_vec = NaN(1,size(delta_protein_snips,3));
count_array = zeros(size(delta_protein_snips(:,:,1)));
test = imgaussfilt(ones(size(r_ref)),1);

for i = 1:numel(max_ind_vec)
%     if sum(sum(~isnan(spot_protein_snips_sm(:,:,i)))) < .5*numel(spot_protein_snips_sm(:,:,i))
%         error('afsasd')
%         continue
%     end
    snip = delta_protein_snips(:,:,i);
    [~, max_ind_vec(i)] = nanmax(delta_protein_snips(:,:,i),[],[1,2],'linear');        
    count_array(max_ind_vec(i)) = count_array(max_ind_vec(i)) + 1;
end

count_array = count_array(2:end-1,2:end-1);
count_array = count_array/sum(count_array(:)) * 100;

mean_spot = nanmean(delta_protein_snips,3);
mean_spot = mean_spot(2:end-1,2:end-1);
mean_spot = mean_spot / sum(mean_spot(:)) * 100;

% make figure comparing mean fluo and count snips
enrichment_source_fit = figure('Position',[0,0,1024,512]);
colormap(hm_cm)
subplot(1,2,1)
imagesc(count_array);
h = colorbar;
ylabel(h,'percent share')
title('distribution of maxima')
% caxis([0 20e-3])


subplot(1,2,2)
imagesc(mean_spot);
h = colorbar;
title('average fluorescence')
ylabel(h,'percent share')
% caxis([0 20e-3])

saveas(enrichment_source_fit,[figPath 'maxima_position_dist_fig.tif'])