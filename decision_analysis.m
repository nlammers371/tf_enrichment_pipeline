% Script to assess time-scale (sample size) needed to distinguish reliably
% between locus and control
clear
close all
% ID variable
project = 'Dl_Venus_snaBAC_mCherry';
ControlType = 'edge';
dropboxFolder = 'E:\Nick\Dropbox (Garcia Lab)\ProcessedEnrichmentData\';
dataPath = [dropboxFolder project '\'];
writePath = ['E:\Nick\Dropbox (Garcia Lab)\LocalEnrichmentFigures\' project '\'];
mkdir(writePath);
% load data
load([dataPath 'nucleus_struct_protein.mat'])

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

% basic info params
pixelSize = nucleus_struct_protein(1).PixelSize;
ROIRadius = .3; % radus (um) of region used to query and compare TF concentrations
distLim = .6;

%% protein concentration vectors using different radii
% Generate distance vector for filtering snip stacks
dist_vec = [nucleus_struct_protein.(['spot_' ControlType '_dist_vec'])];
time_vec = [nucleus_struct_protein.time];
mf_protein_vec = [nucleus_struct_protein.protein];

set_vec = [];
particle_vec = [];
snip_filter = [];
for i = 1:numel(nucleus_struct_protein)
    snip_frame_vec = nucleus_struct_protein(i).snip_frame_vec;  
    nc_frames = nucleus_struct_protein(i).frames;
    snip_filter = [snip_filter ismember(nc_frames,snip_frame_vec)];
end    

% Snip stacks
spot_protein_snips = cat(3,nucleus_struct_protein.spot_protein_snips);
null_protein_snips = cat(3,nucleus_struct_protein.([ControlType '_null_protein_snips']));
spot_mcp_snips = cat(3,nucleus_struct_protein.spot_mcp_snips);
null_mcp_snips = cat(3,nucleus_struct_protein.([ControlType '_null_mcp_snips']));

% Make r reference array
snip_size = size(spot_protein_snips,1);
[y_ref, x_ref] = meshgrid(1:snip_size,1:snip_size);
r_ref = sqrt((x_ref-ceil(snip_size/2)).^2 + (y_ref-ceil(snip_size/2)).^2)*pixelSize;

% Make protein arrays using a range of different ROI sizes
roi_index = pixelSize:pixelSize:5*pixelSize;
spot_protein_array = NaN(size(spot_protein_snips,3),numel(roi_index));
null_protein_array = NaN(size(spot_protein_snips,3),numel(roi_index));

for i = 1:numel(roi_index)   
    r_ft = repmat(r_ref<=roi_index(i),1,1,size(spot_protein_snips,3));
    spot_protein_array(:,i) = reshape(nansum(nansum(r_ft.*spot_protein_snips,1),2),1,[]) / sum(r_ref(:)<=roi_index(i));
    null_protein_array(:,i) = reshape(nansum(nansum(r_ft.*null_protein_snips,1),2),1,[]) / sum(r_ref(:)<=roi_index(i));
end

% apply distance filter
dist_ft = dist_vec(snip_filter==1)*pixelSize >= distLim;
spot_protein_array = spot_protein_array(dist_ft,:);
null_protein_array = null_protein_array(dist_ft,:);

%% Simulate to estimate error rate as function of sample size
cm = jet(128);
n_sim = 1000;
n_draws = 100;
n_vec = 1:n_draws;
sample_index = 1:size(spot_protein_array,1);

mean_deltas = NaN(n_draws,numel(roi_index));
success_rates = NaN(n_draws,numel(roi_index));

for r = 1:numel(roi_index)
    delta_mat = NaN(n_draws,n_sim);  
    for n = 1:n_sim
        s_ids = randsample(sample_index,n_draws,false);
        delta = spot_protein_array(s_ids,r)-null_protein_array(s_ids,r);
        delta_mat(:,n) = cumsum(delta) ./ n_vec';
    end
    success_rate = nanmean(delta_mat > 0,2);
    delta_mean = nanmean(delta_mat,2);
    mean_deltas(:,r) = delta_mean;
    success_rates(:,r) = success_rate;

    decision_fig = figure;
    yyaxis left
    hold on
    plot(delta_mat,'-','Color',[cm(10,:) .02])
    p1 = plot(delta_mean,'-','Color',cm(10,:),'LineWidth',1.2);
    ylabel('locus minus control')
    ylim([-2000 5000])
    yyaxis right
    p2 = plot(success_rate,'Color','black','LineWidth',1.5);
    ylabel('success rate')
    xlim([1 n_draws])
    legend([p1 p2],'average difference','success rate')
    xlabel('number of independent samples')
    grid on
    title(['Decision Dynamics (ROI = ' num2str(round(roi_index(r),2)) ')'])
    saveas(decision_fig,[writePath 'decision_plot_roi' num2str(round(roi_index(r),2)) '.png']);
end
%%
% make plot comparing convergence rates
success_rate_fig = figure;
hold on
lgd_str = {};
for i = 1:numel(roi_index)
    plot(success_rates(:,i),'LineWidth',1.5)
    lgd_str = [lgd_str{:} {['roi=' num2str(round(roi_index(i),2))]}];   
end
grid on
legend(lgd_str{:})
xlabel('number of independent samples')
ylabel('success rate')

title('Assessing relative information content for different ROIs')
xlim([0 30])
saveas(success_rate_fig,[writePath 'success_rates_fig.png'])