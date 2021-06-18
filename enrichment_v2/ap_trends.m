clear
close all

% designate project
projectName = 'Bcd-GFP_hbMS2-mCh_Airy_fast_int';%'Bcd-GFP_hbP2P-mCh';%
liveProject = LiveEnrichmentProject(projectName);
resultsRoot = [liveProject.dataPath filesep];
FigurePath = liveProject.figurePath;

% load data
load([resultsRoot 'spot_struct.mat'])
load([resultsRoot 'spot_struct_protein.mat'])
load([resultsRoot 'proteinSamplingInfo.mat'])
% load([resultsRoot 'snip_data.mat'])

%% calculate AP-dependent profiles for (1) nuclear concentration, (2) local
% concentration, (3) and mean spot fluorescence

% specify basic qc constraints
DistLim = 0.8;
PixelSize = liveProject.includedExperiments{1}.pixelSize_um;
nAPBins = 10;
nBcdBins = 10;
nBoots = 100;
minDP = 50;

% extract basic vectors
set_vec = [spot_struct_protein.setID];
n_vec = [spot_struct_protein.N];
set_filter = set_vec < 3 & [spot_struct_protein.TraceQCFlag];%true(size(set_vec));
spot_edge_dist_vec = [spot_struct_protein(set_filter).spot_edge_dist_vec]*PixelSize;
spot_edge_dist_flags = spot_edge_dist_vec>=DistLim;

x_pos_vec = [spot_struct_protein(set_filter).xPosParticle];
y_pos_vec = [spot_struct_protein(set_filter).yPosParticle];
ap_pos_vec = [spot_struct_protein(set_filter).APPosParticle];
spot_protein_vec = [spot_struct_protein(set_filter).spot_protein_vec];
ctrl_protein_vec = [spot_struct_protein(set_filter).edge_null_protein_vec];
nucleus_protein_vec = [spot_struct_protein(set_filter).nuclear_protein_vec];
% nucleus_protein_vec_alt = [spot_struct_protein(set_filter).nuclear_protein_vec];
fluo_vec = [spot_struct_protein(set_filter).fluo];
time_vec = [spot_struct_protein(set_filter).time]/60;

% apply edge filtering
ap_pos_vec_ft = ap_pos_vec(spot_edge_dist_flags);
x_pos_vec_ft = x_pos_vec(spot_edge_dist_flags);
y_pos_vec_ft = y_pos_vec(spot_edge_dist_flags);
spot_protein_vec_ft = spot_protein_vec(spot_edge_dist_flags);
ctrl_protein_vec_ft = ctrl_protein_vec(spot_edge_dist_flags);
nucleus_protein_vec_ft = nucleus_protein_vec(spot_edge_dist_flags);
fluo_vec_ft = fluo_vec(spot_edge_dist_flags);
time_vec_ft = time_vec(spot_edge_dist_flags);

% calculate bounds for AP-dependent calculations
minAP = 18;%prctile(ap_pos_vec_ft,0.5);
maxAP = 44;%prctile(ap_pos_vec_ft,99.5);
minBcd = prctile(nucleus_protein_vec_ft,0.5);
maxBcd = prctile(nucleus_protein_vec_ft,99.5);
maxTime = 30;
minTime = 5;
dT = 5;

% generate grouping vector
ap_bin_vec = linspace(minAP,maxAP,nAPBins+1);
bcd_bin_vec = linspace(minBcd,maxBcd,nBcdBins+1);
time_bin_vec = linspace(minTime,maxTime,5);
nTimeBins = length(time_bin_vec);

% make plot axes
ap_plot = ap_bin_vec(1:end-1) + diff(ap_bin_vec)/2;
bcd_plot = bcd_bin_vec(1:end-1) + diff(bcd_bin_vec)/2;
time_plot = round(time_bin_vec(1:end-1) + diff(time_bin_vec)/2,1);

% initialize bootstrap arrays
nucleus_protein_array = NaN(nTimeBins,nAPBins,nBoots);

delta_rel_protein_array_ap = NaN(nTimeBins,nAPBins,nBoots);
ctrl_protein_array_ap = NaN(nTimeBins,nAPBins,nBoots);
delta_protein_array_ap = NaN(nTimeBins,nAPBins,nBoots);
fluo_array_ap = NaN(nTimeBins,nAPBins,nBoots);

delta_rel_protein_array_bcd = NaN(nTimeBins,nAPBins,nBoots);
ctrl_protein_array_bcd = NaN(nTimeBins,nAPBins,nBoots);
delta_protein_array_bcd = NaN(nTimeBins,nAPBins,nBoots);
fluo_array_bcd = NaN(nTimeBins,nAPBins,nBoots);

for t = 1:length(time_bin_vec)-1
    for a = 1:length(ap_bin_vec)-1
        % generate AP filter
        ap_filter = ap_pos_vec_ft>=ap_bin_vec(a) & ap_pos_vec_ft<ap_bin_vec(a+1);
        bcd_filter = nucleus_protein_vec_ft>=bcd_bin_vec(a) & nucleus_protein_vec_ft<bcd_bin_vec(a+1);
        time_filter = time_vec_ft>=time_bin_vec(t) & time_vec_ft<time_bin_vec(t+1);
        ap_time_filter =  ap_filter & time_filter;
        bcd_time_filter =  bcd_filter & time_filter;
        
        
        % generate indexin g vec
        ap_indices = find(ap_time_filter);
        bcd_indices = find(bcd_time_filter);
        for n = 1:nBoots
            if sum(ap_time_filter) > 25
                % randomly select bootstrap indices
                ap_boot_indices = randsample(ap_indices,length(ap_indices),true);
                % now, calculate mean value for each parameter of interest
                nucleus_protein_array(t,a,n) = nanmean(nucleus_protein_vec_ft(ap_boot_indices));   

    %                 delta_rel_protein_array(t,a,n) = nanmean((spot_protein_vec_ft(boot_indices)-ctrl_protein_vec_ft(boot_indices))./...
    %                                               ctrl_protein_vec_ft(boot_indices));

                ctrl_protein_array_ap(t,a,n) = nanmean(ctrl_protein_vec_ft(ap_boot_indices));
                delta_protein_array_ap(t,a,n) = nanmean(spot_protein_vec_ft(ap_boot_indices)-ctrl_protein_vec_ft(ap_boot_indices));
                delta_rel_protein_array_ap(t,a,n) = delta_protein_array_ap(t,a,n)./ctrl_protein_array_ap(t,a,n);
                fluo_array_ap(t,a,n) = nanmean(fluo_vec_ft(ap_boot_indices));
            end
            
            if sum(bcd_time_filter) > 25
                % randomly select bootstrap indices
                bcd_boot_indices = randsample(bcd_indices,length(bcd_indices),true);
                % now, calculate mean value for each parameter of interest                

                ctrl_protein_array_bcd(t,a,n) = nanmean(ctrl_protein_vec_ft(bcd_boot_indices));
                delta_protein_array_bcd(t,a,n) = nanmean(spot_protein_vec_ft(bcd_boot_indices)-ctrl_protein_vec_ft(bcd_boot_indices));
                delta_rel_protein_array_bcd(t,a,n) = delta_protein_array_ap(t,a,n)./ctrl_protein_array_bcd(t,a,n);
                fluo_array_bcd(t,a,n) = nanmean(fluo_vec_ft(bcd_boot_indices));
            end
        end
    end
end

% plot ap trends
MarkerSize = 75;

% close all
c_temp = linspace(0,1,length(time_bin_vec));
c_ticks = c_temp(1:end-1) + diff(c_temp)/2;%c_ticks;
%%%%%%%%%%%%%%%%%%%%%
%%% Spot intensity
%%%%%%%%%%%%%%%%%%%%%
mean_fluo_array = nanmean(fluo_array_ap,3);
ste_fluo_array = nanstd(fluo_array_ap,[],3);

fluo_fig = figure;
hold on
cm1 = flipud(brewermap(length(time_plot),'Spectral'));
colormap(cm1);
for t = 1:length(time_plot)
    errorbar(ap_plot,mean_fluo_array(t,:),ste_fluo_array(t,:),'CapSize',0,'Color','k','LineWidth',1)
    scatter(ap_plot,mean_fluo_array(t,:),MarkerSize,'MarkerFaceColor',cm1(t,:),'MArkerEdgeColor','k')
end
% formating
set(gca,'Fontsize',14);
set(gca,'Color',[228,221,209]/255) 
grid on
fluo_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
% make colorbar
h = colorbar;
h.Ticks = c_ticks;
h.TickLabels = time_plot;
% axis labels
xlabel('AP position (% embryo length)')
ylabel('hbP2P spot intensity (au)')
ylabel(h,'minutes into nc14')
% axis limits
xlim([ap_plot(1)-1 ap_plot(end)+1])
% save
saveas(fluo_fig,[FigurePath 'ap_fluo_trends.png'])

%%%%%%%%%%%%%%%%%%%%%
%%% Bcd enrichment
%%%%%%%%%%%%%%%%%%%%%
mean_delta_array = nanmean(delta_protein_array_ap,3);
ste_delta_array = nanstd(delta_protein_array_ap,[],3);

enrichment_fig = figure;
hold on
colormap(cm1);
for t = 1:length(time_plot)
    errorbar(ap_plot, mean_delta_array(t,:), ste_delta_array(t,:),'CapSize',0,'Color','k','LineWidth',0.75)
    scatter(ap_plot, mean_delta_array(t,:),MarkerSize,'MarkerFaceColor',cm1(t,:),'MarkerEdgeColor','k')
end
% plot reference line at 0
plot([0 ap_plot 100], zeros(1,length(ap_plot)+2),'--k');

% formating
set(gca,'Fontsize',14);
set(gca,'Color',[228,221,209]/255) 
grid on
enrichment_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
% make colorbar
h = colorbar;
h.Ticks = c_ticks;
h.TickLabels = time_plot;
% axis labels
xlabel('AP position (% embryo length)')
ylabel('Bcd-GFP enrichment (au)')
ylabel(h,'minutes into nc14')
% axis limits
ylim([-30 75])
xlim([ap_plot(1)-1 ap_plot(end)+1])
% save
saveas(enrichment_fig,[FigurePath 'ap_enrichment_trends.png'])

%%%%%%%%%%%%%%%%%%%%%
%%% "nuclear" bcd
%%%%%%%%%%%%%%%%%%%%%
mean_nucleus_array = nanmean(nucleus_protein_array,3);
ste_nucleus_array = nanstd(nucleus_protein_array,[],3);

nucleus_fig = figure;
hold on
colormap(cm1);
for t = 1:length(time_plot)
    errorbar(ap_plot, mean_nucleus_array(t,:), ste_nucleus_array(t,:),'CapSize',0,'Color','k','LineWidth',0.75)
    scatter(ap_plot, mean_nucleus_array(t,:),MarkerSize,'MarkerFaceColor',cm1(t,:),'MarkerEdgeColor','k')
end

% formating
set(gca,'Fontsize',14);
set(gca,'Color',[228,221,209]/255) 
grid on
nucleus_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
% make colorbar
h = colorbar;
h.Ticks = c_ticks;
h.TickLabels = time_plot;
% axis labels
xlabel('AP position (% embryo length)')
ylabel('Bcd-GFP in nucleus (au)')
ylabel(h,'minutes into nc14')
% axis limits
xlim([ap_plot(1)-1 ap_plot(end)+1])
% save
saveas(nucleus_fig,[FigurePath 'ap_nucleus_concentration_trends.png'])

%%%%%%%%%%%%%%%%%%%%%
%%% Bcd at spot
%%%%%%%%%%%%%%%%%%%%%
mean_rel_array = nanmean(delta_rel_protein_array_ap,3);
ste_rel_array = nanstd(delta_rel_protein_array_ap,[],3);

spot_fig = figure;
hold on
colormap(cm1);
for t = 1:length(time_plot)
    errorbar(ap_plot, mean_rel_array(t,:), ste_rel_array(t,:),'CapSize',0,'Color','k','LineWidth',0.75)
    scatter(ap_plot, mean_rel_array(t,:),MarkerSize,'MarkerFaceColor',cm1(t,:),'MarkerEdgeColor','k')
end

% formating
set(gca,'Fontsize',14);
set(gca,'Color',[228,221,209]/255) 
grid on
spot_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
% make colorbar
h = colorbar;
h.Ticks = c_ticks;
h.TickLabels = time_plot;
% axis labels
xlabel('AP position (% embryo length)')
ylabel('relative Bcd enrichment(au)')
ylabel(h,'minutes into nc14')
% axis limits
xlim([ap_plot(1)-1 ap_plot(end)+1])
% save
saveas(spot_fig,[FigurePath 'ap_relative_enrichment_trends.png'])

%%%%%%%%%%%%%%%%%%
%%% Make comparison plot
%%%%%%%%%%%%%%%%%%
t_index = 2;

comparison_fig = figure;
hold on
cm2 = brewermap([],'Set2');

% plot enrichment
errorbar(ap_plot, mean_delta_array(t_index,:), ste_delta_array(t_index,:),'CapSize',0,'Color','k','LineWidth',1)
scatter(ap_plot, mean_delta_array(t_index,:),MarkerSize,'MarkerFaceColor',cm2(2,:),'MarkerEdgeColor','k')

ylabel('Bcd-GFP enrichment at spot (au)')
ax = gca;
% ax.YColor = cm2(2,:);

% now plot spot fluorescence
yyaxis right
errorbar(ap_plot, mean_fluo_array(t_index,:), ste_fluo_array(t_index,:),'CapSize',0,'Color','k','LineWidth',1)
scatter(ap_plot, mean_fluo_array(t_index,:),MarkerSize,'MarkerFaceColor',cm2(3,:),'MarkerEdgeColor','k')

ax = gca;
ax.YAxis(1).Color = cm2(2,:);
ax.YAxis(2).Color = cm2(3,:);

% formating
set(gca,'Fontsize',14);
set(gca,'Color',[228,221,209]/255) 
grid on
comparison_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

% axis labels
xlabel('AP position (% embryo length)')
ylabel('hbP2P spot intensity (au)')

% axis limits
xlim([ap_plot(1)-1 ap_plot(end)+1])

% save
saveas(comparison_fig,[FigurePath 'comparison_plot.png'])


% plot bcd trends

close all

%%%%%%%%%%%%%%%%%%%%%
%%% Spot intensity
%%%%%%%%%%%%%%%%%%%%%

mean_fluo_array_bcd = nanmean(fluo_array_bcd,3);
ste_fluo_array_bcd = nanstd(fluo_array_bcd,[],3);

fluo_fig = figure;
hold on
cm1 = flipud(brewermap(length(time_plot),'Spectral'));
colormap(cm1);
for t = 1:length(time_plot)
    errorbar(bcd_plot,mean_fluo_array_bcd(t,:),ste_fluo_array_bcd(t,:),'CapSize',0,'Color','k','LineWidth',1)
    scatter(bcd_plot,mean_fluo_array_bcd(t,:),MarkerSize,'MarkerFaceColor',cm1(t,:),'MArkerEdgeColor','k')
end
% formating
set(gca,'Fontsize',14);
set(gca,'Color',[228,221,209]/255) 
grid on
fluo_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
% make colorbar
h = colorbar;
h.Ticks = c_ticks;
h.TickLabels = time_plot;
% axis labels
xlabel('Bcd-GFP in nucleus (au)')
ylabel('hbP2P spot intensity (au)')
ylabel(h,'minutes into nc14')
% axis limits
xlim([bcd_plot(1)-25 bcd_plot(end)+25])
% save
saveas(fluo_fig,[FigurePath 'bcd_fluo_trends.png'])

%%%%%%%%%%%%%%%%%%%%%
% Bcd enrichment
%%%%%%%%%%%%%%%%%%%%%
mean_delta_array_bcd = nanmean(delta_protein_array_bcd,3);
ste_delta_array_bcd = nanstd(delta_protein_array_bcd,[],3);

enrichment_fig = figure;
hold on
colormap(cm1);
for t = 1:length(time_plot)
    errorbar(bcd_plot, mean_delta_array_bcd(t,:), ste_delta_array_bcd(t,:),'CapSize',0,'Color','k','LineWidth',0.75)
    scatter(bcd_plot, mean_delta_array_bcd(t,:),MarkerSize,'MarkerFaceColor',cm1(t,:),'MarkerEdgeColor','k')
end
% plot reference line at 0
plot([0 bcd_plot bcd_plot(end)+100], zeros(1,length(ap_plot)+2),'--k');

% formating
set(gca,'Fontsize',14);
set(gca,'Color',[228,221,209]/255) 
grid on
enrichment_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
% make colorbar
h = colorbar;
h.Ticks = c_ticks;
h.TickLabels = time_plot;
% axis labels
xlabel('Bcd-GFP in nucleus (au)')
ylabel('Bcd-GFP enrichment at spot (au)')
ylabel(h,'minutes into nc14')
% axis limits
% ylim([-30 75])
xlim([bcd_plot(1)-25 bcd_plot(end)+25])
% save
saveas(enrichment_fig,[FigurePath 'bcd_enrichment_trends.png'])


%%%%%%%%%%%%%%%%%%%%%
% Bcd at spot
%%%%%%%%%%%%%%%%%%%%%
mean_rel_array_bcd = nanmean(delta_rel_protein_array_bcd,3);
ste_rel_array_bcd = nanstd(delta_rel_protein_array_bcd,[],3);

spot_fig = figure;
hold on
colormap(cm1);
for t = 1:length(time_plot)
    errorbar(bcd_plot, mean_rel_array_bcd(t,:), ste_rel_array_bcd(t,:),'CapSize',0,'Color','k','LineWidth',0.75)
    scatter(bcd_plot, mean_rel_array_bcd(t,:),MarkerSize,'MarkerFaceColor',cm1(t,:),'MarkerEdgeColor','k')
end

% formating
set(gca,'Fontsize',14);
set(gca,'Color',[228,221,209]/255) 
grid on
spot_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
% make colorbar
h = colorbar;
h.Ticks = c_ticks;
h.TickLabels = time_plot;
% axis labels
xlabel('Bcd-GFP in nucleus (au)')
ylabel('relative Bcd enrichment(au)')
ylabel(h,'minutes into nc14')
% axis limits
xlim([bcd_plot(1)-25 bcd_plot(end)+25])
% save
saveas(spot_fig,[FigurePath 'bcd_rel_enrichment_trends.png'])


%%
plot_index = 453;
qc_indices = find(set_filter);
close all

figure;
hold on
plot(spot_struct_protein(qc_indices(plot_index)).time,spot_struct_protein(qc_indices(plot_index)).fluo,'-o')
yyaxis right
pt_vec = spot_struct_protein(qc_indices(plot_index)).spot_protein_vec;
plot(spot_struct_protein(qc_indices(plot_index)).time,pt_vec,'-o')
% plot(spot_struct_protein(qc_indices(plot_index)).time,spot_struct_protein(qc_indices(plot_index)).serial_null_protein_vec,'--s','Color','y')


