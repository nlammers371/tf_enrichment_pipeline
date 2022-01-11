clear
close all


dataRoot = 'C:\Users\nlamm\Dropbox (Personal)\ProcessedEnrichmentData\';
FigurePath = 'C:\Users\nlamm\Dropbox (Personal)\LocalEnrichmentFigures\PipelineOutput\Bcd-GFP_at_hbP2P_comparisons\';
if ~exist(dataRoot,'dir')
    dataRoot = 'S:\Nick\Dropbox\ProcessedEnrichmentData\';
    FigurePath = 'S:\Nick\Dropbox\LocalEnrichmentFigures\PipelineOutput\Bcd-GFP_at_hbP2P_comparisons\';
end    
mkdir(FigurePath);

% designate projects to examine
projectNameCell = {'Bcd-GFP_hbMS2-mCh_AiryHP','Bcd-GFP-McpMcherry-hbP2P-delta6'};%'Bcd-GFP_hbP2P-mCh';%
                 
master_struct = struct;                 
% load all the relecant projects
for p = 1:length(projectNameCell)
    projectName = projectNameCell{p};    
%     liveProject = LiveEnrichmentProject(projectName);
    resultsRoot = [dataRoot filesep projectName filesep];    

    % load data
    load([resultsRoot 'spot_struct.mat'])
    load([resultsRoot 'spot_struct_protein.mat'])
    load([resultsRoot 'proteinSamplingInfo.mat'])
    load([resultsRoot 'snip_data.mat'])
    
    master_struct(p).spot_struct = spot_struct;
    master_struct(p).spot_struct_protein = spot_struct_protein;
    master_struct(p).proteinSamplingInfo = proteinSamplingInfo;
    master_struct(p).snip_data = snip_data;
    
    % get project info
    projectName = projectNameCell{p}; 
    liveProject = LiveEnrichmentProject(projectName);
    master_struct(p).PixelSize = liveProject.includedExperiments{1}.pixelSize_um;
    master_struct(p).liveProject = liveProject;
end    
% load([resultsRoot 'snip_data.mat'])

%% calculate average enrichment for constrained window of time and space
sets_to_use_cell = {[1:4],[1 2 3 4]};
norm_factor_cell = {[1 1 4/3 4/3] [1 4/3 1 4/3]};

% specify basic qc constraints
DistLim = 0.8; % minimum distance to nuclear boundary we'll permit
nAPBins = 15;
nTimeBins = 65;
ap_bin_vec = linspace(0.1,0.5,nAPBins+1);
time_bin_vec = linspace(-20,45,nTimeBins+1);
% time_bin_vec = -20:5:45;
% nTimeBins = length(time_bin_vec)-1;
nBoots = 100;
minDP = 50;
filtWidth = 3;
filtSigma = 2;
imageFilter=fspecial('gaussian',filtWidth,filtSigma);

ap_window = [0.2 0.3];
nc13_window = [-15 -5];
nc14_window = [10 20];

try
    parpool(24);
catch
    % do nothing
end

fluo_mean_vec = NaN(2,length(master_struct));
fluo_ste_vec = NaN(2,length(master_struct));
delta_mean_vec = NaN(2,length(master_struct));
delta_ste_vec = NaN(2,length(master_struct));
nuc_mean_vec = NaN(2,length(master_struct));
nuc_ste_vec = NaN(2,length(master_struct));


fluo_mean_array = NaN(nTimeBins, nAPBins, length(master_struct));
fluo_ste_array = NaN(nTimeBins, nAPBins, length(master_struct));

fluo_mean_array_full = NaN(nTimeBins, nAPBins, length(master_struct));
fluo_ste_array_full = NaN(nTimeBins, nAPBins, length(master_struct));

nuc_bcd_mean_array = NaN(nTimeBins, nAPBins, length(master_struct));
nuc_bcd_ste_array = NaN(nTimeBins, nAPBins, length(master_struct));

nuc_bcd_mean_array_full = NaN(nTimeBins, nAPBins, length(master_struct));
nuc_bcd_ste_array_full = NaN(nTimeBins, nAPBins, length(master_struct));

ctrl_bcd_mean_array = NaN(nTimeBins, nAPBins, length(master_struct));
ctrl_bcd_ste_array = NaN(nTimeBins, nAPBins, length(master_struct));

spot_bcd_mean_array = NaN(nTimeBins, nAPBins, length(master_struct));
spot_bcd_ste_array = NaN(nTimeBins, nAPBins, length(master_struct));

delta_bcd_mean_array = NaN(nTimeBins, nAPBins, length(master_struct));
delta_bcd_ste_array = NaN(nTimeBins, nAPBins, length(master_struct));

% Loop through each project and extract key data vectors
for  p = 1:length(master_struct)
    spot_struct_protein = master_struct(p).spot_struct_protein;
%     set_vec = [spot_struct_protein.setID];
    traceQCVec = [spot_struct_protein.N]>=7;
    spot_struct = master_struct(p).spot_struct;
    particleIDVec = [spot_struct.particleID];
    
    % "raw" vectors
    spot_protein_vec = [spot_struct_protein(traceQCVec).spot_protein_vec];
    ctrl_protein_vec = [spot_struct_protein(traceQCVec).edge_null_protein_vec];
    delta_protein_vec = spot_protein_vec-ctrl_protein_vec;
    bcd_vec = [spot_struct_protein(traceQCVec).nuclear_protein_vec];
    fluo_vec = [spot_struct_protein(traceQCVec).fluo];
    time_vec = [spot_struct_protein(traceQCVec).time];
    ap_vec = [spot_struct_protein(traceQCVec).APPosParticle]/100;
    
    % full vectors
    bcd_vec_full = [spot_struct.rawNCProteinInterp];
    fluo_vec_full = [spot_struct.fluoInterp];
%     fluo_vec_full(isnan(fluo_vec_full)) = 0;
    ap_vec_full = [spot_struct.APPosNucleusInterp];
    time_vec_full = [spot_struct.timeInterp];
      
    % generate id vector to allow us to map back to individual spots
    set_vec = [spot_struct.setID];
    nc_vec = [spot_struct.nc];
    set_index = unique(set_vec);
    zero_time_index = NaN(size(set_index));
    offset_norm_index = norm_factor_cell{p};
    for s = sets_to_use_cell{p}
        time_temp = [spot_struct(nc_vec==14&set_vec==set_index(s)).time];
        ap_temp = [spot_struct(nc_vec==14&set_vec==set_index(s)).APPosNucleus];
        offset_temp = [spot_struct(nc_vec==14&set_vec==set_index(s)).fluo];
        zero_time_index(s) = nanmin(time_temp);
%         time_round = (time_temp-zero_time_index(s))/60;
%         ft = time_round>10&time_round<15 &ap_temp>=0.3&ap_temp<=0.35;
%         offset_norm_index(s) = nanmean(offset_temp(ft));
    end
    master_struct(p).fluo_mean = offset_norm_index;
    zero_time_vec = [];    
    fluo_norm_vec = [];    
    for i = find(traceQCVec)
        setID = spot_struct_protein(i).setID;
        zero_time_vec = [zero_time_vec repelem(zero_time_index(set_index==setID),length(spot_struct_protein(i).time))];
        fluo_norm_vec = [fluo_norm_vec repelem(offset_norm_index(set_index==setID),length(spot_struct_protein(i).time))];
    end
    
    zero_time_full_vec = [];    
    fluo_norm_full_vec = [];
    for i = 1:length(spot_struct)
        setID = spot_struct(i).setID;
        zero_time_full_vec = [zero_time_full_vec repelem(zero_time_index(set_index==setID),length(spot_struct(i).timeInterp))];
        fluo_norm_full_vec = [fluo_norm_full_vec repelem(offset_norm_index(set_index==setID),length(spot_struct(i).timeInterp))];
    end  
    
    % apply corrections
    time_vec_full = time_vec_full - zero_time_full_vec;
    time_vec = time_vec - zero_time_vec;
    
    fluo_vec_full = fluo_vec_full .* fluo_norm_full_vec;
    fluo_vec = fluo_vec .* fluo_norm_vec;                    
    
    fluo_temp = NaN(size(fluo_mean_array_full,1),size(fluo_mean_array_full,2),nBoots);
    bcd_temp = NaN(size(fluo_mean_array_full,1),size(fluo_mean_array_full,2),nBoots);    
    parfor n = 1:nBoots
        boot_indices = randsample(1:length(time_vec_full),length(time_vec_full),true);
        [N,~,~, timeID_full, apID_full] = histcounts2(time_vec_full(boot_indices)/60, ap_vec_full(boot_indices), time_bin_vec, ap_bin_vec);
        keep_indices = find(apID_full~=0&timeID_full~=0);
        lin_ind_full = sub2ind(size(N),timeID_full(keep_indices),apID_full(keep_indices));
        lin_index = unique(lin_ind_full);
        
        % calculate group means
        temp_fluo = fluo_temp(:,:,n);        
        temp_fluo(lin_index) = grpstats(fluo_vec_full(boot_indices(keep_indices)),lin_ind_full,'mean');                
        fluo_temp(:,:,n) = nanconv(temp_fluo,imageFilter,'nanout');
        
        temp_bcd = bcd_temp(:,:,n);
        temp_bcd(lin_index) = grpstats(bcd_vec_full(boot_indices(keep_indices)),lin_ind_full,'mean');        
        bcd_temp(:,:,n) = nanconv(temp_bcd,imageFilter,'nanout');                  
        
    end
        
    % record
    fluo_mean_array_full(:,:,p) = nanmean(fluo_temp,3);
    fluo_ste_array_full(:,:,p) = nanstd(fluo_temp,[],3);
    nuc_bcd_mean_array_full(:,:,p) = nanmean(bcd_temp,3);    
    nuc_bcd_ste_array_full(:,:,p) = nanstd(bcd_temp,[],3);    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % now filtered quantities
    % check distance from nucleus boundary
    spot_edge_dist_vec = [spot_struct_protein(traceQCVec).spot_edge_dist_vec]*master_struct(p).PixelSize;
    se_indices = find(spot_edge_dist_vec>=DistLim);
          
    fluo_temp = NaN(size(fluo_mean_array,1),size(fluo_mean_array,2),nBoots);
    bcd_temp = NaN(size(fluo_mean_array,1),size(fluo_mean_array,2),nBoots);    
    spot_temp = NaN(size(fluo_mean_array,1),size(fluo_mean_array,2),nBoots);    
    ctrl_temp = NaN(size(fluo_mean_array,1),size(fluo_mean_array,2),nBoots);    
    delta_temp = NaN(size(fluo_mean_array,1),size(fluo_mean_array,2),nBoots);    
    
    fluo_coarse_temp = NaN(nBoots,2);
    delta_coarse_temp = NaN(nBoots,2);
    nuc_coarse_temp = NaN(nBoots,2);
    
    parfor n = 1:nBoots
        boot_indices = randsample(se_indices,length(se_indices),true);
        [N,~,~, timeID, apID] = histcounts2(time_vec(boot_indices)/60, ap_vec(boot_indices), time_bin_vec, ap_bin_vec);
        keep_indices = find(apID~=0&timeID~=0);
        lin_ind = sub2ind(size(N),timeID(keep_indices),apID(keep_indices));
        lin_index = unique(lin_ind);
        
        % calculate group means
        temp_fluo = fluo_temp(:,:,n);
        temp_fluo(lin_index) = grpstats(fluo_vec(boot_indices(keep_indices)),lin_ind,'mean');        
        fluo_temp(:,:,n) = nanconv(temp_fluo,imageFilter,'nanout');
        
        temp_bcd = bcd_temp(:,:,n);
        temp_bcd(lin_index) = grpstats(bcd_vec(boot_indices(keep_indices)),lin_ind,'mean');        
        bcd_temp(:,:,n) = nanconv(temp_bcd,imageFilter,'nanout');
        
        temp_spot = spot_temp(:,:,n);
        temp_spot(lin_index) = grpstats(spot_protein_vec(boot_indices(keep_indices)),lin_ind,'mean');        
        spot_temp(:,:,n) = nanconv(temp_spot,imageFilter,'nanout');
        
        temp_ctrl = ctrl_temp(:,:,n);
        temp_ctrl(lin_index) = grpstats(ctrl_protein_vec(boot_indices(keep_indices)),lin_ind,'mean');        
        ctrl_temp(:,:,n) = nanconv(temp_ctrl,imageFilter,'nanout');
        
        delta_temp(:,:,n) = nanconv(temp_spot-temp_ctrl,imageFilter,'nanout');
        
        % calculate cours-grained quantities
        ap_filter = ap_vec(boot_indices)>= ap_window(1) & ap_vec(boot_indices) <= ap_window(2);
        nc13_filter = time_vec(boot_indices)/60>=nc13_window(1) & time_vec(boot_indices)/60<=nc13_window(2);
        nc14_filter = time_vec(boot_indices)/60>=nc14_window(1) & time_vec(boot_indices)/60<=nc14_window(2);
        
        fluo_coarse_temp(n,:) = [nanmean(fluo_vec(boot_indices(ap_filter&nc13_filter))) nanmean(fluo_vec(boot_indices(ap_filter&nc14_filter)))];
        delta_coarse_temp(n,:) = [nanmean(delta_protein_vec(boot_indices(ap_filter&nc13_filter))) nanmean(delta_protein_vec(boot_indices(ap_filter&nc14_filter)))];
        nuc_coarse_temp(n,:) = [nanmean(bcd_vec(boot_indices(ap_filter&nc13_filter))) nanmean(bcd_vec(boot_indices(ap_filter&nc14_filter)))];
    end
        
    % record
    fluo_mean_array(:,:,p) = nanmean(fluo_temp,3);
    fluo_ste_array(:,:,p) = nanstd(fluo_temp,[],3);
    nuc_bcd_mean_array(:,:,p) = nanmean(bcd_temp,3);    
    nuc_bcd_ste_array(:,:,p) = nanstd(bcd_temp,[],3);   
    spot_bcd_mean_array(:,:,p) = nanmean(spot_temp,3);
    spot_bcd_ste_array(:,:,p) = nanstd(spot_temp,[],3);
    ctrl_bcd_mean_array(:,:,p) = nanmean(ctrl_temp,3);    
    ctrl_bcd_ste_array(:,:,p) = nanstd(ctrl_temp,[],3);   
    delta_bcd_mean_array(:,:,p) = nanmean(delta_temp,3);
    delta_bcd_ste_array(:,:,p) = nanstd(delta_temp,[],3);
            
    fluo_mean_vec(:,p) = nanmean(fluo_coarse_temp)';
    fluo_ste_vec(:,p) = nanstd(fluo_coarse_temp)';
    delta_mean_vec(:,p) = nanmean(delta_coarse_temp)';
    delta_ste_vec(:,p) = nanstd(delta_coarse_temp)';
    nuc_mean_vec(:,p) = nanmean(nuc_coarse_temp)';
    nuc_ste_vec(:,p) = nanstd(nuc_coarse_temp)';
    
end

%%  Make time trend plots
markerSize = 75;
close all
plot_ap_indices = 5:1:9;
time_axis = time_bin_vec(1:end-1) + diff(time_bin_vec)/2;
ap_axis = round(100*(ap_bin_vec(1:end-1) + diff(ap_bin_vec)/2));

for a = 1:length(plot_ap_indices)
    delta_fig = figure;
    s = [];
    cmap = brewermap([],'Set2');
    hold on 
    errorbar(time_axis,delta_bcd_mean_array(:,plot_ap_indices(a),1),delta_bcd_ste_array(:,plot_ap_indices(a),1),'-','Color','k','Capsize',0)
    s(1) = scatter(time_axis,delta_bcd_mean_array(:,plot_ap_indices(a),1),markerSize','MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k');
    
    errorbar(time_axis,delta_bcd_mean_array(:,plot_ap_indices(a),2),delta_bcd_ste_array(:,plot_ap_indices(a),2),'-','Color','k','Capsize',0)
    s(2) = scatter(time_axis,delta_bcd_mean_array(:,plot_ap_indices(a),2),markerSize','MarkerFaceColor',cmap(5,:),'MarkerEdgeColor','k');
    
    plot(time_axis,repelem(0,length(time_axis)),'--k')
    
    legend(s,'WT','\Delta 6','Location','northeast')
    
    xlabel('minutes from 13th mitosis')
    ylabel('excess Bcd at locus (au)')
    set(gca,'Fontsize',14)
    set(gca,'Color',[228,221,209]/255) 
    grid on
    delta_fig.InvertHardcopy = 'off';
    set(gcf,'color','w');
    xlim([-15 25])
    saveas(delta_fig,[FigurePath 'delta_vs_time_ap' num2str(ap_axis(a)) '.png'])
    
    
    fluo_fig = figure;
    s = [];
    cmap = brewermap([],'Set2');
    hold on 
    errorbar(time_axis,fluo_mean_array(:,plot_ap_indices(a),1),fluo_ste_array(:,plot_ap_indices(a),1),'-','Color','k','Capsize',0)
    s(1) = scatter(time_axis,fluo_mean_array(:,plot_ap_indices(a),1),markerSize,'s','MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k');
    
    errorbar(time_axis,fluo_mean_array(:,plot_ap_indices(a),2),fluo_ste_array(:,plot_ap_indices(a),2),'-','Color','k','Capsize',0)
    s(2) = scatter(time_axis,fluo_mean_array(:,plot_ap_indices(a),2),markerSize,'s','MarkerFaceColor',cmap(5,:),'MarkerEdgeColor','k');
    
    legend(s,'WT','\Delta 6','Location','northeast')
    
    xlabel('minutes from 13th mitosis')
    ylabel('average spot fluorescence (au)')
    set(gca,'Fontsize',14)
    set(gca,'Color',[228,221,209]/255) 
    grid on
    fluo_fig.InvertHardcopy = 'off';
    set(gcf,'color','w');
    xlim([-15 25])
    saveas(fluo_fig,[FigurePath 'fluo_vs_time_ap' num2str(ap_axis(a)) '.png'])
    
end    
    
 
%% make AP trend plots

plot_time_indices = find(ismember(round(time_axis),[-6 12]));
close all
for a = 1%1:length(plot_time_indices)
    delta_fig = figure;
    s = [];
    cmap = brewermap(3,'Spectral');
    hold on 
    errorbar(ap_axis,delta_bcd_mean_array(plot_time_indices(1),:,1),delta_bcd_ste_array(plot_time_indices(1),:,1),'-','Color','k','Capsize',0)
    s(1) = scatter(ap_axis,delta_bcd_mean_array(plot_time_indices(1),:,1),markerSize','MarkerFaceColor',cmap(1,:),'MarkerEdgeColor','k');
    
    errorbar(ap_axis,delta_bcd_mean_array(plot_time_indices(2),:,1),delta_bcd_ste_array(plot_time_indices(2),:,1),'-','Color','k','Capsize',0)
    s(2) = scatter(ap_axis,delta_bcd_mean_array(plot_time_indices(2),:,1),markerSize','MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k');
    
    plot(ap_axis,repelem(0,length(ap_axis)),'--k')
    
    legend(s,'WT (nc13)','WT (nc14)','Location','southwest')
    
    xlabel('AP position (% embryo length)')
    ylabel('excess Bcd at locus (au)')
    set(gca,'Fontsize',14)
    set(gca,'Color',[228,221,209]/255) 
    grid on
    delta_fig.InvertHardcopy = 'off';
    set(gcf,'color','w');
    xlim([20 50])
    saveas(delta_fig,[FigurePath 'delta_vs_ap.png'])
    
    
    fluo_fig = figure;
    s = [];
    cmap = brewermap(3,'Spectral');
    hold on 
    errorbar(ap_axis,fluo_mean_array(plot_time_indices(1),:,1),fluo_ste_array(plot_time_indices(1),:,1),'-','Color','k','Capsize',0)
    s(1) = scatter(ap_axis,fluo_mean_array(plot_time_indices(1),:,1),markerSize','MarkerFaceColor',cmap(1,:),'MarkerEdgeColor','k');
    
    errorbar(ap_axis,fluo_mean_array(plot_time_indices(2),:,1),fluo_ste_array(plot_time_indices(2),:,1),'-','Color','k','Capsize',0)
    s(2) = scatter(ap_axis,fluo_mean_array(plot_time_indices(2),:,1),markerSize','MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k');
    
    legend(s,'WT (nc13)','WT (nc14)','Location','northeast')
    
    xlabel('AP position (% embryo length)')
    ylabel('average spot fluorescence (au)')
    set(gca,'Fontsize',14)
    set(gca,'Color',[228,221,209]/255) 
    grid on
    fluo_fig.InvertHardcopy = 'off';
    set(gcf,'color','w');
    xlim([20 50])
    saveas(fluo_fig,[FigurePath 'fluo_vs_ap_time' num2str(time_axis(a)) '.png'])
    
end    