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
projectNameCell = {'Bcd-GFP_hbP2P38F1-mCh_Leica','Bcd-GFP_hbP2Pvk33-mCh_Leica','Bcd-GFP_hbP2Pvk33-mCh_AiryHP',...
                   'Bcd-GFP_hbP2P-mCh','Bcd-GFP_snailMS2-mCh_AiryHP','Bcd-GFP_hbMS2-mCh_Airy_fast',...
                   'Bcd-GFP_hbMS2-mCh_NoAiry_02','Bcd-GFP_hbMS2-mCh_AiryscanTest_','Bcd-GFP_snaBAC-mCh','Bcd-GFP-McpMcherry-hbP2P-delta6'};%'Bcd-GFP_hbP2P-mCh';%
                 
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

% specify basic qc constraints
DistLim = 0.8; % minimum distance to nuclear boundary we'll permit
refAPWindow = [0.2 ,0.3]*100; % control for AP position
refTimeWindow = [5, 20]*60; % control for time period
nBoots = 100;
minDP = 50;

% Loop through each project and extract key data vectors
for  p = 1:length(master_struct)
    spot_struct_protein = master_struct(p).spot_struct_protein;
    traceQCVec = [spot_struct_protein.TraceQCFlag];
    spot_struct = master_struct(p).spot_struct;
    particleIDVec = [spot_struct.particleID];
    
    % basic protein vectors
    spot_protein_vec = [spot_struct_protein(traceQCVec).spot_protein_vec];
    ctrl_protein_vec = [spot_struct_protein(traceQCVec).edge_null_protein_vec];
    nucleus_protein_vec = [spot_struct_protein(traceQCVec).nuclear_protein_vec];
  
    % spot intensity
    fluo_vec = [spot_struct_protein(traceQCVec).fluo];
    
    % generate id vector to allow us to map back to individual spots
    spot_id_vec = [];
%     nc_lin_vec = [];
    nc_sub_lin_vec = [];
    nuclear_cycle_vec = [];
    for i = find(traceQCVec)
        spot_id_vec = [spot_id_vec repelem(i,length(spot_struct_protein(i).nuclear_protein_vec))];
        nc_sub_lin_vec = [nc_sub_lin_vec 1:length(spot_struct_protein(i).nuclear_protein_vec)];
        % add nuclear cycle info
        alt_ind = particleIDVec==spot_struct_protein(i).particleID;
        nuclear_cycle_vec = [nuclear_cycle_vec repelem(spot_struct(alt_ind).nc,length(spot_struct_protein(i).nuclear_protein_vec))];
    end
    % time and position vectors (if present)
    ap_flag = true;
    if isfield(spot_struct_protein,'APPosNucleus') && ~ismember(p,[4,5,9,3]) %NL: these sets do not have reliable AP info
        ap_pos_vec = [spot_struct_protein(traceQCVec).APPosParticle];
    else
        ap_flag = false;
        ap_pos_vec = repelem(mean(refAPWindow),length(fluo_vec));
    end
    time_vec = [spot_struct_protein(traceQCVec).time];
    
    % check distance from nucleus boundary
    spot_edge_dist_vec = [spot_struct_protein(traceQCVec).spot_edge_dist_vec]*master_struct(p).PixelSize;
    spot_edge_dist_flags = spot_edge_dist_vec>=DistLim;
    
    % generate AP and time-related flags
    ap_keep_flags = ap_pos_vec>=refAPWindow(1)&ap_pos_vec<=refAPWindow(2);
    time_keep_flags = time_vec>=refTimeWindow(1)&ap_pos_vec<=refTimeWindow(2);
    
    % store filtered longform vecors
    master_struct(p).spot_protein_vec = spot_protein_vec(spot_edge_dist_flags&ap_keep_flags&time_keep_flags);
    master_struct(p).ctrl_protein_vec = ctrl_protein_vec(spot_edge_dist_flags&ap_keep_flags&time_keep_flags);
    master_struct(p).nucleus_protein_vec = nucleus_protein_vec(spot_edge_dist_flags&ap_keep_flags&time_keep_flags);
    
    bcd_denominator = nanmean(master_struct(p).nucleus_protein_vec);
    
    master_struct(p).fluo_vec = fluo_vec(spot_edge_dist_flags&ap_keep_flags&time_keep_flags);
    master_struct(p).time_vec = time_vec(spot_edge_dist_flags&ap_keep_flags&time_keep_flags);
    master_struct(p).ap_pos_vec = ap_pos_vec(spot_edge_dist_flags&ap_keep_flags&time_keep_flags);
    master_struct(p).spot_id_vec = spot_id_vec(spot_edge_dist_flags&ap_keep_flags&time_keep_flags);  
    master_struct(p).nc_sub_lin_vec = nc_sub_lin_vec(spot_edge_dist_flags&ap_keep_flags&time_keep_flags);  
    master_struct(p).nuclear_cycle_vec = nuclear_cycle_vec(spot_edge_dist_flags&ap_keep_flags&time_keep_flags);  
    master_keep_indices = find(spot_edge_dist_flags&ap_keep_flags&time_keep_flags);
    
    master_struct(p).ap_flag = ap_flag;
    
    % store locus-level enrichment stats
    
    master_struct(p).delta_protein_vec = (master_struct(p).spot_protein_vec-master_struct(p).ctrl_protein_vec)/bcd_denominator;
    spot_index = unique(master_struct(p).spot_id_vec);
    master_struct(p).mean_delta_vec = NaN(size(spot_index));
    master_struct(p).mean_spot_protein_vec = NaN(size(spot_index));
    master_struct(p).mean_ctrl_protein_vec = NaN(size(spot_index));
    
    for s = 1:length(spot_index)
        master_struct(p).mean_delta_vec(s) = nanmean(master_struct(p).delta_protein_vec(master_struct(p).spot_id_vec==spot_index(s)));
        master_struct(p).mean_spot_protein_vec(s) = nanmean(master_struct(p).spot_protein_vec(master_struct(p).spot_id_vec==spot_index(s)));
        master_struct(p).mean_ctrl_protein_vec(s) = nanmean(master_struct(p).ctrl_protein_vec(master_struct(p).spot_id_vec==spot_index(s)));
    end        
    master_struct(p).bcd_denominator = bcd_denominator;
    
    % calculate the mean and standard error across spots
    mean_delta_array = NaN(1,nBoots);
    boot_options = 1:length(master_struct(p).mean_delta_vec);
    for n = 1:nBoots
        boot_indices = randsample(boot_options,length(boot_options),true);
        mean_delta_array(n) = nanmean(master_struct(p).mean_delta_vec(boot_indices));
    end
    master_struct(p).delta_mean = nanmean(mean_delta_array);
    master_struct(p).delta_ste = nanstd(mean_delta_array);
    
    % Generate filtered snips
    snip_data = master_struct(p).snip_data;
    
    spot_snips = cat(3,snip_data.spot_mcp_snips);
    ctrl_snips = cat(3,snip_data.edge_control_mcp_snips);
    
    nc_lin_vec_snip = [snip_data.nc_lin_index_vec];
    nc_sub_lin_vec_snip = [snip_data.nc_sub_index_vec];
    nc_master_vec = [snip_data.nc_master_vec];
    
    snip_filter = ismember([nc_lin_vec_snip' nc_sub_lin_vec_snip'],[master_struct(p).spot_id_vec' master_struct(p).nc_sub_lin_vec'],'rows');
    
    master_struct(p).spot_snips = spot_snips(:,:,snip_filter);
    master_struct(p).ctrl_snips = ctrl_snips(:,:,snip_filter);
    
end

%%  Make some box plots
markerSize = 100;
close all
master_indices = [4 9 8 6 3 5 1 2 10];
% master_indices_3 = 3;
color_index_3 = [2 3 8 4 2 3 4 2 5];
shape_index_3 = [1 1 3 3 3 3 2 2 3];
shape_cell = {'o','s','d'};
x_string_cell = {'@hb-vk33-Leica (2019)','@snail-Leica (2019)','@hb-38F1-980-SlowAiry','@hb-38F1-980-FastAiry',...
                   '@hb-vk33-980-FastAiry','@snail-980-FastAiry','@hb-38F1-Leica (2021)','@hb-vk33-Leica (2021)',...
                   '@hb-38F1-D6-980-FastAiry'};
                 
instruction_cell = {1:2,1:3,1:4,1:5,1:6,1:8,[4,6,9]};
for a = 1:length(instruction_cell)
    
    local_indices = instruction_cell{a};
    iter_indices = master_indices(local_indices);
    x_string_iter = x_string_cell(local_indices);
    
    % generate normalized arrays                
    
    scatter_fig = figure;
    hold on

    cmap = brewermap([],'Set2');
    errorbar(1:length(iter_indices),[master_struct(iter_indices).delta_mean],[master_struct(iter_indices).delta_ste],'o','Color','k','LineWidth',1.5)
    for i = 1:length(iter_indices)
        scatter(i,master_struct(iter_indices(i)).delta_mean,markerSize,shape_cell{shape_index_3(local_indices(i))},...
                  'MarkerFaceColor',cmap(color_index_3(local_indices(i)),:),'MarkerEdgeColor','k')
    end
    
    plot(0:25,zeros(1,26),'--k')
    
    ylabel({'Bcd-GFP enrichment' ; '(Bcd_{spot}-Bcd_{ctrl})/Bcd_{nuc}'})
    
    if a < length(instruction_cell)
        xlim([0.5 length(master_indices)+.5])        
    else
        xlim([0.5 length(x_string_iter)+.5])
    end
    ylim([-0.03 0.09])
    set(gca,'xtick',1:length(x_string_iter),'xticklabels',x_string_iter)
    set(gca,'FontSize',12)
    set(gca,'Color',[228,221,209]/255) 
    grid on
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.XAxis(1).Color = 'k';
    scatter_fig.InvertHardcopy = 'off';
    set(gcf,'color','w');

    
    xtickangle(-30)
    saveas(scatter_fig,[FigurePath 'scatter_' num2str(iter_indices) '.png'])
end

