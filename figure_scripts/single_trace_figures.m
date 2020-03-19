% Script to generate figures comparing MS2 fitting methodologies
clear 
close all
addpath('../utilities')
% set ID variables
DropboxFolder = 'S:\Nick\Dropbox\';
project_cell = {'Dl-Ven_snaBAC-mCh_v4','Dl-Ven_snaBAC-mCh_v4'};%,'Dl-Ven_snaBAC-mCh_F-F-F_v1'};
title_cell = {'OG (old)','OG (new)'};
fluo_dim_vec = [2,3];
protein_dim_vec = [3,3];
type_name = 'OG_only';
% Params
K = 3;
w = 7;

% load data for each project
master_struct = struct;
for p = 1:numel(project_cell)
    project = project_cell{p};
    fluo_dim = fluo_dim_vec(p);
    protein_dim = protein_dim_vec(p);
    % set write paths
    [~, DataPath, FigureRoot] =   header_function(DropboxFolder, project); 
    
    % intermediate HMM results
    load([DataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '_f' num2str(fluo_dim) 'D_p' num2str(protein_dim) 'D_dt.mat'])
    master_struct(p).hmm_struct = hmm_input_output;
    clear hmm_input_output;
    
    % final results
    load([DataPath 'hmm_input_output_results_w' num2str(w) '_K' num2str(K) '_f' num2str(fluo_dim) 'D_p' num2str(protein_dim) 'D.mat'])
    master_struct(p).results_struct = results_struct;
    clear results_struct;

end

%% make figure directory
FigPath = [FigureRoot '\spot_fit_comparisons\' type_name '\'];
mkdir(FigPath)


% create analysis structure containing list of burst features
analysis_struct = struct;
Tres = 20; % seconds
window_size = 15;
time_axis = (-window_size:window_size)*Tres / 60;

min_pause_len = 5; % minimum length of preceding OFF period (in time steps)
max_pause_len = 1000;
min_burst_len = 2;
max_burst_len = 1000;

for p = 1:numel(project_cell)
    
    % extract relevant arrays from project 
    results_struct = master_struct(p).results_struct;
    analysis_struct(p).lag_dur_vec = results_struct.lag_dur_vec;
    analysis_struct(p).lead_dur_vec = results_struct.lead_dur_vec;
    analysis_struct(p).hmm_array_dm = results_struct.hmm_array;
    analysis_struct(p).hmm_array_dm = analysis_struct(p).hmm_array_dm ./ nanstd(analysis_struct(p).hmm_array_dm);
    fluo_array_dm = results_struct.fluo_array  - nanmean(results_struct.fluo_array,2);
    analysis_struct(p).fluo_array_dm = fluo_array_dm / nanstd(fluo_array_dm(:));
    analysis_struct(p).mf_array_dm = results_struct.mf_array - nanmean(results_struct.mf_array,2);     
    analysis_struct(p).time_vec = results_struct.center_time_vec/60;
    analysis_struct(p).spot_array_dm = results_struct.spot_array_dm;
    analysis_struct(p).virtual_array_dm = results_struct.virtual_array_dm;
    analysis_struct(p).feature_sign_vec = results_struct.feature_sign_vec;
    analysis_struct(p).particle_id_vec = results_struct.particle_id_vec;
    analysis_struct(p).particle_sub_id_vec = results_struct.particle_sub_id_vec;
    
    % generate basic filter for target locus and computational controls
    analysis_struct(p).burst_ft = results_struct.feature_sign_vec == 1&results_struct.lead_dur_vec>=min_pause_len&results_struct.lead_dur_vec<=max_pause_len...
        &results_struct.lag_dur_vec>=min_burst_len&results_struct.lag_dur_vec<=max_burst_len;%  
    % less restirctive filter taking all bursts
    analysis_struct(p).burst_ft_all = results_struct.feature_sign_vec == 1&results_struct.lead_dur_vec>=1;
    
    % get list of particle ids amd event times
    % record sampling vector
    analysis_struct(p).sample_options = find(analysis_struct(p).burst_ft);
    % particle ids
    analysis_struct(p).surge_particles = analysis_struct(p).particle_id_vec(analysis_struct(p).burst_ft);
    analysis_struct(p).surge_sub_particles = analysis_struct(p).particle_sub_id_vec(analysis_struct(p).burst_ft);
%     % event times
%     analysis_struct(p).surge_burst_times = analysis_struct(p).center_time(analysis_struct(p).burst_ft);    
end

% get unique list of particles with at least 1 burst used in surge analysis
[surge_particle_array, indices] = unique([[analysis_struct.surge_particles]'  [analysis_struct.surge_sub_particles]'],'rows');
%% iterate through traces and make plots
hmm_particle_vec = [master_struct(1).hmm_struct.ParticleID];
hmm_particle_sub_vec = NaN(size(hmm_particle_vec));
for i = 1:numel(hmm_particle_vec)
    hmm_particle_sub_vec(i) = sum(hmm_particle_vec(1:i)==hmm_particle_vec(i));
end

%  calculate mapping between fluo vectors 
fluo_all1 = [master_struct(1).hmm_struct.fluo];
fluo_all2 = [master_struct(2).hmm_struct.fluo];

fluo_factor = fluo_all1' \ fluo_all2';


rng(341);
n_samples = 250;
sample_vec = 1:size(surge_particle_array,1);
plot_indices = randsample(sample_vec,n_samples,false);
cmap = brewermap(9,'Set2');
for p = plot_indices
    % find ids that correspond to this particle
    pt_id = surge_particle_array(p,1);
    pt_sub_id = surge_particle_array(p,2);
    
    hmm_filter = hmm_particle_vec==pt_id & hmm_particle_sub_vec==pt_sub_id;
    
    % extract basic time, fluo, and protein fields
    time = master_struct(1).hmm_struct(hmm_filter).time_raw/60;
    time_interp = master_struct(1).hmm_struct(hmm_filter).time/60;
    
    fluo1 = master_struct(1).hmm_struct(hmm_filter).fluo_raw*fluo_factor;
    fluo1_interp = master_struct(1).hmm_struct(hmm_filter).fluo*fluo_factor;
    
    fluo2 = master_struct(2).hmm_struct(hmm_filter).fluo_raw;
    fluo2_interp = master_struct(2).hmm_struct(hmm_filter).fluo;
    
    protein = master_struct(1).hmm_struct(hmm_filter).spot_protein_dt; % note: should be same for both projects
    protein = protein - nanmin(protein);
    
    % get HMM info (project 1)
    all_burst_ft1 = master_struct(1).hmm_struct(hmm_filter).z_diff_vec>0;
    all_burst_times1 = analysis_struct(1).time_vec(all_burst_ft1);
    
    surge_burst_ft1 = analysis_struct(1).burst_ft & analysis_struct(1).particle_id_vec == pt_id & analysis_struct(1).particle_sub_id_vec == pt_sub_id;
    surge_burst_times1 = analysis_struct(1).time_vec(surge_burst_ft1);
    
     % get HMM info (project 2)
    all_burst_ft2 = master_struct(2).hmm_struct(hmm_filter).z_diff_vec>0;
    all_burst_times2 = analysis_struct(2).time_vec(all_burst_ft2);
    if max(find(all_burst_ft2)) > max(find(~isnan(fluo2_interp)))
        error('wtf')
    end
    surge_burst_ft2 = analysis_struct(2).burst_ft & analysis_struct(2).particle_id_vec == pt_id & analysis_struct(2).particle_sub_id_vec == pt_sub_id;
    surge_burst_times2 = analysis_struct(2).time_vec(surge_burst_ft2);
       
    % make figure
    burst_surge_fig = figure('Visible','off');
    hold on   
    
    yyaxis left
%     p_pt = plot(time_interp,protein,'--','Color',cmap(2,:),'LineWidth',1);
    p_pt = area(time_interp,protein,'FaceColor',cmap(2,:),'EdgeAlpha',0,'FaceAlpha',0.25);
    ax = gca;
    ax.YColor = cmap(2,:);
    ylabel('Dorsal al locus (detrended)')
    
    yyaxis right
    % fluo trends
    p2D = plot(time_interp,fluo1_interp,'-','Color','black','LineWidth',1.5);
%     scatter(time,fluo1,20,'MarkerFaceColor','black','MarkerEdgeAlpha',0)
        
    p3D = plot(time_interp,fluo2_interp,'-','Color',cmap(3,:),'LineWidth',1.5);
%     scatter(time,fluo2,20,'MarkerFaceColor',cmap(8,:),'MarkerEdgeAlpha',0)
    
    % plot inferred burst events
    ax = gca;
    y_lim = ax.YLim;
    ylim([0 max(y_lim)])
    y_lim = ax.YLim;   
    ax.YColor = 'black';
    xlim([min(time_interp) max(time_interp)])
%     
%     plot([all_burst_times2;all_burst_times2],[repelem(y_lim(1),numel(all_burst_times2)) ; ...
%         repelem(y_lim(2),numel(all_burst_times2))],'-','Color',[cmap(3,:) 1])
%     plot([all_burst_times1;all_burst_times1],[repelem(y_lim(1),numel(all_burst_times1)) ; ...
%         repelem(y_lim(2),numel(all_burst_times1))],'--','Color',[0 0 0 1])
    
    plot([surge_burst_times2;surge_burst_times2],[repelem(y_lim(1),numel(surge_burst_times2)) ; ...
        repelem(y_lim(2),numel(surge_burst_times2))],'-','Color',[cmap(3,:) 1],'LineWidth',1.5)
    plot([surge_burst_times1;surge_burst_times1],[repelem(y_lim(1),numel(surge_burst_times1)) ; ...
        repelem(y_lim(2),numel(surge_burst_times1))],'--','Color',[0 0 0 1],'LineWidth',1.5)
    
    ylabel('MS2 spot intensity')
    chi=get(gca, 'Children');
    set(gca, 'Children',flipud(chi));
    legend([p2D p3D p_pt],'2D method','3D method','Dorsal')    
    xlabel('minutes into nc14')
    set(gca,'Fontsize',14);
    box on;
    
    saveas(burst_surge_fig,[FigPath 'trace_pt' num2str(pt_id*1e5) '_sub' num2str(pt_sub_id) '.png'])
    close all
    
    
%     burst_surge_fig = figure('Visible','off');
%     hold on   
%     
%     yyaxis right
%     p_pt = plot(time_interp,protein,'--','Color',cmap(2,:),'LineWidth',1);
%     ax = gca;
%     ax.YColor = cmap(2,:);
%     ylabel('Dorsal al locus (detrended)')
%     
%     yyaxis left
%     % fluo trends
%     p2D = plot(time_interp,fluo1_interp,'-','Color','black','LineWidth',1.5);
% %     scatter(time,fluo1,20,'MarkerFaceColor','black','MarkerEdgeAlpha',0)
%         
%     p3D = plot(time_interp,fluo2_interp,'-','Color',cmap(3,:),'LineWidth',1.5);
% %     scatter(time,fluo2,20,'MarkerFaceColor',cmap(8,:),'MarkerEdgeAlpha',0)
%     
%     % plot inferred burst events
%     ax = gca;
%     y_lim = ax.YLim;
%     ylim([0 max(y_lim)])
%     y_lim = ax.YLim;   
%     ax.YColor = 'black';   
%     
%     ylabel('MS2 spot intensity')
%            
%     legend([p2D p3D p_pt],'2D method','3D method','Dorsal')    
%     xlabel('minutes into nc14')
%     set(gca,'Fontsize',14);
%     box on;
%     
%     saveas(burst_surge_fig,[FigPath 'simp_trace_pt' num2str(pt_id*1e5) '_sub' num2str(pt_sub_id) '.png'])
    close all
end