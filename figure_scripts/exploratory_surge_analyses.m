% Script to generate figures establishing presence of enrichment bursts at
% start of transcription bursts
clear 
close all
addpath(genpath('../utilities'));

% set ID variables
targetProjectName = '2xDl-Ven_snaBAC-mCh';
controlProjectName = '2xDl-Ven_hbP2P-mCh';

projectName_cell = {targetProjectName controlProjectName};% targetProjectOrig};


% load data
master_struct = struct;
for i = 1:length(projectName_cell)
  
  % set write paths
  liveProject = LiveEnrichmentProject(projectName_cell{i});
  resultsRoot = [liveProject.dataPath filesep];
  resultsDir = [resultsRoot 'cpHMM_results' filesep];
  
  % load data
  load([resultsDir 'hmm_input_output_results.mat'])
  fieldNames = fieldnames(results_struct);
  for f = 1:length(fieldNames)
      master_struct(i).(fieldNames{f}) = results_struct.(fieldNames{f});
  end
  clear results_struct;
  
  % make figure directory
  if strcmpi(projectName_cell{i},targetProjectName)
      FigurePath = [liveProject.figurePath 'input_output' filesep];
      mkdir(FigurePath)
  end
end


Tres = 20; % seconds

%  determine snip size
n_col = size(master_struct(1).hmm_array,2);
window_size = floor(n_col/2);
time_axis = (-window_size:window_size)*Tres/60;

% set basic analyisis parameters
nBoots = 100; % number of bootstrap samples to use
min_pause_len = 5; % minimum length of preceding OFF period (in time steps)
max_pause_len = 100;
min_burst_len = 3;
max_burst_len = 1000;


for i = 1:length(master_struct)
  % generate basic filter for target locus and computational controls
  master_struct(i).burst_ft = master_struct(i).feature_sign_vec == 1&master_struct(i).lead_dur_vec>=...
    min_pause_len&master_struct(i).lead_dur_vec<=max_pause_len...
      & master_struct(i).lag_dur_vec>=min_burst_len&master_struct(i).lag_dur_vec<=max_burst_len;%&target_swap_qc&target_virtual_qc;; % filter for rise events
    
  sample_options = find(master_struct(i).burst_ft);

  % target
  master_struct(i).burst_rise_spot_array = NaN(nBoots,n_col);
  master_struct(i).burst_rise_delta_array = NaN(nBoots,n_col);
  master_struct(i).burst_rise_delta_vec = NaN(nBoots,n_col);
  master_struct(i).burst_rise_virt_array = NaN(nBoots,n_col);
  master_struct(i).burst_rise_hmm_array = NaN(nBoots,n_col);
  
  % random
  master_struct(i).burst_rise_spot_rand_array = NaN(nBoots,n_col);
  master_struct(i).burst_rise_delta_rand_array = NaN(nBoots,n_col);
  master_struct(i).burst_rise_delta_rand_vec = NaN(nBoots,1);
  master_struct(i).burst_rise_virt_rand_array = NaN(nBoots,n_col);
  master_struct(i).burst_rise_hmm_rand_array = NaN(nBoots,n_col);
  
  % take bootstrap samples
  for n = 1:nBoots
    
      % draw random subsample for qualifying events
      s_ids_target = randsample(sample_options,numel(sample_options),true);    
      master_struct(i).burst_rise_hmm_array(n,:) = nanmean(master_struct(i).hmm_array(s_ids_target,:));
      master_struct(i).burst_rise_spot_array(n,:) = nanmean(master_struct(i).spot_array_dt(s_ids_target,:));
      master_struct(i).burst_rise_nucleus(n,:) = nanmean(nanmean(master_struct(i).nucleus_array(s_ids_target,:)));
      master_struct(i).burst_rise_virt_array(n,:) = nanmean(master_struct(i).virtual_array_dt(s_ids_target,:));    
      
      delta_array = master_struct(i).spot_array_dm(s_ids_target,:)-master_struct(i).virtual_array_dm(s_ids_target,:);
      master_struct(i).burst_rise_delta_array(n,:) = nanmean(delta_array);    
      master_struct(i).burst_rise_delta_vec(n) = nanmean(nanmean(delta_array(:,window_size+1:window_size+3),2)-nanmean(delta_array(:,window_size-2:window_size),2));    
      
      % as a control draw random subsample from all events
      s_ids_random = randsample(1:length(master_struct(i).lead_dur_vec),length(sample_options),true);    
      master_struct(i).burst_rise_hmm_rand_array(n,:) = nanmean(master_struct(i).hmm_array(s_ids_random,:));
      master_struct(i).burst_rise_spot_rand_array(n,:) = nanmean(master_struct(i).spot_array_dt(s_ids_random,:));      
      master_struct(i).burst_rise_virt_rand_array(n,:) = nanmean(master_struct(i).virtual_array_dt(s_ids_random,:));     
      
      delta_array_rand = master_struct(i).spot_array_dm(s_ids_random,:)-master_struct(i).virtual_array_dm(s_ids_random,:);
      master_struct(i).burst_rise_delta_rand_array(n,:) = nanmean(delta_array_rand);
      master_struct(i).burst_rise_delta_rand_vec(n) = nanmean(nanmean(delta_array_rand(:,window_size+1:window_size+3),2)-nanmean(delta_array_rand(:,window_size-2:window_size),2));    
  end

  % HMM trends
  
  % actual events
  master_struct(i).burst_rise_hmm_mean = nanmean(master_struct(i).burst_rise_hmm_array);
  master_struct(i).burst_rise_hmm_ste = nanstd(master_struct(i).burst_rise_hmm_array);
  
  % randomized events
  master_struct(i).burst_rise_hmm_rand_mean = nanmean(master_struct(i).burst_rise_hmm_rand_array);
  master_struct(i).burst_rise_hmm_rand_ste = nanstd(master_struct(i).burst_rise_hmm_rand_array);
  
  %%% calculate mean and standard error for spot
  
  % actual events
  master_struct(i).burst_rise_spot_mean = nanmean(master_struct(i).burst_rise_spot_array);
  master_struct(i).burst_rise_spot_ste = nanstd(master_struct(i).burst_rise_spot_array);
  master_struct(i).br_spot_ub = master_struct(i).burst_rise_spot_mean + master_struct(i).burst_rise_spot_ste;
  master_struct(i).br_spot_lb = master_struct(i).burst_rise_spot_mean - master_struct(i).burst_rise_spot_ste;

  master_struct(i).burst_rise_delta_mean = nanmean(master_struct(i).burst_rise_delta_array);
  master_struct(i).burst_rise_delta_ste = nanstd(master_struct(i).burst_rise_delta_array);
  master_struct(i).br_delta_ub = master_struct(i).burst_rise_delta_mean + master_struct(i).burst_rise_delta_ste;
  master_struct(i).br_delta_lb = master_struct(i).burst_rise_delta_mean - master_struct(i).burst_rise_delta_ste;
  
  % random events
  master_struct(i).burst_rise_spot_rand_mean = nanmean(master_struct(i).burst_rise_spot_rand_array);
  master_struct(i).burst_rise_spot_rand_ste = nanstd(master_struct(i).burst_rise_spot_rand_array);
  master_struct(i).br_spot_rand_ub = master_struct(i).burst_rise_spot_rand_mean + master_struct(i).burst_rise_spot_rand_ste;
  master_struct(i).br_spot_rand_lb = master_struct(i).burst_rise_spot_rand_mean - master_struct(i).burst_rise_spot_rand_ste;
  
  master_struct(i).burst_rise_delta_rand_mean = nanmean(master_struct(i).burst_rise_delta_rand_array);
  master_struct(i).burst_rise_delta_rand_ste = nanstd(master_struct(i).burst_rise_delta_rand_array);
  
  %%% calculate mean and standard error for virtual spot
  
  % actual events
  master_struct(i).burst_rise_virt_mean = nanmean(master_struct(i).burst_rise_virt_array);
  master_struct(i).burst_rise_virt_ste = nanstd(master_struct(i).burst_rise_virt_array);
  master_struct(i).br_virt_ub = master_struct(i).burst_rise_virt_mean + master_struct(i).burst_rise_virt_ste;
  master_struct(i).br_virt_lb = master_struct(i).burst_rise_virt_mean - master_struct(i).burst_rise_virt_ste;
  
  % random events
  master_struct(i).burst_rise_virt_rand_mean = nanmean(master_struct(i).burst_rise_virt_rand_array);
  master_struct(i).burst_rise_virt_rand_ste = nanstd(master_struct(i).burst_rise_virt_rand_array);
  master_struct(i).br_virt_rand_ub = master_struct(i).burst_rise_virt_rand_mean + master_struct(i).burst_rise_virt_rand_ste;
  master_struct(i).br_virt_rand_lb = master_struct(i).burst_rise_virt_rand_mean - master_struct(i).burst_rise_virt_rand_ste;
end



%% make figure
close all

basic_surge_fig = figure;
cmap1 = brewermap([],'Set2');

% snail activity
yyaxis right
plot(time_axis,master_struct(1).burst_rise_hmm_mean,'--','LineWidth',2,'Color','black');
ylabel('snail transcription (au)')
ylim([-2 20])
% set(gca,'ytick',.1:.2:1.1)
ax = gca;
ax.YColor = 'black';

% Dorsal activity
yyaxis left
hold on

% bio control
fill([time_axis fliplr(time_axis)],[master_struct(2).br_spot_ub fliplr(master_struct(2).br_spot_lb)],cmap1(5,:),'FaceAlpha',.5,'EdgeAlpha',0)
p3 = plot(time_axis,master_struct(2).burst_rise_spot_mean,'-','Color',cmap1(5,:),'LineWidth',2);
ylabel('de-trended Dl concentration (au)')

% locus
fill([time_axis fliplr(time_axis)],[master_struct(1).br_spot_ub fliplr(master_struct(1).br_spot_lb)],cmap1(2,:),'FaceAlpha',.5,'EdgeAlpha',0)
p2 = plot(time_axis,master_struct(1).burst_rise_spot_mean,'-','Color',cmap1(2,:),'LineWidth',2);

% locus (randomized)
fill([time_axis fliplr(time_axis)],[master_struct(1).br_virt_ub fliplr(master_struct(1).br_virt_lb)],cmap1(3,:),'FaceAlpha',.5,'EdgeAlpha',0)
p1 = plot(time_axis,master_struct(1).burst_rise_virt_mean,'-','Color',cmap1(3,:),'LineWidth',2);

% ylim([-20 25])
% set(gca,'ytick',-20:5:25)
ax = gca;
ax.YColor = 'black';%cmap1(2,:);
% grid on
xlabel('offset (minutes)')
legend([p2 p3 p1],'Dl at {\it snail} locus', 'Dl at {\it hbP2P}','virtual spot', 'Location','southwest');
set(gca,'Fontsize',14,'xtick',-4:2:4)
chH = get(gca,'Children');
set(gca,'Children',flipud(chH));
% ylim([-.05 .07])
set(gca,    'Box','off',...
            'Color',[228,221,209]/255,...            
            'TickLength',[0.02,0.05])    
basic_surge_fig.Color = 'white';        
basic_surge_fig.InvertHardcopy = 'off';
% save
saveas(basic_surge_fig,[FigurePath 'snail_w_hbP2P_control.tif'])
saveas(basic_surge_fig,[FigurePath 'snail_w_hbP2P_control.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compare to original results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:2
    delta_surge_fig = figure;
    cmap1 = brewermap([],'Set2');

    % snail activity
    yyaxis right
    plot(time_axis,master_struct(1).burst_rise_hmm_mean,'--','LineWidth',2,'Color','black');
    ylabel('snail transcription (au)')
    ylim([-2 20])
    % set(gca,'ytick',.1:.2:1.1)
    ax = gca;
    ax.YColor = 'black';

    % Dorsal activity
    yyaxis left
    hold on

    % bio control
    if i == 2
        fill([time_axis fliplr(time_axis)],[master_struct(2).br_delta_ub fliplr(master_struct(2).br_delta_lb)],cmap1(5,:),'FaceAlpha',.5,'EdgeAlpha',0)
        p3 = plot(time_axis,master_struct(2).burst_rise_delta_mean,'-','Color',cmap1(5,:),'LineWidth',2);
    end
    ylabel('relative Dl concentration (au)')

    % locus
    fill([time_axis fliplr(time_axis)],[master_struct(1).br_delta_ub fliplr(master_struct(1).br_delta_lb)],cmap1(2,:),'FaceAlpha',.5,'EdgeAlpha',0)
    p2 = plot(time_axis,master_struct(1).burst_rise_delta_mean,'-','Color',cmap1(2,:),'LineWidth',2);

    % % locus (randomized)
    % p1 = plot(time_axis,master_struct(1).burst_rise_spot_rand_mean,'-','Color',cmap1(1,:),'LineWidth',2);

    % ylim([-20 25])
    % set(gca,'ytick',-20:5:25)
    ax = gca;
    ax.YColor = 'black';%cmap1(2,:);
    % grid on
    xlabel('offset (minutes)')
    if i == 2
        legend([p2 p3],'Dl at {\it snail} locus', 'Dl at {\it hbP2P}', 'Location','southwest');
    end
    set(gca,'Fontsize',14,'xtick',-4:2:4)
    chH = get(gca,'Children');
    set(gca,'Children',flipud(chH));
    % ylim([-.05 .07])
    set(gca,    'Box','off',...
                'Color',[228,221,209]/255,...            
                'TickLength',[0.02,0.05])    
    delta_surge_fig.Color = 'white';        
    delta_surge_fig.InvertHardcopy = 'off';
    if i == 2
        % save
        saveas(delta_surge_fig,[FigurePath 'dl_delta_snail_w_hbP2P_control.tif'])
        saveas(delta_surge_fig,[FigurePath 'dl_delta_snail_w_hbP2P_control.pdf'])
    else
        saveas(delta_surge_fig,[FigurePath 'dl_delta_snail.tif'])
        saveas(delta_surge_fig,[FigurePath 'dl_delta_snail.pdf'])
    end
end        