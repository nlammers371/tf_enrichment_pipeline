% Script to conduct exploratory analyses regarding correlation between
% Dorsal concentration and snail transcription
clear 
close all
addpath(genpath('../utilities'));

% set ID variables
targetProjectName = 'Bcd-GFP_hbP2P-mCh';%'2xDl-Ven_snaBAC-mCh';
controlProjectName = 'Bcd-GFP_snaBAC-mCh';%'2xDl-Ven_hbP2P-mCh';

projectName_cell = {targetProjectName controlProjectName};% targetProjectOrig};


% load data
master_input_output = cell(1,2);
master_protein = cell(1,2);
for i = 1:length(projectName_cell)
  
  % set write paths
  liveProject = LiveEnrichmentProject(projectName_cell{i});
  resultsRoot = [liveProject.dataPath filesep];
  resultsDir = [resultsRoot 'cpHMM_results' filesep];
  
  % load data
%   load([resultsDir 'hmm_input_output.mat'])  
%   master_input_output{i} = hmm_input_output;
%   clear hmm_input_output;
  
  load([resultsRoot 'spot_struct_protein.mat'])  
  master_protein{i} = spot_struct_protein;
  clear spot_struct_protein;
  
  % make figure directory
  if strcmpi(projectName_cell{i},targetProjectName)
      FigurePath = [liveProject.figurePath 'input_output' filesep];
      mkdir(FigurePath)
  end
end

%% %%%%%%%%%%%%%%%%% Plot local enrichment vs nuclear Dl %%%%%%%%%%%%%%%%%%
% get pixel size
liveExperiment = liveProject.includedExperiments{1};
FrameInfo = getFrameInfo(liveExperiment);
PixelSize = FrameInfo(1).PixelSize;

% set parameters
n_bins = 25;
n_boots = 100;

% extract vectors
target_dist_flags = [master_protein{1}.spot_edge_dist_vec]*PixelSize<0.8 | [master_protein{1}.time]< 900 | [master_protein{1}.time]>1200;
target_ctrl_protein_vec = [master_protein{1}.edge_null_protein_vec];
target_spot_protein_vec = [master_protein{1}.spot_protein_vec];
target_nuclear_protein_vec = [master_protein{1}.nuclear_protein_vec];
target_fluo_vec = [master_protein{1}.fluo];
target_delta_vec = target_spot_protein_vec - target_ctrl_protein_vec;

control_dist_flags = [master_protein{2}.spot_edge_dist_vec]*PixelSize<0.8  | [master_protein{2}.time]< 900 & [master_protein{2}.time]>1200;
control_ctrl_protein_vec = [master_protein{2}.edge_null_protein_vec];
control_spot_protein_vec = [master_protein{2}.spot_protein_vec];
control_nuclear_protein_vec = [master_protein{2}.nuclear_protein_vec];
control_fluo_vec = [master_protein{2}.fluo];
control_delta_vec = control_spot_protein_vec - control_ctrl_protein_vec;

% create bins
protein_bins = linspace(prctile([target_nuclear_protein_vec(~target_dist_flags)],1),prctile([target_nuclear_protein_vec(~target_dist_flags)],99),n_bins + 1);
protein_bins_plot = protein_bins(1:end-1) + diff(protein_bins)/2;
% initialize bootstrap arrays
target_delta_array = NaN(n_boots,n_bins);
control_delta_array = NaN(n_boots,n_bins);

target_fluo_array = NaN(n_boots,n_bins);
control_fluo_array = NaN(n_boots,n_bins);

for b = 1:n_bins-1
    % find qualifying observations
    target_ids = find(target_nuclear_protein_vec<protein_bins(b+1)&target_nuclear_protein_vec>=protein_bins(b)&~target_dist_flags);
    control_ids = find(control_nuclear_protein_vec<protein_bins(b+1)&control_nuclear_protein_vec>=protein_bins(b)&~control_dist_flags);
    for n = 1:n_boots
        % draw sample
        target_boot_ids = randsample(target_ids,length(target_ids),true);
        control_boot_ids = randsample(control_ids,length(control_ids),true);
        % record enrichment info
        target_delta_array(n,b) = nanmean(target_delta_vec(target_boot_ids));
        control_delta_array(n,b) = nanmean(control_delta_vec(control_boot_ids));
        % record  fluo info
        target_fluo_array(n,b) = nanmean(target_fluo_vec(target_boot_ids));
        control_fluo_array(n,b) = nanmean(control_fluo_vec(control_boot_ids));
    end
end
        
 
%% make figures
close all

e_fig1 = figure;

cmap1 = brewermap([],'Set2');
hold on

plot([-.01+protein_bins_plot(1) .01+protein_bins_plot(end)],zeros(1,2),'--k')

e1 = errorbar(protein_bins_plot,nanmean(target_delta_array),nanstd(target_delta_array),'-','Color','k');
e1.CapSize = 0;
s1 = scatter(protein_bins_plot,nanmean(target_delta_array),'MarkerFaceColor',cmap1(2,:),'MarkerEdgeColor','k');

e2 = errorbar(protein_bins_plot,nanmean(control_delta_array),nanstd(control_delta_array),'-','Color','k');
e2.CapSize = 0;
s2 = scatter(protein_bins_plot,nanmean(control_delta_array),'MarkerFaceColor',cmap1(5,:),'MarkerEdgeColor','k');


ax = gca;
ax.YColor = 'black';%cmap1(2,:);
% grid on
xlabel('nuclear [Dl] (au)')
ylabel('Dl enrichment (au)')
legend([s1 s2],'{\it snail} locus', '{\it hbP2P} locus', 'Location','southwest');
set(gca,'Fontsize',14,'xtick',-4:2:4)

xlim([protein_bins_plot(1)-.01 .01+protein_bins_plot(end)])
set(gca,    'Box','off',...
            'Color',[228,221,209]/255,...            
            'TickLength',[0.02,0.05])    
e_fig1.Color = 'white';        
e_fig1.InvertHardcopy = 'off';
% save
saveas(e_fig1,[FigurePath 'enrichment_vs_nuclear_dl.tif'])
saveas(e_fig1,[FigurePath 'enrichment_vs_nuclear_dl.pdf'])
      


e_fig2 = figure;
hold on

% plot([-.01+protein_bins_plot(1) .01+protein_bins_plot(end)],zeros(1,2),'--k')

e1 = errorbar(protein_bins_plot,nanmean(target_fluo_array),nanstd(target_fluo_array),'-','Color','k');
e1.CapSize = 0;
s1 = scatter(protein_bins_plot,nanmean(target_fluo_array),'MarkerFaceColor',cmap1(2,:),'MarkerEdgeColor','k');

e2 = errorbar(protein_bins_plot,nanmean(control_fluo_array),nanstd(control_fluo_array),'-','Color','k');
e2.CapSize = 0;
s2 = scatter(protein_bins_plot,nanmean(control_fluo_array),'MarkerFaceColor',cmap1(5,:),'MarkerEdgeColor','k');


ax = gca;
ax.YColor = 'black';%cmap1(2,:);
% grid on
xlabel('nuclear [Dl] (au)')
ylabel('spot intensity (au)')
legend([s1 s2],'{\it snail} locus', '{\it hbP2P} locus', 'Location','northwest');
set(gca,'Fontsize',14,'xtick',-4:2:4)

xlim([protein_bins_plot(1)-.01 .01+protein_bins_plot(end)])
set(gca,    'Box','off',...
            'Color',[228,221,209]/255,...            
            'TickLength',[0.02,0.05])    
e_fig2.Color = 'white';        
e_fig2.InvertHardcopy = 'off';
% save
saveas(e_fig2,[FigurePath 'fluo_vs_nuclear_dl.tif'])
saveas(e_fig2,[FigurePath 'fluo_vs_nuclear_dl.pdf'])


close all

e_fig3 = figure;
hold on

% plot([-.01+protein_bins_plot(1) .01+protein_bins_plot(end)],zeros(1,2),'--k')
t_fluo_mean = nanmean(target_fluo_array);
t_fluo_ub =  nanstd(target_fluo_array);
t_fluo_lb = nanstd(target_fluo_array);
t_delta_mean = nanmean(target_delta_array);
t_delta_ub = nanstd(target_delta_array);
t_delta_lb = nanstd(target_delta_array);

c_fluo_mean = nanmean(control_fluo_array);
c_fluo_ub = nanstd(control_fluo_array);
c_fluo_lb = nanstd(control_fluo_array);
c_delta_mean = nanmean(control_delta_array);
c_delta_ub = nanstd(control_delta_array);
c_delta_lb = nanstd(control_delta_array);


e1 = errorbar(t_delta_mean,t_fluo_mean,t_fluo_lb,t_fluo_ub,t_delta_lb,t_delta_ub,'o','Color','k');
e1.CapSize = 0;
s1 = scatter(t_delta_mean,t_fluo_mean,'MarkerFaceColor',cmap1(2,:),'MarkerEdgeColor','k');

e2 = errorbar(c_delta_mean,c_fluo_mean,c_fluo_lb,c_fluo_ub,c_delta_lb,c_delta_ub,'o','Color','k');
e2.CapSize = 0;
s2 = scatter(c_delta_mean,c_fluo_mean,'MarkerFaceColor',cmap1(5,:),'MarkerEdgeColor','k');

ax = gca;
ax.YColor = 'black';%cmap1(2,:);
xlim([-.1 .35])
% grid on
xlabel('Dl enrichment (au)')
ylabel('spot intensity (au)')
legend([s1 s2],'{\it snail} locus', '{\it hbP2P} locus', 'Location','northwest');
set(gca,'Fontsize',14)

set(gca,    'Box','off',...
            'Color',[228,221,209]/255,...            
            'TickLength',[0.02,0.05])    
e_fig3.Color = 'white';        
e_fig3.InvertHardcopy = 'off';
% save
saveas(e_fig3,[FigurePath 'fluo_vs_enrichment.tif'])
saveas(e_fig3,[FigurePath 'fluo_vs_enrichment.pdf'])

%%
e_fig4 = figure;
hold on


e1 = errorbar(t_delta_mean,t_fluo_mean,t_fluo_lb,t_fluo_ub,t_delta_lb,t_delta_ub,'o','Color','k');
e1.CapSize = 0;
s1 = scatter(t_delta_mean,t_fluo_mean,'MarkerFaceColor',cmap1(2,:),'MarkerEdgeColor','k');

ax = gca;
ax.YColor = 'black';%cmap1(2,:);
xlim([.05 .35])
% grid on
xlabel('Dl enrichment at {\it snail} (au)')
ylabel('spot intensity (au)')
% legend([s1 s2],'{\it snail} locus', '{\it hbP2P} locus', 'Location','northwest');
set(gca,'Fontsize',14)

set(gca,    'Box','off',...
            'Color',[228,221,209]/255,...            
            'TickLength',[0.02,0.05])    
e_fig4.Color = 'white';        
e_fig4.InvertHardcopy = 'off';
% save
saveas(e_fig4,[FigurePath 'fluo_vs_enrichment_target_only.tif'])
saveas(e_fig4,[FigurePath 'fluo_vs_enrichment_target_only.pdf'])

%% %%%%%%%%%%%%%%%%%%%%% fit linear model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pd_vec_mean = target_ctrl_protein_vec(~target_dist_flags);
pd_vec_local = target_spot_protein_vec(~target_dist_flags);
y_vec = target_fluo_vec(~target_dist_flags);

mdl1 = fitlm(pd_vec_mean',y_vec')

mdl2 = fitlm([pd_vec_local'],y_vec')

mdl3 = fitlm([pd_vec_mean' pd_vec_local'],y_vec') 



% mdl3 = fitlm([y_vec'],pd_vec_local')
% 
% mdl4 = fitlm([pd_vec_mean'],pd_vec_local')

%% %%%%%%%%%%%%%%%%%%%%% cross-correlations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_lags = 15;
cov_type = 'coeff';

% target locus
target_local_vs_r = NaN(length(master_input_output{1}),2*n_lags+1);
target_ctrl_vs_r = NaN(length(master_input_output{1}),2*n_lags+1);
target_local_vs_f = NaN(length(master_input_output{1}),2*n_lags+1);
target_ctrl_vs_f = NaN(length(master_input_output{1}),2*n_lags+1);

target_f_vs_f = NaN(length(master_input_output{1}),2*n_lags+1);

iter = 1;
for i = 1:length(master_input_output{1})
    % extract vectors
    target_qc_flags = master_input_output{1}(i).dt_filter_gap;
    target_ctrl_protein_vec = master_input_output{1}(i).serial_null_protein_vec;
    target_spot_protein_vec = master_input_output{1}(i).spot_protein_vec;
%     target_nuclear_protein_vec = [master_input_output{1}.nuclear_protein_vec];
    target_fluo_vec = master_input_output{1}(i).fluo;
    target_r_vec = master_input_output{1}(i).r_vec_soft';
    
    % find contiguous segments that are sufficiently long to use for xcov
    % calculation
    flag_indices = unique([1 find(target_qc_flags) length(target_qc_flags)]);
    long_stretches = find(diff(flag_indices)-1>n_lags);
    for j = 1:length(long_stretches)
        indices = flag_indices(long_stretches(j))+1:flag_indices(long_stretches(j)+1)-1;
        target_local_vs_r(iter,:) = xcov(target_spot_protein_vec(indices),target_r_vec(indices),n_lags,cov_type);
        target_ctrl_vs_r(iter,:) = xcov(target_ctrl_protein_vec(indices),target_r_vec(indices),n_lags,cov_type);
        
        target_local_vs_f(iter,:) = xcov(target_spot_protein_vec(indices),target_fluo_vec(indices),n_lags,cov_type);
        target_ctrl_vs_f(iter,:) = xcov(target_ctrl_protein_vec(indices),target_fluo_vec(indices),n_lags,cov_type);
        
        target_f_vs_f(iter,:) = xcov(target_fluo_vec(indices),target_fluo_vec(indices),n_lags,cov_type);
        
        iter = iter + 1;
    end
end

% target locus
control_local_vs_r = NaN(length(master_input_output{2}),2*n_lags+1);
control_local_vs_f = NaN(length(master_input_output{2}),2*n_lags+1);

iter = 1;
for i = 1:length(master_input_output{2})
    % extract vectors
    control_qc_flags = master_input_output{2}(i).dt_filter_gap;    
    control_spot_protein_vec = master_input_output{2}(i).spot_protein_vec;
    control_fluo_vec = master_input_output{2}(i).fluo;
    control_r_vec = master_input_output{2}(i).r_vec_soft;
    
    % find contiguous segments that are sufficiently long to use for xcov
    % calculation
    flag_indices = unique([1 find(control_qc_flags) length(control_qc_flags)]);
    long_stretches = find(diff(flag_indices)-1>n_lags);
    for j = 1:length(long_stretches)
        indices = flag_indices(long_stretches(j))+1:flag_indices(long_stretches(j)+1)-1;
        control_local_vs_r(iter,:) = xcov(control_spot_protein_vec(indices),control_r_vec(indices),n_lags,cov_type);                
        control_local_vs_f(iter,:) = xcov(control_spot_protein_vec(indices),control_fluo_vec(indices),n_lags,cov_type);        
        
        iter = iter + 1;
    end
end

% get mean and standard deviation
t_fluo_boot = bootstrp(100,@mean,target_local_vs_f);
t_fluo_mean = mean(t_fluo_boot);
t_fluo_ub =  t_fluo_mean + nanstd(t_fluo_boot);
t_fluo_lb = t_fluo_mean - nanstd(t_fluo_boot);

t_r_boot = bootstrp(100,@mean,target_local_vs_r);
t_r_mean = nanmean(t_r_boot);
t_r_ub = t_r_mean + nanstd(t_r_boot);
t_r_lb = t_r_mean - nanstd(t_r_boot);

c_fluo_boot = bootstrp(100,@nanmean,target_ctrl_vs_f);
c_fluo_mean = nanmean(c_fluo_boot);
c_fluo_ub = c_fluo_mean + nanstd(c_fluo_boot);
c_fluo_lb = c_fluo_mean - nanstd(c_fluo_boot);

c_r_boot = bootstrp(100,@nanmean,target_ctrl_vs_r);
c_r_mean = nanmean(c_r_boot);
c_r_ub = c_r_mean + nanstd(c_r_boot);
c_r_lb = c_r_mean - nanstd(c_r_boot);

cc_fluo_boot = bootstrp(100,@mean,control_local_vs_f);
cc_fluo_mean = nanmean(cc_fluo_boot);
cc_fluo_ub = cc_fluo_mean + nanstd(cc_fluo_boot);
cc_fluo_lb = cc_fluo_mean - nanstd(cc_fluo_boot);

cc_r_boot = bootstrp(100,@mean,control_local_vs_r);
cc_r_mean = nanmean(cc_r_boot);
cc_r_ub = cc_r_mean + nanstd(cc_r_boot);
cc_r_lb = cc_r_mean - nanstd(cc_r_boot);

lag_vec = (-n_lags:n_lags)*20/60;



e_fig5 = figure;
hold on

fill([lag_vec fliplr(lag_vec)],[t_fluo_ub fliplr(t_fluo_lb)],cmap1(2,:),'FaceAlpha',.3,'EdgeAlpha',0)
fill([lag_vec fliplr(lag_vec)],[c_fluo_ub fliplr(c_fluo_lb)],cmap1(3,:),'FaceAlpha',.3,'EdgeAlpha',0)
fill([lag_vec fliplr(lag_vec)],[cc_fluo_ub fliplr(cc_fluo_lb)],cmap1(5,:),'FaceAlpha',.3,'EdgeAlpha',0)

% target fluo
p1 = plot(lag_vec,t_fluo_mean,'-','Color',cmap1(2,:),'LineWidth',2);

% control spot
p2 = plot(lag_vec,c_fluo_mean,'-','Color',cmap1(3,:),'LineWidth',2);

% control locus
p3 = plot(lag_vec,cc_fluo_mean,'-','Color',cmap1(5,:),'LineWidth',2);

legend([p1 p3 p2],'Dl at {\it snail} locus', 'Dl at {\it hbP2P}','virtual spot', 'Location','southwest');
ax = gca;
ax.YColor = 'black';%cmap1(2,:);
% grid on
xlabel('offset (minutes)')
ylabel('cross-covariance (fluo vs. local protein)')
set(gca,'Fontsize',14,'xtick',-5:5)

set(gca,    'Box','off',...
            'Color',[228,221,209]/255,...            
            'TickLength',[0.02,0.05])    
e_fig5.Color = 'white';        
e_fig5.InvertHardcopy = 'off';

% save
saveas(e_fig5,[FigurePath 'fluo_xcorr.tif'])
saveas(e_fig5,[FigurePath 'fluo_xcorr.pdf'])


e_fig6 = figure;
hold on

fill([lag_vec fliplr(lag_vec)],[t_r_ub fliplr(t_r_lb)],cmap1(2,:),'FaceAlpha',.3,'EdgeAlpha',0)
fill([lag_vec fliplr(lag_vec)],[c_r_ub fliplr(c_r_lb)],cmap1(3,:),'FaceAlpha',.3,'EdgeAlpha',0)
fill([lag_vec fliplr(lag_vec)],[cc_r_ub fliplr(cc_r_lb)],cmap1(5,:),'FaceAlpha',.3,'EdgeAlpha',0)

% target fluo
p1 = plot(lag_vec,t_r_mean,'-','Color',cmap1(2,:),'LineWidth',2);

% control spot
p2 = plot(lag_vec,c_r_mean,'-','Color',cmap1(3,:),'LineWidth',2);

% control locus
p3 = plot(lag_vec,cc_r_mean,'-','Color',cmap1(5,:),'LineWidth',2);

legend([p1 p3 p2],'Dl at {\it snail} locus', 'Dl at {\it hbP2P}','virtual spot', 'Location','southwest');
ax = gca;
ax.YColor = 'black';%cmap1(2,:);
% grid on
xlabel('offset (minutes)')
ylabel('cross-covariance (r vs. local protein)')
set(gca,'Fontsize',14,'xtick',-5:5)

set(gca,    'Box','off',...
            'Color',[228,221,209]/255,...            
            'TickLength',[0.02,0.05])    
e_fig6.Color = 'white';        
e_fig6.InvertHardcopy = 'off';

% save
saveas(e_fig6,[FigurePath 'r_xcorr.tif'])
saveas(e_fig6,[FigurePath 'r_xcorr.pdf'])


%%

mdl_a = fitlm(spot_target',r_target')
mdl_b = fitlm(control_target',r_target')
mdl_c = fitlm([control_target' spot_target'],r_target')
mdl_d = fitlm([imgaussfilt(-control_target' + spot_target',3)],imgaussfilt(r_target',3))