% script to plot results from cpHMM inference
clear
close all
addpath(genpath('utilities'))

% get path
baseProject = 'NSv1';
liveProject = LiveEnrichmentProject(baseProject);
% make output paths for data and figures
slashesFig = strfind(liveProject.figurePath,'\');
FigurePath = [liveProject.figurePath(1:slashesFig(end-1)) 'eve2_recon_analyses' filesep 'cpHMM' filesep];
mkdir(FigurePath)

% projectNameCell = {'EveGtSL','EveGtSL-S1Null','EveWt','EveS1Null'};%};
projectNameCell = {'MSE-WT','NSv1','Rand1','Rand4'};%};
% resultsRoot = 'S:\Nick\Dropbox\InductionLogic\';
w = 7;
K = 3;

% initialize structure to store all results
master_struct = struct;

% loop through projects and compile results
for p = 1:length(projectNameCell)
    
    % set project to analyze 
    projectName = projectNameCell{p};

    % get path to results    
    liveProject = LiveEnrichmentProject(projectName);
    resultsDir = [liveProject.dataPath 'cpHMM_results' filesep];
    
    % get list of all inference subdirectories. By default, we'll generate
    % summaries for all non-empty inference sub-directory
    infDir = dir([resultsDir 'w' num2str(w) '*']);
    if length(infDir) > 1
        error('multiple inf directories found')
    end
    % iterate through the directories and compile the results 
    resultsPath = [infDir(1).folder filesep infDir(1).name filesep];

    % load inference options file 
    load([resultsPath 'inferenceOptions.mat'])
    if length(fieldnames(inferenceOptions)) == 1
      inferenceOptions = inferenceOptions.inferenceOptions;
    end

    % extract key inference parameters
    Tres = inferenceOptions.Tres;
    nStates = inferenceOptions.nStates;

    % load
    load([infDir(1).folder filesep 'compiledResults_' infDir(1).name '.mat'],'compiledResults')
    
    % add to structure
    master_struct(p).compiledResults = compiledResults;
    master_struct(p).projectName = projectNameCell{p};
end
    
%% Make plots
close all
stripe_center_ind = 2;
legend_str_short = {'MSE2','NS v1', 'Rand v1','Rand v4'};
% initiation rate
init_fig = figure;
cmap1 = brewermap([],'Set2');
hold on

% individual result scatters
for p = 1:length(master_struct)
    init_results = master_struct(p).compiledResults.init_results{stripe_center_ind}*60 / 1e4;    
    outlier_flags = master_struct(p).compiledResults.outlier_flags{stripe_center_ind};
    scatter(repelem(p,sum(~outlier_flags)),init_results(~outlier_flags),50,...
                'MarkerFaceColor',cmap1(p,:),'MarkerEdgeColor','k',...
                'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25);
end        

% means
for p = 1:length(master_struct)
    init_mean = master_struct(p).compiledResults.init_vec_mean(stripe_center_ind) / 1e4;    
    scatter(p,init_mean,75,'s',...
                'MarkerFaceColor',cmap1(p,:),'MarkerEdgeColor','k',...
                'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);
end

ylabel('initiation rate (AU per minute)')
% ylabel(y_label)
set(gca,'Fontsize',14);
grid on
xlim([0.5 4.5])
ylim([0 15])
% legend(p,legend_str{:},'Location',loc_str{param_index})

set(gca,'xtick',1:length(projectNameCell),'xticklabels',legend_str_short)
xtickangle(-30) 
set(gca,'Color',[228,221,209]/255) 

init_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(init_fig,[FigurePath 'init_scatter.png'])
saveas(init_fig,[FigurePath 'init_scatter.pdf'])

%% burst frequency
freq_fig = figure;
cmap1 = brewermap([],'Set2');
hold on

% individual result scatters
for p = 1:length(master_struct)
    freq_results = master_struct(p).compiledResults.freq_results{stripe_center_ind}*60;    
    outlier_flags = master_struct(p).compiledResults.outlier_flags{stripe_center_ind};
    scatter(repelem(p,sum(~outlier_flags)),freq_results(~outlier_flags),50,...
                'MarkerFaceColor',cmap1(p,:),'MarkerEdgeColor','k',...
                'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25);
end        

% means
for p = 1:length(master_struct)
    freq_mean = master_struct(p).compiledResults.freq_vec_mean(stripe_center_ind);    
    scatter(p,freq_mean,75,'s',...
                'MarkerFaceColor',cmap1(p,:),'MarkerEdgeColor','k',...
                'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);
end

ylabel('burst frequency (events per minute)')
% ylabel(y_label)
set(gca,'Fontsize',14);
grid on
xlim([0.5 4.5])
% ylim([0 15])
% legend(p,legend_str{:},'Location',loc_str{param_index})

set(gca,'xtick',1:length(projectNameCell),'xticklabels',legend_str_short)
xtickangle(-30) 
set(gca,'Color',[228,221,209]/255) 
ylim([0 1.5])
freq_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(freq_fig,[FigurePath 'freq_scatter.png'])
saveas(freq_fig,[FigurePath 'freq_scatter.pdf'])

%% burst duration
dur_fig = figure;
cmap1 = brewermap([],'Set2');
hold on

% individual result scatters
for p = 1:length(master_struct)
    dur_results = master_struct(p).compiledResults.dur_results{stripe_center_ind}/60;    
    outlier_flags = master_struct(p).compiledResults.outlier_flags{stripe_center_ind};
    scatter(repelem(p,sum(~outlier_flags)),dur_results(~outlier_flags),50,...
                'MarkerFaceColor',cmap1(p,:),'MarkerEdgeColor','k',...
                'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25);
end        

% means
for p = 1:length(master_struct)
    dur_mean = master_struct(p).compiledResults.dur_vec_mean(stripe_center_ind);    
    scatter(p,dur_mean,75,'s',...
                'MarkerFaceColor',cmap1(p,:),'MarkerEdgeColor','k',...
                'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);
end
            
ylabel('burst duration (minutes)')
% ylabel(y_label)
set(gca,'Fontsize',14);
grid on
xlim([0.5 4.5])
% ylim([0 15])
% legend(p,legend_str{:},'Location',loc_str{param_index})

set(gca,'xtick',1:length(projectNameCell),'xticklabels',legend_str_short)
xtickangle(-30) 
set(gca,'Color',[228,221,209]/255) 
ylim([0 1.5])
dur_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(dur_fig,[FigurePath 'dur_scatter.png'])
saveas(dur_fig,[FigurePath 'dur_scatter.pdf'])            