% script to plot results from cpHMM inference
clear
close all
addpath(genpath('utilities'))

projectNameCell = {'20200807_WT', '20200807_opto_chronic'};
DataRoot = 'S:\Nick\Dropbox\InductionLogic\';

% make figure path
FigPath = [DataRoot 'comparison_figures' filesep];
mkdir(FigPath);

% Set basic plotting parameters
xVar = 'time';

% useful colors
MarkerSize = 100;
blue = [115 143 193]/256;
purple = [171 133 172]/256;
red = [213 108 85]/256;

% set axes
dur_lims = [0 3.5];
freq_lims = [0 3];
init_lims = [0 18]*1e4;

% initialize plotting structure
plotting_struct = struct;

for p = 1:length(projectNameCell)
    % set project to analyze 
    projectName = projectNameCell{p};

    % get path to results    
    resultsDir = [DataRoot projectName filesep 'cpHMM_results' filesep];        
    
    % get list of projects
    resultList = dir([resultsDir '*result*']);            
    
    for r = 1%:length(resultList)
      
        % load data
        load([resultsDir filesep resultList(r).name]);
                  
        % transfer results
        plotting_struct(p).init_vec_mean = compiledResults.init_vec_mean;
        plotting_struct(p).init_vec_ste = compiledResults.init_vec_ste;
        
        plotting_struct(p).freq_vec_mean = compiledResults.freq_vec_mean;
        plotting_struct(p).freq_vec_ste = compiledResults.freq_vec_ste;
        
        plotting_struct(p).dur_vec_mean = compiledResults.dur_vec_mean;
        plotting_struct(p).dur_vec_ste = compiledResults.dur_vec_ste;
        
        plotting_struct(p).time_vec_mean = compiledResults.time_vec_mean;
        plotting_struct(p).time_vec_ste = compiledResults.time_vec_ste;
        
        plotting_struct(p).fluo_vec_mean = compiledResults.fluo_mean;
        plotting_struct(p).fluo_vec_ste = compiledResults.fluo_ste;
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make figures        
close all

r_trend = figure;
hm_cm = brewermap(9,'Set2');
colormap(hm_cm);

hold on
for i = 1:length(plotting_struct)  
    e = errorbar([1 2]+(i-1.5)*.05,plotting_struct(i).init_vec_mean,plotting_struct(i).init_vec_ste,'o--','Color','black','LineWidth',1);          
%     e.CapSize = 0;
end
s = [];
for i = 1:length(plotting_struct)    
    s(i) = scatter([1 2]+(i-1.5)*.05,plotting_struct(i).init_vec_mean,MarkerSize,'o','MarkerFaceColor',hm_cm(i+1,:),'MarkerEdgeColor','black');
end        
grid on
%         xlim(x_lim)
% ylim([50 95])
% xlabel('time') % NL: need to make this dynamic
set(gca,'xtick',[1 2]+(i-1.5)*.05,'xticklabel',{'before','after'})
ylabel('burst amplitude (au/min)')
legend(s,'Sox2 (wildtype control)','Sox2 (opto-chronic)','Location','southwest')
title(['Burst Amplitude (r)'])
set(gca,'Fontsize',14)
% ylim([2000 3000]);
xlim([.8 2.2])
box on
saveas(r_trend,[FigPath, 'burst_amp.tif'])
saveas(r_trend,[FigPath, 'burst_amp.pdf'])

%%
dur_trend = figure;
hm_cm = brewermap(9,'Set2');
colormap(hm_cm);

hold on
for i = 1:length(plotting_struct)  
    e = errorbar([1 2]+(i-1.5)*.05,plotting_struct(i).dur_vec_mean,plotting_struct(i).dur_vec_ste,'o--','Color','black','LineWidth',1);          
%     e.CapSize = 0;
end
s = [];
for i = 1:length(plotting_struct)    
    s(i) = scatter([1 2]+(i-1.5)*.05,plotting_struct(i).dur_vec_mean,MarkerSize,'o','MarkerFaceColor',hm_cm(i+1,:),'MarkerEdgeColor','black');
end        
grid on
%         xlim(x_lim)
% ylim([50 95])
% xlabel('time') % NL: need to make this dynamic
set(gca,'xtick',[1 2]+(i-1.5)*.05,'xticklabel',{'before','after'})
ylabel('burst duration (min)')
legend(s,'Sox2 (wildtype control)','Sox2 (opto-chronic)','Location','southwest')
title(['Burst Duration (1/k_{off})'])
set(gca,'Fontsize',14)
% ylim([2000 3000]);
xlim([.8 2.2])
box on
saveas(dur_trend,[FigPath, 'burst_dur.tif'])
saveas(dur_trend,[FigPath, 'burst_dur.pdf'])

%%

freq_trend = figure;

hold on
for i = 1:length(plotting_struct)  
    e = errorbar([1 2]+(i-1.5)*.05,plotting_struct(i).freq_vec_mean,plotting_struct(i).freq_vec_ste,'o--','Color','black','LineWidth',1);          
%     e.CapSize = 0;
end
s = [];
for i = 1:length(plotting_struct)    
    s(i) = scatter([1 2]+(i-1.5)*.05,plotting_struct(i).freq_vec_mean,MarkerSize,'o','MarkerFaceColor',hm_cm(i+1,:),'MarkerEdgeColor','black');
end        
grid on
%         xlim(x_lim)
% ylim([50 95])
% xlabel('time') % NL: need to make this dynamic
set(gca,'xtick',[1 2]+(i-1.5)*.05,'xticklabel',{'before','after'})
ylabel('burst frequency (1/min)')
legend(s,'Sox2 (wildtype control)','Sox2 (opto-chronic)','Location','southwest')
title(['Burst Frequency (k_{on})'])
set(gca,'Fontsize',14)
% ylim([2000 3000]);
xlim([.8 2.2])
box on
saveas(freq_trend,[FigPath, 'burst_freq.tif'])
saveas(freq_trend,[FigPath, 'burst_freq.pdf'])