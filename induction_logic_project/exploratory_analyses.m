% script to explore trends in opto data
clear
close all

% set basic paths
DataRoot = 'C:\Users\nlamm\Dropbox (Personal)\InductionLogic\';
% DataRoot = 'S:\Nick\Dropbox\InductionLogic\';
project = '20201117_v2';
% project = '20200807_opto_chronic';
DataPath = [DataRoot project filesep];
% make figure path
FigurePath = [DataPath 'exploratory_figures' filesep];
mkdir(FigurePath);

% get list of data sets
SpotsFile = dir([DataPath '*spots_only.xlsx']);
ProteinFile = dir([DataPath '*protein_only.xlsx']);
SpotsFileWT = dir([DataPath '*WT.xlsx']);

% specify time res
dT = 30;

% load spots dataset (WT)
raw_spots_table_wt = readtable([DataPath SpotsFileWT(1).name]);
raw_spots_array_wt = raw_spots_table_wt{:,:};

% load spots dataset
raw_spots_table = readtable([DataPath SpotsFile(1).name]);
raw_spots_array = raw_spots_table{:,:};
time_vec = (0:dT:dT*size(raw_spots_array,1)-1)/60;

% load protein dataset
raw_protein_table = readtable([DataPath ProteinFile(1).name]);
raw_protein_array = raw_protein_table{:,:};

% It looks like there are 10 spot frames for every protein frame
proteinFrames = 1:10:241;
raw_protein_array_interp = interp1(proteinFrames,raw_protein_array,proteinFrames(1):proteinFrames(end));

%% %%%%%% Assess significance of shift during time period of interest %%%%%

protein_figure = figure; 
plot(time_vec,nanmean(raw_protein_array_interp,2));
xlabel('time (minutes)')
ylabel('[YAP] (au)')
set(gca,'Fontsize',14)
grid on

%%

% set sampling parameters
windowSize = 10;
pauseSize = 5;
trueCenterIndex = find(time_vec==23.5);

% calculate percent change at perturbation point
preFluo = mean(mean(raw_spots_array(trueCenterIndex-pauseSize-windowSize:trueCenterIndex-pauseSize,:),2));
postFluo = mean(mean(raw_spots_array(trueCenterIndex+pauseSize:trueCenterIndex+pauseSize+windowSize,:),2));
actualDelta = (postFluo-preFluo)/preFluo;

% now randomly sample WT changepoints
nTracesPerSample = size(raw_spots_array,2);
nBoots = 1e3;

% falseCenterOptions = windowSize+pauseSize+1:length(time_vec)-windowSize-pauseSize;
traceOptions = 1:size(raw_spots_array_wt,2);
falseDeltas = NaN(1,nBoots);
falseCenterTimes = NaN(1,nBoots);

trueCenterOptions = [trueCenterIndex trueCenterIndex];
trueDeltas = NaN(1,nBoots);
trueCenterTimes = NaN(1,nBoots);

for n = 1:nBoots
    % sample time
    falseCenterIndex = randsample(trueCenterOptions,1);
    falseCenterTimes(n) = time_vec(falseCenterIndex);
    % sample traces 
    falseBootTraces = raw_spots_array_wt(:,randsample(traceOptions,nTracesPerSample,true));
    % calculate metrics
    preFluo = mean(mean(falseBootTraces(falseCenterIndex-pauseSize-windowSize:falseCenterIndex-pauseSize,:),2));
    postFluo = mean(mean(falseBootTraces(falseCenterIndex+pauseSize:falseCenterIndex+pauseSize+windowSize,:),2));
    falseDeltas(n) = 100*(postFluo-preFluo)/preFluo;
    
    % sample time
    trueCenterIndex = randsample(trueCenterOptions,1);
    trueCenterTimes(n) = time_vec(trueCenterIndex);
    % sample traces 
    trueBootTraces = raw_spots_array(:,randsample(1:nTracesPerSample,nTracesPerSample,true));
    % calculate metrics
    preFluo = mean(mean(trueBootTraces(trueCenterIndex-pauseSize-windowSize:trueCenterIndex-pauseSize,:),2));
    postFluo = mean(mean(trueBootTraces(trueCenterIndex+pauseSize:trueCenterIndex+pauseSize+windowSize,:),2));
    trueDeltas(n) = 100*(postFluo-preFluo)/preFluo;
    
end

% make figure
close all
deltaFig = figure;
cmap = brewermap(9,'Paired');
hold on

deltaBins = 100*linspace(-.5,.5,50);

histogram(trueDeltas,deltaBins,'FaceColor',cmap(3,:),'Normalization','probability')
histogram(falseDeltas,deltaBins,'FaceColor',cmap(2,:),'Normalization','probability')

xlabel('percent change')
ylabel('share of total')
set(gca,'Fontsize',14)
grid on
legend('actual perturbation time','WT control')
set(gca,'xtick',-50:10:50)

saveas(deltaFig,[FigurePath 'opto_vs_ctrl_hist.png'])


scatter_fig = figure;
hold on
% OPTO
scatter(1.5+rand(size(trueDeltas))*.25-.125, trueDeltas, 'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k',...
              'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',0)
            
errorbar(1.5,mean(trueDeltas),std(trueDeltas),'o','Color','k','LineWidth',1.5)

sTrue = scatter(1.5, mean(trueDeltas),75 , 'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k',...
        'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);

      
% CONTROL      
scatter(1+rand(size(falseDeltas))*.25-.125, falseDeltas, 'MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k',...
        'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',0)
errorbar(1,mean(falseDeltas),std(falseDeltas),'o','Color','k','LineWidth',1.5)
sFalse = scatter(1, mean(falseDeltas),75, 'MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k',...
        'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);
      
      
set(gca,'xtick',[1 1.5],'xticklabels',{'wt control', 'opto response'})      
grid on
ylabel('percent change')
set(gca,'Fontsize',14)
xlim([.75 1.75])
saveas(scatter_fig,[FigurePath 'opto_vs_ctrl_scatter.png'])
% %%
% pegTime = 150;
% window = 45;
% beforeIndices = pegTime-window:pegTime-5;
% afterIndices = pegTime+5:pegTime+window;
% 
% before_fluo = nanmean(raw_spots_array(beforeIndices,:));
% after_fluo = nanmean(raw_spots_array(afterIndices,:));
% 
% before_fluo_wt = nanmean(raw_spots_array_wt(beforeIndices,:));
% after_fluo_wt = nanmean(raw_spots_array_wt(afterIndices,:));
% % 
% before_protein = nanmean(raw_protein_array_interp(beforeIndices,:));
% after_protein = nanmean(raw_protein_array_interp(afterIndices,:));
% 
% testFilter = after_fluo == 0 & before_fluo ~= 0; 
% testFilterWT = after_fluo_wt == 0 & before_fluo_wt ~= 0; 
% sum(testFilter)
% sum(testFilterWT)
% 
% close all
% figure;
% scatter(before_protein(~testFilter),after_protein(~testFilter))
% hold on
% scatter(before_protein(testFilter),after_protein(testFilter),'filled')
% 
% figure;
% scatter(before_fluo,after_fluo)
