clc
clear
close all


%% calculate background
% load data

load('S:\Jake\Dropbox\OptoDropbox\2021-07-12-optoknirps_vasa-eYFP_background_snapshot\CompiledNuclei.mat')
mean(MeanVectorAll)

% calculate background
eYFP_mean = mean(MeanVectorAll)

%% plot background trend from another dataset
load('S:\Jake\Dropbox\OptoDropbox\2021-07-12-optoknirps_vasa-eYFP_background_embryo1\CompiledNuclei.mat')

fig = figure;
hold on
errorbar(APbinID*100,MeanVectorAP(31,:),SDVectorAP(31,:),'- .','MarkerSize',20)
%plot(APbinID*100,MeanVectorAP(31,:),'- .','MarkerSize',20)
errorbar(APbinID*100,MeanVectorAP(60,:),SDVectorAP(31,:),'- .','MarkerSize',20)
%plot(APbinID*100,MeanVectorAP(60,:),'- .','MarkerSize',20)
errorbar(APbinID*100,MeanVectorAP(90,:),SDVectorAP(31,:),'- .','MarkerSize',20)
%plot(APbinID*100,MeanVectorAP(90,:),'- .','MarkerSize',20)

legend('10 min','20 min','30 min')
xlabel('AP position (% embryo length)')
xlim([40 80])
ylim([0 6E5])
pbaspect([3 2 1])


