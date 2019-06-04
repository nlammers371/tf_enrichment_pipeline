clear
close all;

hbProject = 'Dl-Ven_hbP2P-mCh';
snaProject = 'Dl-Ven_snaBAC-mCh';
twoSpotProject = '_Archive\Dl_Venus_hbP2PsnaBAC_MCPmCherry';
dropboxFolder = 'E:\Nick\Dropbox (Garcia Lab)\';
rawPath = 'E:\LocalEnrichment\Data\PreProcessedData\';
proteinChannel = 1;

twoSpotDataPath = [dropboxFolder '\ProcessedEnrichmentData\' twoSpotProject '/'];
writePath = 'E:\Nick\Dropbox (Garcia Lab)\LocalEnrichmentFigures\twoSpot\';

%% Load hb, sna, and twoSpot datasets
% Filter for qc_flag & against unique nuclei
load([dropboxFolder '\ProcessedEnrichmentData\' hbProject '\nucleus_struct.mat'],'nucleus_struct')
nucleus_struct_filtered_hb = nucleus_struct([nucleus_struct.qc_flag]==1);
[nucleus_struct_filtered_hb.locus_identity] = deal('hbP2P');


load([dropboxFolder '\ProcessedEnrichmentData\' snaProject '\nucleus_struct.mat'],'nucleus_struct')
nucleus_struct_filtered_sna = nucleus_struct([nucleus_struct.qc_flag]==1);
[nucleus_struct_filtered_sna.locus_identity] = deal('snaBAC');

load([twoSpotDataPath '\nucleus_struct.mat'],'nucleus_struct')
nucleus_struct_filtered_twoSpot = nucleus_struct([nucleus_struct.qc_flag]==1);
[nucleus_struct_filtered_twoSpot.locus_identity] = deal('twoSpot');

%% Make cell arrays containing the identities of the loci to be used as
% nominal variables
nucleus_struct_locusID_hb = {nucleus_struct_filtered_hb.locus_identity}';
nucleus_struct_locusID_sna = {nucleus_struct_filtered_sna.locus_identity}';

% Calculate metrics for use in the regression
[nucleus_struct_metrics_hb, tOff_hb, fluoStartEndRatio_hb, fluoMean_hb, ...
    fluoVar_hb, fluoMax_hb, index_hb, fluoMaxTime_hb, fluoCenterMass_hb... %, mhmmLogLDiff_hb
    ] = calcRegressionMetrics(nucleus_struct_filtered_hb);
[nucleus_struct_metrics_sna, tOff_sna, fluoStartEndRatio_sna, fluoMean_sna, ...
    fluoVar_sna, fluoMax_sna, index_sna, fluoMaxTime_sna, ...
    fluoCenterMass_sna] = ... %, mhmmLogLDiff_sna
    calcRegressionMetrics(nucleus_struct_filtered_sna);
% [nucleus_struct_metrics_twoSpot, tOff_twoSpot, fluoStartEndRatio_twoSpot,...
%     fluoMean_twoSpot, fluoVar_twoSpot, fluoMax_twoSpot, index_twoSpot, ...
%     fluoMaxTime_twoSpot, fluoCenterMass_twoSpot, mhmmLogLDiff_twoSpot] = ...
%     calcRegressionMetrics(nucleus_struct_filtered_twoSpot);

% Trim & combine metrics to make training and testing datasets
hbTrainDataLength = int16(numel(nucleus_struct_filtered_hb)/2) - 1;
snaTrainDataLength = int16(numel(nucleus_struct_filtered_sna)/2) - 1;
hbTestDataLength = numel(nucleus_struct_filtered_hb) - hbTrainDataLength;
snaTestDataLength = numel(nucleus_struct_filtered_sna) - snaTrainDataLength;
% Ensure hb & sna datasets are the same size - replicate the hb training
% data until it matches the amount of sna training data
hbsnaTrainRemain = mod(snaTrainDataLength,hbTrainDataLength);
hbsnaTrainQuotient = (snaTrainDataLength-hbsnaTrainRemain)/hbTrainDataLength;
nucleus_struct_metrics_train_hb = nucleus_struct_metrics_hb(1:hbTrainDataLength,:);
nucleus_struct_locusID_train_hb = nucleus_struct_locusID_hb(1:hbTrainDataLength,:);
for i = 1:(hbsnaTrainQuotient-1)
    nucleus_struct_metrics_train_hb = cat(1, nucleus_struct_metrics_train_hb,nucleus_struct_metrics_hb(1:hbTrainDataLength,:));
    nucleus_struct_locusID_train_hb = cat(1, nucleus_struct_locusID_train_hb,nucleus_struct_locusID_train_hb(1:hbTrainDataLength,:));
end
if hbsnaTrainRemain ~= 0
    nucleus_struct_metrics_train_hb = cat(1, nucleus_struct_metrics_train_hb,nucleus_struct_metrics_hb(1:hbsnaTrainRemain,:));
    nucleus_struct_locusID_train_hb = cat(1, nucleus_struct_locusID_train_hb,nucleus_struct_locusID_hb(1:hbsnaTrainRemain,:));
end

% Create unpaired training & testing structures
nucleus_struct_metrics_train = cat(1, nucleus_struct_metrics_train_hb, nucleus_struct_metrics_sna(1:snaTrainDataLength,:)); 
nucleus_struct_metrics_test = cat(1, nucleus_struct_metrics_hb(hbTrainDataLength+1:end,:), nucleus_struct_metrics_sna(snaTrainDataLength+1:end,:)); 
nucleus_struct_locusID_train = cat(1, nucleus_struct_locusID_train_hb, nucleus_struct_locusID_sna(1:snaTrainDataLength)); 
nucleus_struct_locusID_test = cat(1, nucleus_struct_locusID_hb(hbTrainDataLength+1:end), nucleus_struct_locusID_sna(snaTrainDataLength+1:end)); 

% Create paired training & testing structures
pairedTrainSplit = int16(snaTrainDataLength/2) - 1;
metrics_spotAisHb_train = nucleus_struct_metrics_train_hb(1:pairedTrainSplit,:) - nucleus_struct_metrics_sna(1:pairedTrainSplit,:);
metrics_spotBisHb_train = nucleus_struct_metrics_sna(pairedTrainSplit+1:snaTrainDataLength,:) - nucleus_struct_metrics_train_hb(pairedTrainSplit+1:snaTrainDataLength,:);
metrics_paired_train = cat(1, metrics_spotAisHb_train, metrics_spotBisHb_train);
locusID_paired_train(1:pairedTrainSplit,1) = {'Spot A is hbP2P'};
locusID_paired_train(pairedTrainSplit+1:snaTrainDataLength,1) = {'Spot B is hbP2P'};
pairedTestLength = min(hbTestDataLength, snaTestDataLength);
pairedTestSplit = int16(pairedTestLength/2) - 1;
metrics_spotAisHb_test = nucleus_struct_metrics_hb(hbTrainDataLength+1:(hbTrainDataLength+1+pairedTestSplit),:)...
                            - nucleus_struct_metrics_sna(snaTrainDataLength+1:(snaTrainDataLength+1+pairedTestSplit),:);
metrics_spotBisHb_test = nucleus_struct_metrics_sna((snaTrainDataLength+1+pairedTestSplit+1):snaTrainDataLength+pairedTestLength,:)...
                            - nucleus_struct_metrics_hb((hbTrainDataLength+1+pairedTestSplit+1):hbTrainDataLength+pairedTestLength,:);
metrics_paired_test = cat(1, metrics_spotAisHb_test, metrics_spotBisHb_test);
locusID_paired_test(1:pairedTestSplit,1) = {'Spot A is hbP2P'};
locusID_paired_test(pairedTestSplit+1:pairedTestLength,1) = {'Spot B is hbP2P'};

% Full datasets for unpaired twoSpot analysis
nucleus_struct_metrics_hbsna = cat(1, nucleus_struct_metrics_hb, nucleus_struct_metrics_sna);
nucleus_struct_locusID_hbsna = cat(1, nucleus_struct_locusID_hb, nucleus_struct_locusID_sna); 

% Full datasets for paired twoSpot analysis

% Plot histograms of the metrics
plotMetrics('tOff', tOff_hb, tOff_sna, ...
    'Time locus turns off (sec from start of nc14)', writePath);
plotMetrics('fluoMaxTime', fluoMaxTime_hb, fluoMaxTime_sna, ...
    'Time of maximum locus fluorescence (sec from start of nc14)', writePath)
plotMetrics('fluoCenterMass', fluoCenterMass_hb, fluoCenterMass_sna, ...
    'Center of mass of fluorescence trace (sec from start of nc14)', ...
    writePath)
plotMetrics('fluoMean', fluoMean_hb, fluoMean_sna, ...
    'Mean locus fluorescence (au)', writePath)
plotMetrics('fluoVar', fluoVar_hb, fluoVar_sna, ...
    'Variance of locus fluorescence (au)', writePath)
% plotMetrics('mhmmLogLDiff', mhmmLogLDiff_hb, mhmmLogLDiff_sna, ...
%     'mHMM inference (au)', writePath)

%% Run multinomial logistic regression
% Unpaired testing
% 1 = hb, 2 = sna
locusID_train = nominal(nucleus_struct_locusID_train);
locusID_train = double(locusID_train);
[B_test,dev_test,stats_test] = mnrfit(nucleus_struct_metrics_train,locusID_train);
prob_test = mnrval(B_test,nucleus_struct_metrics_test);
% Check which parameters were weighted most
B_testWeights = nan(numel(B_test),1);
B_testWeights(1) = B_test(1);
for i = 2:numel(B_testWeights)
    B_testWeights(i) = B_test(i) * nucleus_struct_metrics_train(1,(i-1));
end

% Paired testing
% 1 = Spot A is hbP2P, 2 = Spot B is hbP2P
locusID_paired_train = nominal(locusID_paired_train);
locusID_paired_train = double(locusID_paired_train);
[B_test_paired,dev_test_paired,stats_test_paired] = mnrfit(metrics_paired_train,locusID_paired_train);
prob_test_paired = mnrval(B_test_paired,metrics_paired_test);
% Check which parameters were weighted most
B_testWeights_paired = nan(numel(B_test_paired),1);
B_testWeights_paired(1) = B_test_paired(1);
for i = 2:numel(B_testWeights_paired)
    B_testWeights_paired(i) = B_test_paired(i) * metrics_paired_train(1,(i-1));
end

% % twoSpot Analysis
% % Train regression
% locusID_hbsna = nominal(nucleus_struct_locusID_hbsna);
% locusID_hbsna = double(locusID_hbsna);
% [B_twoSpot,dev_twoSpot,stats_twoSpot] = mnrfit(nucleus_struct_metrics_hbsna,locusID_hbsna);
% % Fit twoSpot data
% prob_twoSpot = mnrval(B_twoSpot,nucleus_struct_metrics_twoSpot);

%% Plot results of testing the regression on known data
cm_plot = jet(128);
probBins = linspace(0,1,11);
% Unpaired
% Known hb loci
hbProbFig = figure;
hbProbCt = histc(prob_test(1:hbTestDataLength,1),probBins); % column 1 contains prob locus is hb
bar(probBins,hbProbCt,'FaceAlpha',.5,'FaceColor',cm_plot(20,:));
title('Particles known to be hb loci')
xlabel('Probability particle is hb locus')
ylabel('Number of particles')
saveas(hbProbFig,[writePath 'regression_test_hb_prob.png']);
% Known sna loci
snaProbFig = figure;
snaProbCt = histc(prob_test(hbTestDataLength+1:end,2),probBins);    % column 2 contains prob locus is sna
bar(probBins,snaProbCt,'FaceAlpha',.5,'FaceColor',cm_plot(60,:));
title('Particles known to be sna loci')
xlabel('Probability particle is sna locus')
ylabel('Number of particles')
saveas(snaProbFig,[writePath 'regression_test_sna_prob.png']);

% Paired
% Known Spot A is hbP2P
hbProbFigPaired = figure;
hbProbCt = histc(prob_test_paired(1:pairedTestSplit,1),probBins); % column 1 contains prob Spot A is hb
bar(probBins,hbProbCt,'FaceAlpha',.5,'FaceColor',cm_plot(20,:));
title('Pairs with Spot A known to be hbP2P')
xlabel('Probability Spot A is hb locus')
ylabel('Number of particles')
saveas(hbProbFigPaired,[writePath 'regression_test_hb_prob_paired.png']);
% Known Spot B is hbP2P
snaProbFigPaired = figure;
snaProbCt = histc(prob_test_paired(pairedTestSplit+1:end,2),probBins);    % column 2 contains prob spot B is sna
bar(probBins,snaProbCt,'FaceAlpha',.5,'FaceColor',cm_plot(60,:));
title('Pairs with Spot B known to be hbP2P')
xlabel('Probability Spot B is hb locus')
ylabel('Number of particles')
saveas(snaProbFigPaired,[writePath 'regression_test_sna_prob_paired.png']);

% % 2spot loci
% twoSpotProbFig = figure;
% twospotProbCt = histc(prob_twoSpot,probBins);
% bar(probBins,twospotProbCt,'FaceAlpha',.5,'FaceColor',cm_plot(120,:));
% title('Particles from 2spot dataset - ID unknown')
% xlabel('Probability particle is sna locus')
% ylabel('Number of particles')
% saveas(twoSpotProbFig,[writePath 'regression_twoSpot_prob.png']);


% ---------------- LOCAL FUNCTION calcRegressionMetrics -----------------
function [nucleus_struct_metrics, tOff, fluoStartEndRatio, fluoMean, ...
          fluoVar, fluoMax, index, fluoMaxTime, fluoCenterMass, ...
          mhmmLogLDiff] = calcRegressionMetrics(nucleus_struct)
    % Calculate the desired metrics to be used in the logistic regression

    % initialize all the arrays to their correct sizes
    tOff = nan(numel(nucleus_struct),1);
    fluoStartEndRatio = nan(numel(nucleus_struct),1);
    fluoMean = nan(numel(nucleus_struct),1);
    fluoVar = nan(numel(nucleus_struct),1);
    fluoMax = nan(numel(nucleus_struct),1);
    index = nan(numel(nucleus_struct),1);
    fluoMaxTime = nan(numel(nucleus_struct),1);
    fluoCenterMass = nan(numel(nucleus_struct),1);
%     mhmmLogLDiff = nan(numel(nucleus_struct),1);
    
    
    numMetrics = 5;
    nucleus_struct_metrics = nan(numel(nucleus_struct),numMetrics);

    for i = 1:numel(nucleus_struct)
        time = nucleus_struct(i).time;
        fluorescenceTrace = nucleus_struct(i).fluo;
        
        fluoDatapoints = find(~isnan(fluorescenceTrace));
%         if (fluoDatapoints(end) - fluoDatapoints(end-1)) > 5
%             lastDatapoint = fluoDatapoints(end-1);
%             tOff(i,1) = time(lastDatapoint);
%         else
            lastDatapoint = fluoDatapoints(end);
            tOff(i,1) = time(lastDatapoint);
%         end
        
        fluoStartEndRatio(i,1) = fluorescenceTrace(1) / fluorescenceTrace(end);
        fluoMean(i,1) = nanmean(fluorescenceTrace);
        fluoVar(i,1) = nanvar(fluorescenceTrace);
        
        [fluoMax(i,1), index(i,1)] = nanmax(fluorescenceTrace);
        fluoMaxTime(i,1) = time(index(i,1));
        
        fluoCenterMass(i,1) = nansum(fluorescenceTrace .* time) / nansum(fluorescenceTrace);
%         mhmmLogLDiff(i,1) = nucleus_struct(i).sna_logL - nucleus_struct(i).hbP2P_logL;
%         p = polyfit(time,fluo,4);
    end  

    % Add metrics to a double matrix such that each column is a different
    % metric and each row is the corresponding particle in
    % nucleus_struct_protein
    nucleus_struct_metrics(:,1) = fluoMaxTime;
    nucleus_struct_metrics(:,2) = tOff;
    nucleus_struct_metrics(:,3) = fluoCenterMass;
    nucleus_struct_metrics(:,4) = fluoMean;
    nucleus_struct_metrics(:,5) = fluoVar;
%     nucleus_struct_metrics(:,6) = mhmmLogLDiff;
%     nucleus_struct_metrics(:,6) = fluoStartEndRatio;

end

% ---------------- LOCAL FUNCTION plotMetrics -----------------
function plotMetrics(metricName, metric_hb, metric_sna, xlabelString, writePath)
    cm_plot = jet(128);    
    bins_metric = linspace(min([metric_hb ; metric_sna]),max([metric_hb ; metric_sna]),10);
    metric_hb_ct = histc(metric_hb, bins_metric);
    metric_sna_ct = histc(metric_sna, bins_metric);

    metricFigure = figure;
    hold on
    b_sna = bar(bins_metric,metric_sna_ct,'FaceAlpha',.5,'FaceColor',cm_plot(60,:));
    b_hb = bar(bins_metric,metric_hb_ct,'FaceAlpha',.5,'FaceColor',cm_plot(20,:));
    title(['hb vs sna loci metrics: ', metricName])
    xlabel(xlabelString)
    ylabel('Number of loci')
    legend ([b_sna b_hb], 'sna','hb')
    hold off
    saveas(metricFigure,[writePath, metricName, '_hbVsSna.png']);
end
