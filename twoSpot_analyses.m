clear

hbProject = 'Dl-Ven x hbP2P';
snaProject = 'Dl-Ven x snaBAC';
twoSpotProject = '_Archive\Dl_Venus_hbP2PsnaBAC_MCPmCherry';
keyword = 'Dl_Venus_hbP2PsnaBAC';
dropboxFolder = 'E:\Nick\Dropbox (Garcia Lab)\';
% includeVec = [1 3 4];
rawPath = 'E:\LocalEnrichment\Data\PreProcessedData\';
proteinChannel = 1;

twoSpotDataPath = [dropboxFolder '\ProcessedEnrichmentData\' twoSpotProject '/'];
writePath = 'E:\Nick\Dropbox (Garcia Lab)\LocalEnrichmentFigures\twoSpot\';

%% Load hb, sna, and twoSpot datasets
load([dropboxFolder '\ProcessedEnrichmentData\' hbProject '\nucleus_struct_protein.mat'],'nucleus_struct_protein')
nucleus_struct_protein_hb = nucleus_struct_protein;
[nucleus_struct_protein_hb.locus_identity] = deal('hbP2P');

load([dropboxFolder '\ProcessedEnrichmentData\' snaProject '\nucleus_struct_protein.mat'],'nucleus_struct_protein')
nucleus_struct_protein_sna = nucleus_struct_protein;
[nucleus_struct_protein_sna.locus_identity] = deal('snaBAC');

load([twoSpotDataPath '\nucleus_struct_protein.mat'],'nucleus_struct_protein')
nucleus_struct_protein_twoSpot = nucleus_struct_protein;
[nucleus_struct_protein_twoSpot.locus_identity] = deal('twoSpot');

%% Make cell arrays containing the identities of the loci to be used as
% nominal variables
nucleus_struct_locusID_hb = {nucleus_struct_protein_hb.locus_identity}';
nucleus_struct_locusID_sna = {nucleus_struct_protein_sna.locus_identity}';

% Calculate metrics for use in the regression
[nucleus_struct_metrics_hb, tOff_hb, fluoStartEndRatio_hb, fluoMean_hb, ...
    fluoVar_hb, fluoMax_hb, index_hb, fluoMaxTime_hb, fluoCenterMass_hb] = ...
    calcRegressionMetrics(nucleus_struct_protein_hb);
[nucleus_struct_metrics_sna, tOff_sna, fluoStartEndRatio_sna, fluoMean_sna, ...
    fluoVar_sna, fluoMax_sna, index_sna, fluoMaxTime_sna, fluoCenterMass_sna] = ... 
    calcRegressionMetrics(nucleus_struct_protein_sna);
[nucleus_struct_metrics_twoSpot, tOff_twoSpot, fluoStartEndRatio_twoSpot,...
    fluoMean_twoSpot, fluoVar_twoSpot, fluoMax_twoSpot, index_twoSpot, ...
    fluoMaxTime_twoSpot, fluoCenterMass_twoSpot] = ...
    calcRegressionMetrics(nucleus_struct_protein_twoSpot);

% Trim & combine metrics to make training and testing datasets
hbTrainDataLength = 300; % = numel(nucleus_struct_protein_hb);
snaTrainDataLength = 300; % = numel(nucleus_struct_protein_sna);
hbTestDataLength = numel(nucleus_struct_protein_hb) - hbTrainDataLength;
snaTestDataLength = numel(nucleus_struct_protein_sna) - snaTrainDataLength;

nucleus_struct_metrics_train = cat(1, nucleus_struct_metrics_hb(1:hbTrainDataLength,:), nucleus_struct_metrics_sna(1:snaTrainDataLength,:)); 
nucleus_struct_metrics_test = cat(1, nucleus_struct_metrics_hb(hbTrainDataLength+1:end,:), nucleus_struct_metrics_sna(snaTrainDataLength+1:end,:)); 
nucleus_struct_locusID_train = cat(1, nucleus_struct_locusID_hb(1:hbTrainDataLength), nucleus_struct_locusID_sna(1:snaTrainDataLength)); 
nucleus_struct_locusID_test = cat(1, nucleus_struct_locusID_hb(hbTrainDataLength+1:end), nucleus_struct_locusID_sna(snaTrainDataLength+1:end)); 


% Full datasets for twoSpot analysis
nucleus_struct_metrics_hbsna = cat(1, nucleus_struct_metrics_hb, nucleus_struct_metrics_sna);
nucleus_struct_locusID_hbsna = cat(1, nucleus_struct_locusID_hb, nucleus_struct_locusID_sna); 


%% Run multinomial logistic regression
% Testing
locusID_train = nominal(nucleus_struct_locusID_train);
locusID_train = double(locusID_train);
[B_test,dev_test,stats_test] = mnrfit(nucleus_struct_metrics_train,locusID_train);
prob_test = mnrval(B_test,nucleus_struct_metrics_test);

% twoSpot Analysis
% Train regression
locusID_hbsna = nominal(nucleus_struct_locusID_hbsna);
locusID_hbsna = double(locusID_hbsna);
[B_twoSpot,dev_twoSpot,stats_twoSpot] = mnrfit(nucleus_struct_metrics_hbsna,locusID_hbsna);
% Fit twoSpot data
prob_twoSpot = mnrval(B_twoSpot,nucleus_struct_metrics_twoSpot);

%% Figures
cm_plot = jet(128);

% Plot results of testing the regression on known data
probBins = linspace(0,1,11);

hbProbFig = figure;
hbProbCt = histc(prob_test(1:243,1),probBins);
bar(probBins,hbProbCt,'FaceAlpha',.5,'FaceColor',cm_plot(20,:));
title('Particles known to be hb loci')
xlabel('Probability particle is hb locus')
ylabel('Number of particles')
ylim([0 250])
saveas(hbProbFig,[writePath 'regression_test_hb_prob.png']);

snaProbFig = figure;
snaProbCt = histc(prob_test(244:end,2),probBins);
bar(probBins,snaProbCt,'FaceAlpha',.5,'FaceColor',cm_plot(60,:));
title('Particles known to be sna loci')
xlabel('Probability particle is sna locus')
ylabel('Number of particles')
saveas(snaProbFig,[writePath 'regression_test_sna_prob.png']);

twoSpotProbFig = figure;
twospotProbCt = histc(prob_twoSpot,probBins);
bar(probBins,twospotProbCt,'FaceAlpha',.5,'FaceColor',cm_plot(120,:));
title('Particles from 2spot dataset - ID unknown')
xlabel('Probability particle is sna locus')
ylabel('Number of particles')
saveas(twoSpotProbFig,[writePath 'regression_twoSpot_prob.png']);

%% Plot the metrics
bins_tOff = linspace(min([tOff_hb ; tOff_sna]),max([tOff_hb ; tOff_sna]),10);
tOff_hb_ct = histc(tOff_hb, bins_tOff);
tOff_sna_ct = histc(tOff_sna, bins_tOff);

tOffFig = figure;
hold on
b_sna = bar(bins_tOff,tOff_sna_ct,'FaceAlpha',.5,'FaceColor',cm_plot(60,:));
b_hb = bar(bins_tOff,tOff_hb_ct,'FaceAlpha',.5,'FaceColor',cm_plot(20,:));
title('hb vs sna loci metrics: toff')
xlabel('Time locus shuts off (sec from start of nc14)')
ylabel('Number of particles')
legend ([b_sna b_hb], 'sna','hb')
hold off
saveas(tOffFig,[writePath 'tOff_hbVsSna.png']);

bins_fluoMaxTime = linspace(min([fluoMaxTime_hb ; fluoMaxTime_sna]),max([fluoMaxTime_hb ; fluoMaxTime_sna]),10);
tOff_hb_ct = histc(fluoMaxTime_hb, bins_fluoMaxTime);
tOff_sna_ct = histc(fluoMaxTime_sna, bins_fluoMaxTime);

fluoMaxTimeFig = figure;
hold on
b_sna = bar(bins_fluoMaxTime,tOff_sna_ct,'FaceAlpha',.5,'FaceColor',cm_plot(60,:));
b_hb = bar(bins_fluoMaxTime,tOff_hb_ct,'FaceAlpha',.5,'FaceColor',cm_plot(20,:));
title('hb vs sna loci metrics: time of max fluorescence')
xlabel('Time of maximum locus fluorescence (sec from start of nc14)')
ylabel('Number of particles')
legend ([b_sna b_hb], 'sna','hb')
hold off
saveas(fluoMaxTimeFig,[writePath 'fluoMaxTime_hbVsSna.png']);

bins_fluoCenterMass = linspace(min([fluoCenterMass_hb ; fluoCenterMass_sna]),max([fluoCenterMass_hb ; fluoCenterMass_sna]),10);
fluoCentMass_hb_ct = histc(fluoCenterMass_hb, bins_fluoCenterMass);
fluoCentMass_sna_ct = histc(fluoCenterMass_sna, bins_fluoCenterMass);

fluoCentMassFig = figure;
hold on
b_sna = bar(bins_fluoCenterMass,fluoCentMass_sna_ct,'FaceAlpha',.5,'FaceColor',cm_plot(60,:));
b_hb = bar(bins_fluoCenterMass,fluoCentMass_hb_ct,'FaceAlpha',.5,'FaceColor',cm_plot(20,:));
title('hb vs sna loci metrics: Fluorescence profile')
xlabel('Center of mass of fluorescence (sec from start of nc14)')
ylabel('Number of particles')
legend ([b_sna b_hb], 'sna','hb')
hold off
saveas(fluoCentMassFig,[writePath 'fluoCenterMass_hbVsSna.png']);

bins_fluoMean = linspace(min([fluoMean_hb ; fluoMean_sna]),max([fluoMean_hb ; fluoMean_sna]),10);
fluoMean_hb_ct = histc(fluoMean_hb, bins_fluoMean);
fluoMean_sna_ct = histc(fluoMean_sna, bins_fluoMean);

fluoMeanFig = figure;
hold on
b_sna = bar(bins_fluoMean,fluoMean_sna_ct,'FaceAlpha',.5,'FaceColor',cm_plot(60,:));
b_hb = bar(bins_fluoMean,fluoMean_hb_ct,'FaceAlpha',.5,'FaceColor',cm_plot(20,:));
title('hb vs sna loci metrics: Mean spot fluorescence')
xlabel('Mean fluorescence (au)')
ylabel('Number of particles')
legend ([b_sna b_hb], 'sna','hb')
hold off
saveas(fluoMeanFig,[writePath 'fluoMean_hbVsSna.png']);

bins_fluoVar = linspace(min([fluoVar_hb ; fluoVar_sna]),max([fluoVar_hb ; fluoVar_sna]),10);
fluoVar_hb_ct = histc(fluoVar_hb, bins_fluoVar);
fluoVar_sna_ct = histc(fluoVar_sna, bins_fluoVar);

fluoVarFig = figure;
hold on
b_sna = bar(bins_fluoVar,fluoVar_sna_ct,'FaceAlpha',.5,'FaceColor',cm_plot(60,:));
b_hb = bar(bins_fluoVar,fluoVar_hb_ct,'FaceAlpha',.5,'FaceColor',cm_plot(20,:));
title('hb vs sna loci metrics: Variance of spot fluorescence')
xlabel('Variance of spot fluorescence (au)')
ylabel('Number of particles')
legend ([b_sna b_hb], 'sna','hb')
hold off
saveas(fluoVarFig,[writePath 'fluoVar_hbVsSna.png']);

% ---------------- LOCAL FUNCTION calcRegressionMetrics -----------------
function [nucleus_struct_metrics, tOff, fluoStartEndRatio, fluoMean, ...
          fluoVar, fluoMax, index, fluoMaxTime, fluoCenterMass] = ...
          calcRegressionMetrics(nucleus_struct_protein)

    numMetrics = 5;

    % Calculate the desired metrics to be used in the logistic regression    
    tOff = nan(numel(nucleus_struct_protein),1);
    fluoStartEndRatio = nan(numel(nucleus_struct_protein),1);
    fluoMean = nan(numel(nucleus_struct_protein),1);
    fluoVar = nan(numel(nucleus_struct_protein),1);
    fluoMax = nan(numel(nucleus_struct_protein),1);
    index = nan(numel(nucleus_struct_protein),1);
    fluoMaxTime = nan(numel(nucleus_struct_protein),1);
    fluoCenterMass = nan(numel(nucleus_struct_protein),1);
    nucleus_struct_metrics = nan(numel(nucleus_struct_protein),numMetrics);

    for i = 1:numel(nucleus_struct_protein)
        if (nucleus_struct_protein(i).frames > 9) & (sum(~isnan(nucleus_struct_protein(i).fluo)) > 40)
            tOff(i,1) = nucleus_struct_protein(i).time(end);
            fluo = nucleus_struct_protein(i).fluo;
            fluoStartEndRatio(i,1) = fluo(1) / fluo(end);
            fluoMean(i,1) = mean(~isnan(fluo));
            fluoVar(i,1) = var(~isnan(fluo));
            [fluoMax(i,1), index(i,1)] = max(~isnan(fluo));
            time = nucleus_struct_protein(i).time;
            fluoMaxTime(i,1) = time(index(i,1));
            fluoCenterMass(i,1) = nansum(fluo .* time) / nansum(fluo);
        end
    end  

    % Add metrics to a double matrix such that each column is a different
    % metric and each row is the corresponding particle in
    % nucleus_struct_protein
    nucleus_struct_metrics(:,1) = fluoMaxTime;
    nucleus_struct_metrics(:,2) = tOff;
    nucleus_struct_metrics(:,3) = fluoCenterMass;
    nucleus_struct_metrics(:,4) = fluoMean;
    nucleus_struct_metrics(:,5) = fluoVar;
%     nucleus_struct_metrics(:,6) = fluoStartEndRatio;

end