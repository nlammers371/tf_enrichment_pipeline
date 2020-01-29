close all

% load('E:\Nick\LivemRNA\Dropbox (Personal)\LocalEnrichmentResults\2019-11-12-120mer-eGFP-lamin_10_250nmZStep_B\Spots_1300neighborhood.mat')

%order of entry - 1300, 1000
gauss3DIntensity(1,:) = horzcat(Spots.Fits.gauss3DIntensity);
gauss3DIntensityRaw(1,:) = horzcat(Spots.Fits.gauss3DIntensityRaw);
gauss3DIntensitySE(1,:) = horzcat(Spots.Fits.gauss3DIntensitySE);

gauss3DIntFilter = gauss3DIntensityRaw(1,:) .* (gauss3DIntensityRaw(1,:) < 300);
gauss3DIntFilter(gauss3DIntFilter==0) = NaN;
gauss3DIntMean = nanmean(gauss3DIntFilter);

nBins = 20;
edges = 0:25:250;

intensityFig = figure(1);
histogram(gauss3DIntensity(1,:),edges)

intRawFig = figure(2);
histogram(gauss3DIntensityRaw(1,:),edges)
title('3D Gaussian Fit Intensity - 120mer-eGFP')
xlabel('Fluorescence Intensity (au)')
ylabel('count')
ylim([0 200])
legend(['Mean = ' num2str(gauss3DIntMean)])

scatterFig = figure(3);
scatter(gauss3DIntensity(1,:),gauss3DIntensityRaw(1,:))
ylim([-50 400])
ylabel('gauss3DIntensity')
xlim([0 400])
xlabel('gauss3DIntensityRaw')


