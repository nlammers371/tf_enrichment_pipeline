clear 
close all

rawDynamicsData = 'E:\LocalEnrichment\Data\RawDynamicsData\';
% date = '2019-06-10';
date = '2019-08-19';
% project = '60mer-eGFP_invitro_invivoSettings_025umzStep7uW_01';
project = '60mer-Ven_1-100000pos_7uW';
rawPath = [rawDynamicsData filesep date filesep project filesep];
segPath = [rawDynamicsData filesep date filesep project filesep 'SegmentedImages' filesep];
resultsPath = ['E:\Nick\LivemRNA\Dropbox\LocalEnrichmentResults\' date '_' project filesep];
figPath = 'E:\Nick\LivemRNA\Dropbox\LocalEnrichmentFigures\calibration\';

mkdir(resultsPath)

tileScan = true;  %toggles between Time Series and Tile Scan settings
showFigures = false;
minRegionSize = 3;  %in 3D, mostly trying to get rid of hot pixels
maxRegionSize = 70; %Arbitrary to avoid huge aggregates

% Positioning limits - set based on max sizes in all 3D found below
minXpos = 4;
minYpos = 5;
minZpos = 4;
% Snip size
xySnipRadius = 7;
zSnipRadius = 4;

fileList = dir([rawPath '*.tif']);

%Get list of unique z-stacks
stackNamesCell = cell(1,numel(fileList));
if tileScan
    % unique stacks from Tile Scan labeled as '...s##...'
    searchTerm = 's\d';
else
    % unique stacks from Time Series labeled as '...t##...'
    searchTerm = 't\d';
end
for i = 1:numel(fileList)
    name = fileList(i).name;
    stopIndex = regexp(name,searchTerm) + 2;
    stackNamesCell{i} = name(1:stopIndex);
end
stackList = unique(stackNamesCell);

% Compile unique z stacks into structure
imageStruct = struct;
for i = 1:numel(stackList)%:9
    currStackName = stackList{i};
    indices = find(contains(stackNamesCell,currStackName));
    imStack = [];
    for j = indices
        im = imread([rawPath fileList(j).name]);
        imStack = cat(3,imStack,im);
    end
    imageStruct(i).imStack = imStack;
    imageStruct(i).source = stackList{i};
end

numZSlices = size(imageStruct(1).imStack, 3);
maxZStack = NaN(1,numel(imageStruct));
maxXPix = NaN(1,numel(imageStruct));
maxYPix = NaN(1,numel(imageStruct));
nSpots = 0;
spotSnips = struct;
% imshow(imageStruct(1).imStack(:,:,11),[])

% Segment spots & fit 3D gaussians to them
h = waitbar(0, 'Finding & fitting spots');
tic
for im = 1:numel(imageStruct)
    waitbar(im/numel(imageStruct), h, ['Finding & fitting spots in frame ' num2str(im) ' of ' num2str(numel(imageStruct))])
    imStackRaw = imageStruct(im).imStack;
    % Remove hot pixels and large aggregates
    imStackFilter = bwareaopen(imStackRaw > 0,minRegionSize,18);  %limiting connectivity so a series of diagonally connected hot pixels isn't considered an object
    imStackClean = immultiply(imStackRaw,imStackFilter);
%     % Estimate background for correction
%     backgroundFluo = nanmean(imStackRaw(~imStackFilter));
%     imshowpair(imStackRaw(:,:,14),imStackClean(:,:,14))
%     pause(1)
    
    % get image dimensions
    xDim = size(imStackFilter,2);
    yDim = size(imStackFilter,1);
    zDim = size(imStackFilter,3);  
    
    % Perform difference of gaussians filtering  
    imGauss1 = imgaussfilt3(imStackClean,1);
    imGauss5 = imgaussfilt3(imStackClean,5);
    imGaussDiff = imGauss1 - imGauss5;
    imProject = max(imGaussDiff,[],3);
    if showFigures
        imagesc(imProject)
        pause(1)
    end
    % Threshold the stack using otsu's method
    thresh = multithresh(imGaussDiff,1);
    imStackFilt = imGaussDiff > thresh;
%     imshowpair(max(imStackRaw,[],3),max(imStackFilt,[],3),'Scaling','independent')
%     pause(1)
    imStackNew = immultiply(imStackRaw, imStackFilt);
%     imshow(max(imStackNew,[],3),[])
%     pause(1)
    imStackLabel = bwconncomp(imStackFilt);
    imProps = regionprops3(imStackLabel,'Centroid','Solidity','Volume','VoxelList');
    
    % Find the x, y, and z sizes of all the spots
    nZSliceVec = NaN(1,numel(imProps.VoxelList));
    nYPixVec = NaN(1,numel(imProps.VoxelList));
    nXPixVec = NaN(1,numel(imProps.VoxelList));
    for j = 1:numel(imProps.VoxelList)
        voxelVec = imProps.VoxelList{j};
        nZSliceVec(j) = numel(unique(voxelVec(:,3)));
        nYPixVec(j) = numel(unique(voxelVec(:,1)));
        nXPixVec(j) = numel(unique(voxelVec(:,2)));
    end
    % Find the max number of pixels occupied by the spots in all dimensions
    maxZStack(im) = max(nZSliceVec);
    maxXPix(im) = max(nXPixVec);
    maxYPix(im) = max(nYPixVec);
   
    % Quality control, removing spots too close to the edge of the image 
    % and spots that are too small
    zPosVec = imProps.Centroid(:,3);
    yPosVec = imProps.Centroid(:,1);
    xPosVec = imProps.Centroid(:,2);
    volVec = [imProps.Volume];
    solVec = [imProps.Solidity];
    keepIndices = find(volVec >= minRegionSize ...%& volVec <= maxRegionSize ...
                        & (zPosVec >= minZpos) & (zPosVec <= (zDim - minZpos + 1)) ...
                        & (yPosVec >= minYpos) & (yPosVec <= (yDim - minYpos + 1)) ...
                        & (xPosVec >= minXpos) & (xPosVec <= (xDim - minXpos + 1)) ... %& (solVec >= .8) ...
                        & (nZSliceVec' > 1));               
    imPropsClean = imProps(keepIndices,:);
    % plot inferred spot centroids over top of original image
    centroidArray = vertcat(imPropsClean.Centroid);
    if showFigures
        figure(4);
        imagesc(max(imStackNew,[],3));
        hold on
        scatter(centroidArray(:,1),centroidArray(:,2),20,'green')
        hold off
    end
    
    if showFigures
        figure(5);
        imagesc(max(imStackNew,[],3));
        hold on
    end
    % Fit a 3D gaussian to snips containing each spot
    for k = 1:size(centroidArray,1)
        yCentroid = round(centroidArray(k,1));
        xCentroid = round(centroidArray(k,2));
        zCentroid = round(centroidArray(k,3));
        % Define a X x Y x Z pixel image around the centroid
        yRange = (yCentroid - xySnipRadius):(yCentroid + xySnipRadius); 
        xRange = (xCentroid - xySnipRadius):(xCentroid + xySnipRadius);
        zRange = (zCentroid - zSnipRadius):(zCentroid + zSnipRadius); 
        % Account for spots on the edge by trimming out of bounds indices
        yRangeTrim = yRange(find(yRange > 0 & yRange <= yDim));
        xRangeTrim = xRange(find(xRange > 0 & xRange <= xDim));
        zRangeTrim = zRange(find(zRange > 0 & zRange <= zDim));
        
        snipRaw = imStackRaw(xRangeTrim,yRangeTrim,zRangeTrim);    % Don't know why I needed to flip x & y here
        snipSeg = imStackNew(xRangeTrim,yRangeTrim,zRangeTrim);    % Don't know why I needed to flip x & y here
        snipRawProj = max(snipRaw,[],3);
        snipSegProj = max(snipSeg,[],3);
        snipFilt = snipSeg > 0;
        
        % Only save snips that contain a single spot (one labeled region)
        if max(max(max(bwlabeln(snipFilt)))) == 1
            nSpots = nSpots + 1;
            spotSnips(nSpots).SnipRaw = snipRaw;
            spotSnips(nSpots).SnipSegmented = snipSeg;
            spotSnips(nSpots).Volume = imProps.Volume(k);
            
            % Plot snips to check
            if showFigures
                figure(5)
                scatter(yCentroid,xCentroid,20,'green')

                figure(6)
                imagesc(snipRawProj)
                spotTitle = ['spot number: ' num2str(nSpots)];
                title(spotTitle)
                drawnow
                pause(0.5)
            end
            
            % Fit 3D gaussians
            [gaussFitRaw, gaussTotalFluoRaw] = fit3DGaussianBeta(snipRaw);        
            [gaussFitSeg, gaussTotalFluoSeg] = fit3DGaussianBeta(snipSeg);
            
            spotSnips(nSpots).GaussFitRaw = gaussFitRaw;
            spotSnips(nSpots).GaussTotalFluoRaw = gaussTotalFluoRaw;
            spotSnips(nSpots).GaussFitSegmented = gaussFitSeg;
            spotSnips(nSpots).GaussTotalFluoSegmented = gaussTotalFluoSeg;
        end
    end
    
    if showFigures
        figure(5)
        hold off
    end
end
toc
close(h);

maxZ = max(maxZStack);
maxX = max(maxXPix);
maxY = max(maxYPix);

%%
gaussAmp = NaN(1,numel(spotSnips)); % amplitudes of the gaussian fits
gaussCenter = NaN(3,numel(spotSnips));  % y,x,and z positions
gaussSigma = NaN(3,numel(spotSnips));   % y,x,and z sigma values
gaussSigCov = NaN(3,numel(spotSnips)); % y,x,and z sigma covariance coefficients
gaussOffset = NaN(1,numel(spotSnips)); % inferred background offset
gaussTotalFluo = NaN(1,numel(spotSnips));% Estimate of total fluorescence using gaussian fit parameters

for i = 1:numel(spotSnips)
    [gaussFitCurr] = [spotSnips(i).GaussFitRaw];
    gaussAmp(1,i) = gaussFitCurr(1);
    gaussCenter(:,i) = gaussFitCurr(2:4);
    gaussSigma(:,i) = gaussFitCurr(5:7);
    gaussSigCov(:,i) = gaussFitCurr(8:10);
    gaussOffset(1,i) = gaussFitCurr(11);
    gaussTotalFluo(1,i) = spotSnips(i).GaussTotalFluoRaw;
end
volume = [spotSnips.Volume];

% meanAmp = nanmean(gaussAmp);
% stdAmp = nanstd(gaussAmp);
% varAmp = nanvar(gaussAmp);

meanFluo = nanmean(gaussTotalFluo);
stdFluo = nanstd(gaussTotalFluo);
varFluo = nanvar(gaussTotalFluo);

spotSnipfile = [resultsPath project '.mat'];
save(spotSnipfile,'spotSnips')

histFig = figure(7);
histogram(gaussTotalFluo)
xlabel('Total Fluorescence from 3D Gaussian Fit (au)')
ylabel('Count')
title([project, ' - Total Fluorescence of Gaussian Spot Fit'], 'Interpreter', 'none')
legend(['mean = ' num2str(meanFluo) '; std = ' num2str(stdFluo)])
% saveas(histFig,[figPath project '_hist.fig'])
saveas(histFig,[figPath project '_hist.png'])

% Plotting to see if the really bright outliers are also the largest spots
% ampVolScatterFig = figure(8);
% scatter(gaussAmp, volume)
% xlabel('Amplitude of 3D Gaussian Fit (au)')
% ylabel('Volume of segmented spot (pixels)')
% saveas(ampVolScatterFig,[figPath project '_ampVolScatter.fig'])
% saveas(ampVolScatterFig,[figPath project '_ampVolScatter.png'])

%% Use 0.2 x 0.2 x 0.5 um disk to calculate calibration
yxRadius = 0.2; %um
pxSize = 0.107; %um
zWidth = 0.5; %um
zStep = 0.5; %um
pxRadius = round(yxRadius / pxSize); % Calculating how many pixels are in the radius

diskTotalFluo_vec = NaN(1,numel(spotSnips));
diskMeanFluo_vec = NaN(1,numel(spotSnips));

for i = 1:numel(spotSnips)
    snipRaw = spotSnips(i).SnipRaw;
    snipSeg = spotSnips(i).SnipSegmented;

    %Find brightess z slice
    meanZStep = mean(snipRaw, [1 2]);
    [~, maxMeanZ] = max(meanZStep);
%     sumZStep(1,1:9) = sum(SnipRaw, [1 2]);
%     [~, maxSumZ] = max(sumZStep);
    brightZSnip = double(snipRaw(:,:,maxMeanZ));
%     imshow(brightZSnip,[])
    
    %Make an image where the center pixel from the gaussian fit is 1
    gaussCenter(1,1:3) = spotSnips(i).GaussFitRaw(2:4); % y, x, z position of center
    gaussCenter = round(gaussCenter); % round to nearest pixel
    centerIm = zeros(size(brightZSnip));
    centerIm(gaussCenter(1,1),gaussCenter(1,2)) = 1;
    
    %Create a mask with the distance from the "center" pixel
    distIm = bwdist(centerIm);
    % Mask to keep only the pixels within the defined disk
    brightZDisk = brightZSnip(distIm < pxRadius);
    
    %Calculate total fluorescence in disk
    diskTotalFluo = sum(brightZDisk);
    diskMeanFluo = mean(brightZDisk);
    diskTotalFluo_vec(i) = diskTotalFluo;
    diskMeanFluo_vec(i) = diskMeanFluo;
end

diskMeanFig = figure(8);
histogram(diskMeanFluo_vec)
xlabel('Mean Fluorescence (au)')
ylabel('Count')
title([project, ' - Mean Fluorescence of 0.2 x 0.2 x 0.5 um disk'], 'Interpreter', 'none')
legend(['mean = ' num2str(nanmean(diskMeanFluo_vec)) '; std = ' num2str(nanstd(diskMeanFluo_vec))])
% saveas(histFig,[figPath project '_hist.fig'])
saveas(diskMeanFig,[figPath project '_diskMean_hist.png'])

diskTotalFig = figure(9);
histogram(diskTotalFluo_vec)
xlabel('Total Fluorescence (au)')
ylabel('Count')
title([project, ' - Total Fluorescence of 0.2 x 0.2 x 0.5 um disk'], 'Interpreter', 'none')
legend(['mean = ' num2str(nanmean(diskTotalFluo_vec)) '; std = ' num2str(nanstd(diskTotalFluo_vec))])
% saveas(histFig,[figPath project '_hist.fig'])
saveas(diskTotalFig,[figPath project '_diskTotal_hist.png'])
