clear 
close all

rawDynamicsData = 'E:\LocalEnrichment\Data\RawDynamicsData\';
date = '2019-06-10';
project = '60mer-eGFP_invitro_injectionSettings_375Hz70uW';
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
for i = 1:numel(stackList)
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
for im = 1:10%numel(imageStruct)
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
            gaussFitRaw = fit3DGaussianBeta(snipRaw);        
            gaussFitSeg = fit3DGaussianBeta(snipSeg);
            
            spotSnips(nSpots).GaussFitRaw = gaussFitRaw;
            spotSnips(nSpots).GaussFitSegmented = gaussFitSeg;
        end
    end
    
    if showFigures
        figure(5)
        hold off
    end
end
close(h);

maxZ = max(maxZStack);
maxX = max(maxXPix);
maxY = max(maxYPix);

%% 
gaussAmp = NaN(1,numel(spotSnips));
for i = 1:numel(spotSnips)
    gaussFitCurr = [spotSnips(i).GaussFitRaw];
%     if gaussFitCurr(1) >= 1 && spotSnips(i).Volume >= 5
        gaussAmp(i) = gaussFitCurr(1);
%     end
end
volume = [spotSnips.Volume];

meanAmp = nanmean(gaussAmp);
stdAmp = nanstd(gaussAmp);
varAmp = nanvar(gaussAmp);

spotSnipfile = [resultsPath project '.mat'];
save(spotSnipfile,'spotSnips')

histFig = figure(7);
histogram(gaussAmp)
xlabel('Amplitude of 3D Gaussian Fit (au)')
ylabel('Count')
title(project, 'Interpreter', 'none')
legend(['mean = ' num2str(meanAmp) '; std = ' num2str(stdAmp)])
% saveas(histFig,[figPath project '_hist.fig'])
% saveas(histFig,[figPath project '_hist.png'])

% Plotting to see if the really bright outliers are also the largest spots
% ampVolScatterFig = figure(8);
% scatter(gaussAmp, volume)
% xlabel('Amplitude of 3D Gaussian Fit (au)')
% ylabel('Volume of segmented spot (pixels)')
% saveas(ampVolScatterFig,[figPath project '_ampVolScatter.fig'])
% saveas(ampVolScatterFig,[figPath project '_ampVolScatter.png'])

