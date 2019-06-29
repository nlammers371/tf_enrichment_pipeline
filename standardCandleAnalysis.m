clear 
close all

rawDynamicsData = 'E:\LocalEnrichment\Data\RawDynamicsData\';
date = '2019-06-10\';
project = '60mer-eGFP_invitro_invivoSettings_025umzStep7uW_01';
rawPath = [rawDynamicsData filesep date filesep project filesep];
segPath = [rawDynamicsData filesep date filesep project filesep 'SegmentedImages' filesep];

tileScan = true;  %toggles between Time Series and Tile Scan settings
showFigures = false;
minRegionSize = 2;  %in 3D, mostly trying to get rid of hot pixels
maxRegionSize = 70; %Arbitrary to avoid huge aggregates

% Positioning limits - set based on max sizes in all 3D found below
minXpos = 4;
minYpos = 5;
minZpos = 4;

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
for im = 1:numel(imageStruct)
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
    
        % perform difference of gaussians filtering  
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

    imStackLabel = bwconncomp(imStackFilt);
    imProps = regionprops3(imStackLabel,'Centroid','Solidity','Volume','VoxelList');
    
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
   
    % Quality control, just position and minimum size for now
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
    showFigures = true;
    if showFigures
        figure(4);
        imagesc(imProject);
        hold on
        scatter(centroidArray(:,1),centroidArray(:,2),20,'green')
    end
    
    for k = size(centroidArray,1)
        nSpots = nSpots + 1;
        centroid = round(centroidArray(k,:));
        xRange = centroid(2)-7:centroid(2)+7;
        yRange = centroid(1)-7:centroid(1)+7;
        zRange = centroid(3)-4:centroid(3)+4;
        spotSnips(nSpots).Snip = imStackRaw(yRange,xRange,zRange);
        
    end
end


maxZ = max(maxZStack);
maxX = max(maxXPix);
maxY = max(maxYPix);

