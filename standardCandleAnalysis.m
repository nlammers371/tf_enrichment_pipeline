clear 
close all

rawDynamicsData = 'E:\LocalEnrichment\Data\RawDynamicsData\';
date = '2019-06-10\';
project = '60mer-eGFP_invitro_invivoSettings_025umzStep7uW_01';
rawPath = [rawDynamicsData filesep date filesep project filesep];
segPath = [rawDynamicsData filesep date filesep project filesep 'SegmentedImages' filesep];

tileScanOn = true;  %toggles between Time Series and Tile Scan settings

fileList = dir([rawPath '*.tif']);

%Get list of unique z-stacks
stackNamesCell = cell(1,numel(fileList));
if tileScanOn
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
        im = imread([dataPath fileList(j).name]);
        imStack = cat(3,imStack,im);
    end
    imageStruct(i).imStack = imStack;
    imageStruct(i).source = stackList{i};
end