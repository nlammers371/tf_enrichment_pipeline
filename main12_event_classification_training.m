clear
close all

% define core ID variables
project = 'Dl-Ven_snaBAC-mCh';
% project = 'Dl-Ven_hbP2P-mCh';
dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
dataPath = [dropboxFolder 'ProcessedEnrichmentData\' project '/'];
% load training data set
load([dataPath 'event_training_data.mat'],'training_struct')
%%
% initialize simple CNN
lgraph = layerGraph;

layers = [
    imageInputLayer([32 32 3],'Name','input')  
    convolution2dLayer(3,16,'Padding','same','Name','conv_1')
    batchNormalizationLayer('Name','BN_1')
    reluLayer('Name','relu_1')];
