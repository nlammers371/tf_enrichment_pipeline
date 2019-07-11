clear
close all

% define core ID variables
project = 'Dl-Ven_snaBAC-mCh';
% project = 'Dl-Ven_hbP2P-mCh';
dropboxFolder =  'E:\Nick\LivemRNA\Dropbox\';
dataPath = [dropboxFolder 'ProcessedEnrichmentData\' project '/'];
% load training data set
load([dataPath 'event_training_data.mat'],'training_struct')
%%

% extract vectors
locus_protein_stack = training_struct.locus_protein_stack;
% locus_protein_stack = locus_protein_stack - nanmean(nanmean(locus_protein_stack,1),2);
% locus_protein_stack = locus_protein_stack ./ nanstd(nanstd(locus_protein_stack,[],1),[],2);
% locus_protein_stack(isnan(locus_protein_stack)) = 0;
% locus_protein_stack = mat2gray(locus_protein_stack);
prev_dur_vec = training_struct.prev_event_dur_vec;
class_vec = training_struct.event_class_vec;
% filter
n_half = 5000;
one_indices = find(class_vec == 1);
zero_indices = find(class_vec == 0);
train_indices = randsample([randsample(one_indices,n_half,true) randsample(zero_indices,n_half,true)],2*n_half,false);
test_indices = randsample([randsample(one_indices,n_half,true) randsample(zero_indices,n_half,true)],2*n_half,false);

n_row = size(locus_protein_stack,1);
n_col = size(locus_protein_stack,2);
n_slice = size(locus_protein_stack,3);
% initialize simple CNN
layers = [
    imageInputLayer([n_row n_col 1],'Name','input')
    
    convolution2dLayer(5,16,'Padding','same','Name','conv_1')
    batchNormalizationLayer('Name','BN_1')
    reluLayer('Name','relu_1')
    
%     maxPooling2dLayer(2,'Stride',2,'Name','pool_1')
    
    convolution2dLayer(3,32,'Padding','same','Name','conv_2')
    batchNormalizationLayer('Name','BN_2')
    reluLayer('Name','relu_2')
%     
%     maxPooling2dLayer(2,'Stride',2,'Name','pool_2')
%     
    convolution2dLayer(4,16,'Padding','same','Name','conv_3')
    batchNormalizationLayer('Name','BN_3')
    reluLayer('Name','relu_3')
    
    fullyConnectedLayer(2,'Name','fc')
    softmaxLayer('Name','softmax')
    classificationLayer('Name','classOutput')];

lgraph = layerGraph(layers);
% skipConv = convolution2dLayer(1,32,'Stride',2,'Name','skipConv');
% lgraph = addLayers(lgraph,skipConv);
% lgraph = connectLayers(lgraph,'relu_1','skipConv');
% lgraph = connectLayers(lgraph,'skipConv','add/in2');

figure
plot(lgraph)

% reshape data
LocusProteinTrain = locus_protein_stack(:,:,2,train_indices);
ClassVecTrain = categorical(double(class_vec(train_indices) == 1));

LocusProteinTest = locus_protein_stack(:,:,2,test_indices);
ClassVecTest = categorical(double(class_vec(test_indices) == 1));

% specify training options
options = trainingOptions('sgdm', ...
    'MaxEpochs',5, ...
    'Shuffle','every-epoch', ...
    'ValidationData',{LocusProteinTest,ClassVecTest}, ...
    'ValidationFrequency',30, ...
    'Verbose',false, ...
    'Plots','training-progress');

% train network
net = trainNetwork(LocusProteinTrain,ClassVecTrain,lgraph,options);

%% Examine network behavior
for i = 1:size(LocusProteinTest,4)
    imagesc(imgaussfilt(LocusProteinTest(:,:,1,i),1))
    pause(1)
end