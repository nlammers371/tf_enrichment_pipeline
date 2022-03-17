% Script to call primary cpHMM wrapper function

clear
close all
warning('off','all') %Shut off Warnings
addpath(genpath('utilities'))

% set basic paths
% DataRoot = 'C:\Users\nlamm\Dropbox (Personal)\InductionLogic\';
DataRoot = 'S:\Nick\Dropbox\InductionLogic\';
customProjectPath = DataRoot;

% structure to store key inference info
inferenceInfo = struct;

% set projects to run
projectNameCell = {'20201117'};

% set inference options
inferenceInfo.ProteinBinFlag = 0;
inferenceInfo.FluoBinFlag = 0;
inferenceInfo.timeBins = {[0 Inf]};
% window_size = 15*60;
% infStartTimes = 5*30:300:85*30;
% inferenceInfo.timeBins = cell(size(infStartTimes));
% for i = 1:length(infStartTimes)
%     inferenceInfo.timeBins(i) = {[infStartTimes(i) infStartTimes(i)+window_size]};
% end
inferenceInfo.apBins = [0 Inf];
inferenceInfo.AdditionalGroupingVariable = '';

% set core model specs
inferenceInfo.modelSpecs.nStates = 3; % number of states in system
inferenceInfo.modelSpecs.nSteps = 5; % number of steps to traverse gene
inferenceInfo.modelSpecs.alphaFrac = 1.302/(1.302 +  1.085/2);

% other info
inferenceInfo.SampleSize = 5000;

% call main inference script
main05_conduct_cpHMM_inference(projectNameCell,'inferenceInfo',inferenceInfo,'customProjectPath',customProjectPath)