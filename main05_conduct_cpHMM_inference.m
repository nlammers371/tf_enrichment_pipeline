% Script to call primary cpHMM wrapper function

clear
close all
warning('off','all') %Shut off Warnings

% set project identifier
projectName = '2xDl-Ven_hbP2P-mCh';

% set core model specs
modelSpecs.nStates = 2; % number of states in system
modelSpecs.nSteps = 4; % number of steps to traverse gene
modelSpecs.alphaFrac = 1302/6444;

% Determine filepaths

% Call main inference function