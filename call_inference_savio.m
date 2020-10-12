% wrapper function to call the inference
clear 
close all

currentDir = pdw;
addpath(genpath(currentDir));

% call the inference function
main05_conduct_cpHMM_inference;