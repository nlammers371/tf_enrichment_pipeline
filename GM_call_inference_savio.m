% wrapper function to call the inference
clear 
close all

inferenceDir = pwd;
inferenceDir = [inferenceDir, filesep];
% move to inference directory
cd('~/repos/tf_enrichment_pipeline/');

% call the inference function
GM_main05_conduct_cpHMM_inference(inferenceDir);