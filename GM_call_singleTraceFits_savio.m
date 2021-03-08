% wrapper function to call the inference
clear 
close all

singleTraceFitsDir = pwd;
singleTraceFitsDir = [singleTraceFitsDir, filesep];
% move to inference directory
cd('~/repos/tf_enrichment_pipeline/');

% call the inference function
GM_main06_incorporate_hmm_results(singleTraceFitsDir);