% Script to compare results across different target genes, TFs, and imaging
% conditions
clear 
close all
% set path to projects
DataPath = '../dat/';
% get list of projects
project_list = dir([DataPath '*_*']);
project_list = project_list([project_list.isdir]);
% specify indices of projects to include
include_vec = [1 5 6 7];
master_struct = struct;
for i = 1:numel(include_vec)
    load([DataPath project_list(include_vec(i)).name '/analysis_data.mat'])
    
    error('asfa')
end
