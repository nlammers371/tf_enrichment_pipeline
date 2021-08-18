clc
clear
close all

%% load data

load('S:\Jake\Dropbox\OptoDropbox\2021-07-12-optoknirps_vasa-eYFP_background_snapshot\CompiledNuclei.mat')
mean(MeanVectorAll)

%% calculate background

eYFP_mean = mean(MeanVectorAll)