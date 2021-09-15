% script to correct stripe position based on the curvature

clear % clear all variables in the workspace
close all % close any figure windows that might be open
clc

addpath(genpath('./lib'))


%% estimated parameter

export_rate = [4.66953396799217 5.13736103812975 5.50427054629742];
import_rate = [0.785442382061301 0.650352743543270 0.733385434567568];

export_half_time = log(2)./export_rate*60;
mean(export_half_time)
std(export_half_time)

import_half_time = log(2)./import_rate*60;
mean(import_half_time)
std(import_half_time)

%% make the plot

x = [export_half_time import_half_time];
g = {'export','export','export','import','import','import'};

fig = figure;
boxplot(x,g);
%boxplot(x,g,'PlotStyle','compact');
set(gca, 'YScale', 'log')

ylim([6 70])
ylabel('half-time (seconds)')
pbaspect([2 3 1])

%% try another version

x = [export_half_time; import_half_time]';
fig = figure;
boxchart(x)
set(gca, 'YScale', 'log')
ylim([6 100])
ylabel('half-time (seconds)')
pbaspect([2 3 1])
 