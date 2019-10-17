% Script to call mRNA profile plot function
clear 
close all

addpath('../../utilities')
project_cell = {'mHMMeve2_weka_inf_2018_05_07','eve7stripes_inf_2018_04_28'};
% specify data types to include
data_type_cell = {'reporter','BAC'};
nc14_qc_flag_vec = [1 1];
ymax_vec = [.11 .14];
% basic plot info 

% %%%%% Main Text %%%%%%%
%% % mRNA profiles
plot_time = 2410; % Reference time for plotting  
half_life_vec = [1 7 15 Inf]*60; % seconds
for i = 1:numel(project_cell)
    for j = 1:numel(half_life_vec)
        make_mRNA_profile_plots(project_cell{i},data_type_cell{i},half_life_vec(j),plot_time, 'ymax', ymax_vec(i))
    end
end

%% % ap trend plots
make_ap_trend_plots(project_cell{1})
%% 
% %%%%% SI %%%%%%%
% %%
 % SI 1: time-dependent mean rate plots
% make_temporal_mean_rate_plots(project_cell{1})
% 
%% % SI 2: On and off times by AP
make_ap_on_off_plots(project_cell{1})

%% SI 3:
make_jack_mRNA_plot(project_cell{1},data_type_cell{1},7*60,plot_time)
%% ap patch movies
PreProcPath = '../../../../figure_data/PreProcessedData/';
PostProcPath = '../../../../figure_data/DynamicsResults/';
Prefix = '2017-06-15-eve2_20sec_10';
nc = 14;
close all
% call function for duration
make_patch_movies(Prefix, PreProcPath, PostProcPath, nc, 'duration','last_time',40)

% call function for fraction on
make_patch_movies(Prefix, PreProcPath, PostProcPath, nc, 'on_only','last_time',40)

% call function for mean rate
make_patch_movies(Prefix, PreProcPath, PostProcPath, nc)