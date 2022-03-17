
clear all
close all
clc


%% initialization

% Define some other colors  
yw = [234 194 100]/255; % yellow
bl = [115 143 193]/255; % blue
gr = [191 213 151]/255; % green

%% Figure: burst frequency and duration

% load inference results
%HIGH_bursting = load('P:\Jake\Pipeline\tf_enrichment_pipeline\opto_knirps_project\bursting_data\compiledResults_w7_K2_p0_ap1_t1_f2D_0.03_15_30_HIGH.mat');
%LOW_bursting = load('P:\Jake\Pipeline\tf_enrichment_pipeline\opto_knirps_project\bursting_data\compiledResults_w7_K2_p0_ap1_t1_f2D_0.03_15_30_LOW.mat');

HIGH_bursting = load('P:\Jake\Pipeline\tf_enrichment_pipeline\opto_knirps_project\bursting_data\compiledResults_w7_K2_p0_ap1_t1_f2D_0.03_15_30_HIGH_test.mat');
LOW_bursting = load('P:\Jake\Pipeline\tf_enrichment_pipeline\opto_knirps_project\bursting_data\compiledResults_w7_K2_p0_ap1_t1_f2D_0.03_15_30_LOW_test.mat');



% data from HIGH
burst_freq_HIGH = HIGH_bursting.compiledResults.freq_vec_mean;
burst_freq_ste_HIGH = HIGH_bursting.compiledResults.freq_vec_ste;
burst_dur_HIGH = HIGH_bursting.compiledResults.dur_vec_mean;
burst_dur_ste_HIGH = HIGH_bursting.compiledResults.dur_vec_ste;
burst_rate_HIGH= HIGH_bursting.compiledResults.init_vec_mean*1e-5;
burst_rate_ste_HIGH = HIGH_bursting.compiledResults.init_vec_ste*1e-5;

% data from LOW
burst_freq_LOW = LOW_bursting.compiledResults.freq_vec_mean;
burst_freq_ste_LOW = LOW_bursting.compiledResults.freq_vec_ste;
burst_dur_LOW = LOW_bursting.compiledResults.dur_vec_mean;
burst_dur_ste_LOW = LOW_bursting.compiledResults.dur_vec_ste;
burst_rate_LOW= LOW_bursting.compiledResults.init_vec_mean*1e-5;
burst_rate_ste_LOW = LOW_bursting.compiledResults.init_vec_ste*1e-5;

%% plot the result
x = [1 2];
data_freq = [burst_freq_HIGH burst_freq_LOW];
err_freq = [burst_freq_ste_HIGH burst_freq_ste_LOW];

data_dur = [burst_dur_HIGH burst_dur_LOW];
err_dur = [burst_dur_ste_HIGH burst_dur_ste_LOW];

data_rate = [burst_rate_HIGH burst_rate_LOW];
err_rate = [burst_rate_ste_HIGH burst_rate_ste_LOW];


fig = figure;
tiledlayout(1,3)

nexttile
bar(x,data_freq)              
hold on

er = errorbar(x,data_freq,err_freq);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

nexttile

bar(x,data_dur)              
hold on

er = errorbar(x,data_dur,err_dur);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 

nexttile

bar(x,data_rate)              
hold on

er = errorbar(x,data_rate,err_rate);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 

%pbaspect([1 1 1])



