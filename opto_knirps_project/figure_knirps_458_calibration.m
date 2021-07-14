% script to correct stripe position based on the curvature

clear % clear all variables in the workspace
close all % close any figure windows that might be open
clc

addpath(genpath('./lib'))
%% Part 0: Initialization
% Specify the dataset to analyze

% We need to specify the prefix of the project we wish to analyze
Prefix = {'2021-06-01-knirps_458_calibration_embryo1','2021-06-01-knirps_458_calibration_embryo2'};

%% Part 1: Add working path and load the data
% We will initialize the paths that are used later

% Specify the position of DynamicsResults folder
DynamicsResultsFolder = 'S:/Jake/Dropbox/OptoDropbox';
%% 

% initialize parameter
fluo_wt = [];
fluo_458 = [];

for k = 1:length(Prefix)
    % Now, load the data
    LoadComNucleiPath = [DynamicsResultsFolder '/' Prefix{k} '/CompiledNuclei.mat']; 
    data.nuclei = load(LoadComNucleiPath);

    % total number of frames
    totalFrame = size(data.nuclei.ElapsedTime,2);

    % Compile information

    CompiledNuclei = data.nuclei.CompiledNuclei;


    for i = 1:length(CompiledNuclei)
        frame_temp = CompiledNuclei(i).Frames;
        fluo_temp = CompiledNuclei(i).FluoMax;

        for j = 1:length(frame_temp)-1
            if (mod(frame_temp(j),2) == 1) && ~isnan(fluo_temp(j)) && ~isnan(fluo_temp(j+1))
                fluo_wt = [fluo_wt fluo_temp(j)];
                fluo_458 = [fluo_458 fluo_temp(j+1)];
            end
        end

    end
end

%% plot the result

ft = fittype('a*x');

% perform fitting
f1 = fit(fluo_wt',fluo_458',ft);
f2 = fit(fluo_wt',fluo_458','poly1');
k = f1.a
x_plot = linspace(0,1.2*max(fluo_wt),1000);


fig = figure;
plot(fluo_wt,fluo_458,'.')
hold on
plot(x_plot,f1(x_plot),'LineWidth',1.5);
xlim([0 2.25E6])
ylim([0 3E6])
xlabel('nuclear fluorescence (AU) with 458nm laser OFF')
ylabel('nuclear fluorescence (AU) with 458nm laser ON')

pbaspect([3 2 1])

%% linear regression result
t = fitlm(fluo_wt,fluo_458)






