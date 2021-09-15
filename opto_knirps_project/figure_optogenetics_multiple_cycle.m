% script to correct stripe position based on the curvature

clear % clear all variables in the workspace
close all % close any figure windows that might be open
clc

addpath(genpath('./lib'))
%% Part 0: Initialization

% We need to specify the prefix of the project we wish to analyze

%optokni_eve4+6_WT
Prefix = '2021-07-05-optoknirps_eve4_6_multiple_export_embryo1';

analysisInit = 0;
analysisFinal = 35;

ap_lim = 0.02; % AP range for analysis, 0.02 seems to be a reasonable number

eYFP_background = 375698.13;


% color to be used
k_green = brighten([38 142 75]/256,.4);
color_green = [38 143 75]/256; % color from Jake
mRNA_red = brighten([212 100 39]/256,.2);

%% Part 1: Add working path and load the data
% We will initialize the paths that are used later

% Specify the position of DynamicsResults folder
%DynamicsResultsFolder = '/Users/jiaxi.zhao/Dropbox/OptoDropbox';
DynamicsResultsFolder = 'S:/Jake/Dropbox/OptoDropbox';
%% 
% Now, load all the data

% Now we'll load the CompiledParticles dataset
% This generates a string that gives the path to CompiledParticles
LoadNucleiPath = [DynamicsResultsFolder '/' Prefix '/' Prefix '_lin.mat'];
LoadComNucleiPath = [DynamicsResultsFolder '/' Prefix '/CompiledNuclei.mat']; 

data.nuclei = load(LoadComNucleiPath);
data.sch = load(LoadNucleiPath);

ElapsedTime = data.nuclei.ElapsedTime-data.nuclei.ElapsedTime(data.nuclei.nc14);

% find frames to analyze
initFrame = find(ElapsedTime<=analysisInit,1,'last');
finalFrame = find(ElapsedTime>analysisFinal,1);

% total number of frames
totalFrame = finalFrame-initFrame + 1;

totalTime = ElapsedTime(initFrame:finalFrame);
%% Part 2: Quality check on nuclei tracking

% Quality check and save the nuclei that passed the test
% Requires the nuclei to have continuous trace until lateral movement
X_pass = [];
Y_pass = [];

X_fail = [];
Y_fail = [];

schPass = [];

schnitzcells = data.sch.schnitzcells;

for i = 1:size(schnitzcells,2)
    index_s = find(schnitzcells(i).frames == initFrame);
    index_f = find(schnitzcells(i).frames == finalFrame);
    fluo_temp = max(schnitzcells(i).Fluo(index_s:index_f,:),[],2);
    result = sum(isnan(fluo_temp));
    % quality check
    if ~isempty(index_s) && ~isempty(index_f) && (index_f-index_s == finalFrame-initFrame)  ...
            && (result == 0)
        schPass = [schPass, i];
        X_pass = [X_pass, schnitzcells(i).cenx(index_s)];
        Y_pass = [Y_pass, schnitzcells(i).ceny(index_s)];
    else
        X_fail = [X_fail, schnitzcells(i).cenx(index_s)];
        Y_fail = [Y_fail, schnitzcells(i).ceny(index_s)];
    end
end

fig = figure(1);
plot(X_pass,Y_pass,'o')
hold on
plot(X_fail,Y_fail,'o')
hold off
axis equal
xlim([0 768])
ylim([0 450])
title('Nuclei quality check')

%% Part 3: Compile the schnitzcells together
% compile based on each frame

nucleiNum = length(schPass);

% Initialize storage
processed_data(1).xcoord = [];
processed_data(1).ycoord = [];

processed_data(1).schnitznum = [];
processed_data(1).NuclearFluo = [];


% Compile all the nuclei in each frame and assign basic info
for i = 1:nucleiNum
    schTemp = schPass(i);
    index_s = find(schnitzcells(schTemp).frames == initFrame);
    index_f = find(schnitzcells(schTemp).frames == finalFrame);
    
    for j = index_s:index_f
        frameTemp = schnitzcells(schTemp).frames(j);
        x_coord = schnitzcells(schTemp).cenx(j);
        y_coord = schnitzcells(schTemp).ceny(j);
        fluo = max(schnitzcells(schTemp).Fluo(j,:));
        try
            processed_data(frameTemp).xcoord = [processed_data(frameTemp).xcoord, x_coord];
            processed_data(frameTemp).ycoord = [processed_data(frameTemp).ycoord, y_coord];
            processed_data(frameTemp).schnitznum = [processed_data(frameTemp).schnitznum, schTemp];
            processed_data(frameTemp).NuclearFluo = [processed_data(frameTemp).NuclearFluo fluo];  
        catch
            processed_data(frameTemp).xcoord = x_coord;
            processed_data(frameTemp).ycoord = y_coord;
            processed_data(frameTemp).schnitznum = schTemp;
            processed_data(frameTemp).NuclearFluo = fluo;
        end
    end
end

% Compile based on each nuclei

for i = 1:nucleiNum
    for j = initFrame:finalFrame
        
        % Transform the actual frame count to analyzed frame count
        frameTemp = j-initFrame+1;
        
        CompiledNucleiData(i).xcoord(frameTemp) = processed_data(j).xcoord(i);
        CompiledNucleiData(i).ycoord(frameTemp) = processed_data(j).ycoord(i);
        
        CompiledNucleiData(i).schnitzNum(frameTemp) = processed_data(j).schnitznum(i);
        CompiledNucleiData(i).nuclearFluo(frameTemp) = processed_data(j).NuclearFluo(i);
        
    end
end
% Calculate average coordinate for each nuclei

for i = 1:nucleiNum
    CompiledNucleiData(i).xcoordMean = mean(CompiledNucleiData(i).xcoord);
    CompiledNucleiData(i).ycoordMean = mean(CompiledNucleiData(i).ycoord);
end


%% Correct AP position

ap_map = load(['ap_map' filesep Prefix '_APmap.mat']);

for i = 1:nucleiNum
    CompiledNucleiData(i).APPos = ap_map.APmap(ceil(CompiledNucleiData(i).ycoordMean),ceil(CompiledNucleiData(i).xcoordMean));
end

%% 

nucleiList = [];

for i = 1:nucleiNum
    
    if (CompiledNucleiData(i).APPos>-ap_lim) & (CompiledNucleiData(i).APPos<ap_lim)
       nucleiList = [nucleiList i]; 
    end
    
end

nuclearFluo_mean = zeros(totalFrame,1);

for i = 1:totalFrame
    
    fluo_temp = [];
    
    for j = 1:length(CompiledNucleiData)
        
        fluo_temp = [fluo_temp CompiledNucleiData(j).nuclearFluo(i)];
        
    end
    
    nuclearFluo_mean(i) = mean(fluo_temp);
    
end

%%

time_vec = ElapsedTime(initFrame:finalFrame);

on_frame = 50;
off_frame = 114;

on_frame_2 = 186;
off_frame_2 = 246;

nuclearFluo_mean(on_frame:off_frame) = convert_from_458(nuclearFluo_mean(on_frame:off_frame));
nuclearFluo_mean(on_frame_2:off_frame_2) = convert_from_458(nuclearFluo_mean(on_frame_2:off_frame_2));
nuclearFluo_mean = nuclearFluo_mean-eYFP_background;

fig = figure;
plot(time_vec,nuclearFluo_mean,'-o')

hold on
%errorbar(time_vec,knirps_vec_ste,'Color','k','CapSize',0);
plot(time_vec,nuclearFluo_mean,'-k','LineWidth',1)
scatter(time_vec,nuclearFluo_mean,50,'MarkerFaceColor',k_green,'MarkerEdgeColor','k')

ylim([1.5E5 6E5])

