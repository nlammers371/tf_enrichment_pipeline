% script to correct stripe position based on the curvature

clear % clear all variables in the workspace
close all % close any figure windows that might be open
clc

addpath(genpath('./lib'))
%% Part 0: Initialization

% We need to specify the prefix of the project we wish to analyze


%%
%Prefix = '2021-07-06-optoknirps_eve4_6_parameter_embryo1';
%on_frame = 60;
%on_fit_end = 78;

%off_frame = 100;
%off_fit_end = 132;

%%
Prefix = '2021-07-07-optoknirps_eve4_6_parameter_embryo2';
on_frame = 43;
on_fit_end = 62; % 62, after 2min; 71, after 3min

off_frame = 84;
off_fit_end = 119; % after 4min

%%
%Prefix = '2021-07-07-optoknirps_eve4_6_parameter_embryo3';
%on_frame = 40;
%on_fit_end = 59;

%off_frame = 81;
%off_fit_end = 113;

%%

analysisInit = 0;
analysisFinal = 17.5;

ap_lim = 0.02; % AP range for analysis, 0.02 seems to be a reasonable number

eYFP_background = 375698.13;

fontSize = 10;

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

%% Calculate mean

nucleiList = [];

for i = 1:nucleiNum
    
    if (CompiledNucleiData(i).APPos>-ap_lim) & (CompiledNucleiData(i).APPos<ap_lim)
       nucleiList = [nucleiList i]; 
    end
    
end

nuclearFluo_mean = zeros(totalFrame,1);
nuclearFluo_std = zeros(totalFrame,1);

for i = 1:totalFrame
    
    fluo_temp = [];
    
    for j = 1:length(nucleiList)
        
        fluo_temp = [fluo_temp CompiledNucleiData(nucleiList(j)).nuclearFluo(i)];
        
    end
    
    nuclearFluo_mean(i) = mean(fluo_temp);
    nuclearFluo_std(i) = std(fluo_temp);
    
end

%% Plot result

time_vec = ElapsedTime(initFrame:finalFrame)-ElapsedTime(on_frame);

nuclearFluo_mean(on_frame:off_frame) = convert_from_458(nuclearFluo_mean(on_frame:off_frame));
nuclearFluo_mean = nuclearFluo_mean-eYFP_background;

fig = figure;

hold on
%errorbar(time_vec,nuclearFluo_mean,nuclearFluo_std,'Color','k','CapSize',0);
plot(time_vec,nuclearFluo_mean,'-k','LineWidth',1)
scatter(time_vec,nuclearFluo_mean,50,'MarkerFaceColor',k_green,'MarkerEdgeColor','k')
xlim([-2 10])
ylim([3.5E5 12E5])

pbaspect([3 2 1])

%% try to interpolate the data

time_interp = -2:1/10:12;

nuclearFluo_interp = interp1(time_vec,nuclearFluo_mean,time_interp);

nuclearFluo_interp(51:61) = movmean(nuclearFluo_interp(51:61),10);

fig = figure;

hold on
plot(time_interp,nuclearFluo_interp,'-k','LineWidth',1)
scatter(time_interp,nuclearFluo_interp,50,'MarkerFaceColor',k_green,'MarkerEdgeColor','k')
xlim([-2 10])
ylim([3.5E5 11E5])

pbaspect([3 2 1])


%% Fit export half-time

X = time_vec(on_frame:on_fit_end)-time_vec(on_frame);
Y = nuclearFluo_mean(on_frame:on_fit_end)';

% Plot the initial data.
fig = figure;
plot(X', Y', 'b*', 'LineWidth', 2, 'MarkerSize', 15);
grid on;

% Convert X and Y into a table, which is the form fitnlm() likes the input data to be in.
tbl = table(X', Y');
% Note how this "x" of modelfun is related to big X and big Y.
% x((:, 1) is actually X and x(:, 2) is actually Y - the first and second columns of the table.
modelfun = @(b,x) b(1) * exp(-b(2)*x(:, 1)) + b(3);  
beta0 = [2E5 0 2E5]; % Guess values to start with.  Just make your best guess.
% Now the next line is where the actual model computation is done.
mdl_export = fitnlm(tbl, modelfun, beta0);
% Now the model creation is done and the coefficients have been determined.

% Extract the coefficient values from the the model object.
% The actual coefficients are in the "Estimate" column of the "Coefficients" table that's part of the mode.
coefficients = mdl_export.Coefficients{:, 'Estimate'}
% Create smoothed/regressed data using the model:
yFitted = coefficients(1) * exp(-coefficients(2)*X) + coefficients(3);
% Now we're done and we can plot the smooth model as a red line going through the noisy blue markers.
hold on;
plot(X, yFitted, 'r-', 'LineWidth', 2);
grid on;
title('Exponential Regression with fitnlm()', 'FontSize', fontSize);
xlabel('X', 'FontSize', fontSize);
ylabel('Y', 'FontSize', fontSize);
legendHandle = legend('Original Y', 'Fitted Y', 'Location', 'north');
legendHandle.FontSize = 30;
formulaString = sprintf('Y = %.3f * exp(-%.3f * X) + %.3f', coefficients(1), coefficients(2), coefficients(3))
text(7, 11, formulaString, 'FontSize', 25, 'FontWeight', 'bold');

% Set up figure properties:
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% Get rid of tool bar and pulldown menus that are along top of figure.
% set(gcf, 'Toolbar', 'none', 'Menu', 'none');
% Give a name to the title bar.
set(gcf, 'Name', 'Demo by ImageAnalyst', 'NumberTitle', 'Off') 

%% Fit import half-time

X = time_vec(off_frame:off_fit_end)-time_vec(off_frame);
Y = nuclearFluo_mean(off_frame:off_fit_end)';

% Plot the initial data.
fig = figure;
plot(X', Y', 'b*', 'LineWidth', 2, 'MarkerSize', 15);
grid on;

% Convert X and Y into a table, which is the form fitnlm() likes the input data to be in.
tbl = table(X', Y');
% Note how this "x" of modelfun is related to big X and big Y.
% x((:, 1) is actually X and x(:, 2) is actually Y - the first and second columns of the table.
modelfun = @(b,x) b(1) * exp(-b(2)*x(:, 1)) + b(3);  
beta0 = [-5E5 3 5E5]; % Guess values to start with.  Just make your best guess.
% Now the next line is where the actual model computation is done.
mdl_import = fitnlm(tbl, modelfun, beta0);
% Now the model creation is done and the coefficients have been determined.

% Extract the coefficient values from the the model object.
% The actual coefficients are in the "Estimate" column of the "Coefficients" table that's part of the mode.
coefficients = mdl_import.Coefficients{:, 'Estimate'}
% Create smoothed/regressed data using the model:
yFitted = coefficients(1) * exp(-coefficients(2)*X) + coefficients(3);
% Now we're done and we can plot the smooth model as a red line going through the noisy blue markers.
hold on;
plot(X, yFitted, 'r-', 'LineWidth', 2);
grid on;
title('Exponential Regression with fitnlm()', 'FontSize', fontSize);
xlabel('X', 'FontSize', fontSize);
ylabel('Y', 'FontSize', fontSize);
legendHandle = legend('Original Y', 'Fitted Y', 'Location', 'north');
legendHandle.FontSize = 30;
formulaString = sprintf('Y = %.3f * exp(-%.3f * X) + %.3f', coefficients(1), coefficients(2), coefficients(3))
text(7, 11, formulaString, 'FontSize', 25, 'FontWeight', 'bold');

% Set up figure properties:
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% Get rid of tool bar and pulldown menus that are along top of figure.
% set(gcf, 'Toolbar', 'none', 'Menu', 'none');
% Give a name to the title bar.
set(gcf, 'Name', 'Demo by ImageAnalyst', 'NumberTitle', 'Off') 


