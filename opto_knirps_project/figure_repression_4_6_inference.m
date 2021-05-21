clear
close all
clc

addpath(genpath('./lib'))

%% Initialization

projectName = 'optokni_eve4+6_MCP-GFP_Homo'; 

liveProject = LiveEnrichmentProject(projectName);
resultsRoot = [liveProject.dataPath filesep 'cpHMM_results' filesep];

% load data
load([resultsRoot 'hmm_input_output.mat'])
FigurePath = [liveProject.figurePath 'mHMM_analysis' filesep];
mkdir(FigurePath)

% color to be used
k_green = brighten([38 142 75]/256,.4);
color_green = [38 143 75]/256; % color from Jake
mRNA_red = brighten([212 100 39]/256,.2);


%% Figure: plot traces and inference results

traceNum = 181;

time = hmm_input_output(traceNum).time/60;
pState = hmm_input_output(traceNum).promoter_state;
Fluo = hmm_input_output(traceNum).fluo;
predFluo = hmm_input_output(traceNum).predicted_fluo;

TraceFig = figure;
yyaxis left
stairs(time,pState-1);
ylim([-0.1 2.1])
yyaxis right
plot(time,Fluo)
hold on
plot(time,predFluo)


