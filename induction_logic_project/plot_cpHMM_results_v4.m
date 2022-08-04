% script to compile results from cpHMM inference
clear
close all
addpath(genpath('utilities'))

% projectNameCell = {'EveGtSL','EveGtSL-S1Null','EveWt','EveS1Null'};%};
doseNameCell = {'20220701_Oct4_dose_dose_response_mid_200-400a',...
                  '20220701_Oct4_dose_dose_response_low_0-200a','20220701_Oct4_dose_dose_response_high_400-1200a'};
low_order_flag = 0;
iter = 1;
% resultsRoot = 'S:\Nick\Dropbox\InductionLogic\';
dose_struct = struct;
for p = 1:length(doseNameCell)
    
    % set project to analyze 
    projectName = doseNameCell{p};

    % get path to results
    try
        liveProject = LiveEnrichmentProject(projectName);
        resultsDir = [liveProject.dataPath 'cpHMM_results' filesep];
    catch
%         resultsRoot = 'S:/Nick/Dropbox/ProcessedEnrichmentData/';
        if isdir('C:\Users\nlamm\Dropbox (Personal)\')
            resultsRoot = 'C:\Users\nlamm\Dropbox (Personal)\ProcessedEnrichmentData';
        else
           resultsRoot = 'S:\Nick\Dropbox (Personal)\ProcessedEnrichmentData';
        end
        resultsDir = [resultsRoot filesep doseNameCell{p} filesep 'cpHMM_results' filesep];
    end
    
    % get list of all inference subdirectories. By default, we'll generate
    % summaries for all non-empty inference sub-directory
    infDirList = dir([resultsDir 'w*']);

    % iterate through the directories and compile the results
    for inf = 1:length(infDirList)
        load([infDirList(inf).folder filesep 'compiledResults_' infDirList(inf).name '.mat'],'compiledResults');
        fnames = fieldnames(compiledResults);
        for f = 1:length(fnames)
            dose_struct(iter).(fnames{f}) = compiledResults.(fnames{f});
        end
        iter = iter + 1;
    end
end

% Now load opto results

% projectNameCell = {'EveGtSL','EveGtSL-S1Null','EveWt','EveS1Null'};%};
optoNameCell = {'20220701_Oct4_opto_opto_contro','20220701_Oct4_opto_opto_LEXY-YA'};

iter = 1;
% resultsRoot = 'S:\Nick\Dropbox\InductionLogic\';
opto_struct = struct;
for p = 1:length(optoNameCell)
    
    % set project to analyze 
    projectName = optoNameCell{p};

    % get path to results
    try
        liveProject = LiveEnrichmentProject(projectName);
        resultsDir = [liveProject.dataPath 'cpHMM_results' filesep];
    catch
%         resultsRoot = 'S:/Nick/Dropbox/ProcessedEnrichmentData/';
        if isdir('C:\Users\nlamm\Dropbox (Personal)\')
            resultsRoot = 'C:\Users\nlamm\Dropbox (Personal)\ProcessedEnrichmentData';
        else
           resultsRoot = 'S:\Nick\Dropbox (Personal)\ProcessedEnrichmentData';
        end
        resultsDir = [resultsRoot filesep optoNameCell{p} filesep 'cpHMM_results' filesep];
    end
    
    % get list of all inference subdirectories. By default, we'll generate
    % summaries for all non-empty inference sub-directory
    infDirList = dir([resultsDir 'w*']);

    % iterate through the directories and compile the results
    for inf = 1:length(infDirList)
        load([infDirList(inf).folder filesep 'compiledResults_' infDirList(inf).name '.mat'],'compiledResults');
        fnames = fieldnames(compiledResults);
        for f = 1:length(fnames)
            opto_struct(iter).(fnames{f}) = compiledResults.(fnames{f});
        end
        iter = iter + 1;
    end
end
