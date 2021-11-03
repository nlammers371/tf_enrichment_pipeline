% script to plot results from cpHMM inference
clear
close all
addpath(genpath('../utilities'))

% projectNameCell = {'EveGtSL','EveGtSL-S1Null','EveWt','EveS1Null'};%};
projectNameCell = {'20210928_Oct4_raw_traces_nz_Oct4_dark_control_raw_trace','20210928_Oct4_raw_traces_nz_Oct4_opto_raw_trace'};%};
geneName = 'Oct4';
infString = 'K2_p0_ap1_t2_f2D';
conditionCell = {'control', 'opto'};
geneCell = {'Oct4','Oct4'};
nz_flag = 1;

    
for p = 1:length(projectNameCell)
    
    % set project to analyze 
    projectName = projectNameCell{p};
    geneName = geneCell{p};
    conditionName = conditionCell{p};
    
    % get path to results
    try
        liveProject = LiveEnrichmentProject(projectName);
        resultsDir = [liveProject.dataPath 'cpHMM_results' filesep];
    catch
%         resultsRoot = 'S:/Nick/Dropbox/ProcessedEnrichmentData/';
        resultsRoot = 'C:\Users\nlamm\Dropbox (Personal)\ProcessedEnrichmentData';
        resultsDir = [resultsRoot filesep projectNameCell{p} filesep 'cpHMM_results' filesep];
    end    
    
    % get list of projects
    resultList = dir([resultsDir '*result*']);            
    name_cell = {resultList.name};
    inf_index = find(contains(name_cell,infString));        
      
    % load data
    load([resultsDir filesep resultList(inf_index).name]);

    timeGroup = compiledResults.timeGroupVec;
    
    % transfer results
    init_vec_mean = compiledResults.init_vec_mean;
    init_vec_ste = compiledResults.init_vec_ste;

    freq_vec_mean = compiledResults.freq_vec_mean;
    freq_vec_ste = compiledResults.freq_vec_ste;

    dur_vec_mean = compiledResults.dur_vec_mean;
    dur_vec_ste = compiledResults.dur_vec_ste;

    time_vec_mean = compiledResults.time_vec_mean;
    time_vec_ste = compiledResults.time_vec_ste;

    fluo_vec_mean = compiledResults.fluo_mean;
    fluo_vec_ste = compiledResults.fluo_ste;            
    
    
    % make output dataset
    result_array = [double(timeGroup') fluo_vec_mean' fluo_vec_ste' init_vec_mean' init_vec_ste' ...
          dur_vec_mean' dur_vec_ste' freq_vec_mean' freq_vec_ste'];
    result_table = array2table(result_array,'VariableNames',{'time_group','spot_intensity_mean','spot_intensity_ste',...
              'initiation_rate_mean','initiation_rate_ste','burst_duration_mean','burst_duration_ste',...
              'burst_frequency_mean','burst_frequency_ste'});
%     result_table.experiment_type = condition_key';
%     result_table = [result_table(:,end) result_table(:,1:end-1)];        
    
    writetable(result_table,[resultsDir geneName '_' conditionName '_nz' num2str(nz_flag) '_inference_results.csv'])
end


