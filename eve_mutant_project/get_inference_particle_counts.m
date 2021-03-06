% script to compile results from cpHMM inference
clear
close all
addpath(genpath('utilities'))

projectNameCell = {'EveGtSL','EveGtSL-S1Null','EveWt','EveS1Null'};%};

% resultsRoot = 'S:\Nick\Dropbox\InductionLogic\';

for p = 1:length(projectNameCell)
    % set project to analyze 
    projectName = projectNameCell{p};

    % get path to results
    if ~exist('resultsRoot','var')
        liveProject = LiveEnrichmentProject(projectName);
        resultsDir = [liveProject.dataPath 'cpHMM_results' filesep];
    else
        resultsDir = [resultsRoot filesep projectNameCell{p} filesep 'cpHMM_results' filesep];
    end
    % get list of all inference subdirectories. By default, we'll generate
    % summaries for all non-empty inference sub-directory
    infDirList = dir([resultsDir 'w*']);

    % iterate through the directories and compile the results
    for inf = 1:length(infDirList)

        resultsPath = [infDirList(inf).folder filesep infDirList(inf).name filesep];

        % load inference options file 
        load([resultsPath 'inferenceOptions.mat'])
        if length(fieldnames(inferenceOptions)) == 1
          inferenceOptions = inferenceOptions.inferenceOptions;
        end

        % extract key inference parameters
        Tres = inferenceOptions.Tres;
        nStates = inferenceOptions.nStates;

        % read in individual inference results
        inf_files = dir([resultsPath 'hmm_results*.mat']);
        inferenceResults = struct;
        iter = 1;
        for i = 1:numel(inf_files)
          try
              load([resultsPath inf_files(i).name])
              if ~output.skip_flag && all(~isnan(output.A))
                  fn_list = fieldnames(output);
                  for f = 1:numel(fn_list)
                      inferenceResults(iter).(fn_list{f}) = output.(fn_list{f});
                  end

                  % extract group number from name
                  gInd = strfind(inf_files(i).name,'group');
                  inferenceResults(iter).groupID = str2double(inf_files(i).name(gInd+5:gInd+7));

                  iter = iter + 1;
              end
          catch
              warning(['Could not load ' inf_files(i).name '. File may be corrupt'])
          end
        end
        
        % get total particle counts and total time point counts by
        % fluorescence bin and additional group
        fluo_bin_vec = inferenceOptions.indexInfo.indexVarArray(:,3);
        group_bin_vec = inferenceOptions.indexInfo.indexVarArray(:,4);
        
        % initialize
        particle_count_vec = NaN(size(group_bin_vec));
        tp_count_vec = NaN(size(group_bin_vec));
        
        % get master indexing vectors
        fluo_bin_index = [inferenceResults.intensityBin];
        group_bin_index = [inferenceResults.additionalBin];
        
        for f = 1:length(fluo_bin_vec)
            % filter
            inf_filter = fluo_bin_index==fluo_bin_vec(f)&group_bin_index==group_bin_vec(f);
            if any(inf_filter)
                % get filtered vectors
                pt_vec_temp = [inferenceResults(inf_filter).particle_ids];
                f_vec_temp = vertcat(inferenceResults(inf_filter).fluo_data);
                % get unique list
                [pt_vec_u,u_ids] = unique(pt_vec_temp);
                f_vec_u = f_vec_temp(u_ids);
                % get counts
                particle_count_vec(f) = length(pt_vec_u);
                tp_count_vec(f) = length([f_vec_u{:}]);
            else
                particle_count_vec(f) = 0;
                tp_count_vec(f) = 0;
            end
        end
            
        count_array = [group_bin_vec fluo_bin_vec particle_count_vec tp_count_vec];
        count_table = array2table(count_array,'VariableNames',{inferenceOptions.AdditionalGroupingVariable, 'fluo_bin','n_nuclei','n_time_points'});
        writetable(count_table,[infDirList(inf).folder filesep 'particle_counts_' infDirList(inf).name '.csv'])
    end
end