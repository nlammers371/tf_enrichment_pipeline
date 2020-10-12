% script to compile results from cpHMM inference
clear
close all
addpath(genpath('../utilities'))

projectNameCell = {'EveS1Null','EveGtSL','EveGtSL-S1Null','EveWT'};


for p = 1:length(projectNameCell)
    % set project to analyze 
    projectName = projectNameCell{p};

    % get path to results
    liveProject = LiveProject(projectName);
    resultsDir = [liveProject.dataPath filesep 'cpHMM_results' filesep];

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
              if ~output.skip_flag
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

        % get unique list of groups
        groupID_vec = [inferenceResults.groupID];
        groupID_index = unique(groupID_vec);

        % initialize structure to store results
        compiledResults = struct;
        compiledResults.inferenceOptions = inferenceOptions;
        compiledResults.groupID_vec = groupID_vec;
        compiledResults.groupID_index = groupID_index;
        
        % calculate average initiation rate, burst freq, and burst duration
        compiledResults.init_vec_mean = NaN(size(groupID_index));
        compiledResults.init_vec_ste = NaN(size(groupID_index));
        compiledResults.freq_vec_mean = NaN(size(groupID_index));
        compiledResults.freq_vec_ste = NaN(size(groupID_index));
        compiledResults.dur_vec_mean = NaN(size(groupID_index));
        compiledResults.dur_vec_ste = NaN(size(groupID_index));
        compiledResults.fluo_mean = NaN(size(groupID_index));
        compiledResults.fluo_ste = NaN(size(groupID_index));
        compiledResults.particle_ids = cell(size(groupID_index));
                
        % iterate through groups
        for g = 1:length(groupID_index)
          
            d_ids = find(groupID_vec==groupID_index(g));
            
            particle_ids_temp = [];
            init_array = NaN(size(d_ids));
            freq_array = NaN(size(d_ids));
            dur_array = NaN(size(d_ids)); 
            fluo_array = NaN(size(d_ids)); 

            for i = 1:length(d_ids)
              
                particle_ids_temp = [particle_ids_temp inferenceResults(d_ids(i)).particle_ids];
                
                % enforce rank ordering of states by initiation rate
                [r,ri] = sort(inferenceResults(d_ids(i)).r);
                A = inferenceResults(d_ids(i)).A_mat(ri,ri);

                % convert to transition rate matrix
                R = logm(A) / Tres;
                if ~isreal(R) || sum(R(:)<0) > nStates
                    out = prob_to_rate_fit_sym(A, Tres, 'gen', .005, 1);            
                    R = out.R_out;     
                end
                [V,D] = eig(R);
                [~,di] = max(diag(D));
                ss_vec = V(:,di) / sum(V(:,di));

                % generate 2 state initiation rates
                init_array(i) = (r(2) * ss_vec(2) + r(3) * ss_vec(3)) / (ss_vec(2)+ss_vec(3));   
                if sum(r.*ss_vec) > 1.1*init_array(i)
                    warning('nontrivial "off" state initiation')
                end            

                % burst freq
                freq_array(i) = -R(1,1);

                % burst dur
                dur_array(i) = -1/R(1,1) *(1/ss_vec(1) - 1);
                
                % fluorescence
                fluo_array(i) = nanmean([inferenceResults(d_ids(i)).fluo_data{:}]);
            end
            
            compiledResults.particle_ids{g} = unique(particle_ids_temp);
            
            % calculate average and ste
            compiledResults.init_vec_mean(g) = nanmean(init_array)*60;
            compiledResults.init_vec_ste(g) = nanstd(init_array)*60;

            compiledResults.freq_vec_mean(g) = nanmean(freq_array)*60;
            compiledResults.freq_vec_ste(g) = nanstd(freq_array)*60;

            compiledResults.dur_vec_mean(g) = nanmean(dur_array)/60;
            compiledResults.dur_vec_ste(g) = nanstd(dur_array)/60;
            
            compiledResults.fluo_mean(g) = nanmean(fluo_array);
            compiledResults.fluo_ste(g) = nanstd(fluo_array);
        end  

        % save
        save([infDirList(inf).folder filesep 'compiledResults_' infDirList(inf).name '.mat'],'compiledResults')
    end
end

            