% script to compile results from cpHMM inference
clear
close all
addpath(genpath('utilities'))

% projectNameCell = {'EveGtSL','EveGtSL-S1Null','EveWt','EveS1Null'};%};
projectNameCell = {'20210928_Oct4_raw_traces_Oct4_dark_control_raw_trace','20210928_Oct4_raw_traces_Oct4_opto_raw_trace'};
low_order_flag = 0;
% resultsRoot = 'S:\Nick\Dropbox\InductionLogic\';

for p = 2%:length(projectNameCell)
    
    % set project to analyze 
    projectName = projectNameCell{p};

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
        resultsDir = [resultsRoot filesep projectNameCell{p} filesep 'cpHMM_results' filesep];
    end
    % get list of all inference subdirectories. By default, we'll generate
    % summaries for all non-empty inference sub-directory
    infDirList = dir([resultsDir 'w*']);

    % iterate through the directories and compile the results
    for inf = 2%:length(infDirList)

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

        % get unique list of groups
        groupID_vec = [inferenceResults.groupID];
        groupID_index = inferenceOptions.indexInfo.indexVecUnique;%unique(groupID_vec);
        n_groups_orig = length(inferenceOptions.indexInfo.intensity_group_vec);
        
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
        compiledResults.sigma_vec_mean = NaN(size(groupID_index));
        compiledResults.sigma_vec_ste = NaN(size(groupID_index));
        compiledResults.r_array_mean = NaN(nStates,length(groupID_index));
        compiledResults.r_array_ste = NaN(nStates,length(groupID_index));
        compiledResults.pi0_array_mean = NaN(nStates,length(groupID_index));
        compiledResults.pi0_array_ste = NaN(nStates,length(groupID_index));
        compiledResults.R_array_mean = NaN(nStates,nStates,length(groupID_index));
        compiledResults.R_array_ste = NaN(nStates,nStates,length(groupID_index));
        compiledResults.fluo_mean = NaN(size(groupID_index));
        compiledResults.fluo_ste = NaN(size(groupID_index));
        compiledResults.time_vec_mean = NaN(size(groupID_index));
        compiledResults.time_vec_ste = NaN(size(groupID_index));
        compiledResults.outlier_frac = NaN(size(groupID_index));
        
        % cell arrays to store full results
        compiledResults.particle_ids = cell(size(groupID_index));
        compiledResults.init_results = cell(size(groupID_index));
        compiledResults.dur_results = cell(size(groupID_index));
        compiledResults.freq_results = cell(size(groupID_index));
        compiledResults.r_results = cell(size(groupID_index));
        compiledResults.R_results = cell(size(groupID_index));
        compiledResults.A_results = cell(size(groupID_index));
        compiledResults.sigma_results = cell(size(groupID_index));
        compiledResults.pi0_results = cell(size(groupID_index));
        compiledResults.logL_results = cell(size(groupID_index));
        compiledResults.outlier_flags = cell(size(groupID_index));
        
        % Pull grouping variable info
        compiledResults.timeBins = inferenceOptions.timeBins;
        compiledResults.apBins = inferenceOptions.apBins;        
        compiledResults.additionalGroupVar = inferenceOptions.AdditionalGroupingVariable;       
        compiledResults.skippedFlag = true(size(compiledResults.additionalGroupVar));
        
        % longform group vectors (1 per inference group)
        compiledResults.apGroupVec = inferenceOptions.indexInfo.ap_group_vec;
        compiledResults.timeGroupVec = inferenceOptions.indexInfo.time_group_vec;
        compiledResults.additionalGroupVec = inferenceOptions.indexInfo.additional_group_vec;         
        
        if inferenceOptions.FluoBinFlag
          compiledResults.spot_intensity_vec = inferenceOptions.indexInfo.intensity_value_vec;%(ismember(1:n_groups_orig,groupID_index));
        elseif inferenceOptions.ProteinBinFlag
          compiledResults.protein_intensity_vec = inferenceOptions.indexInfo.intensity_value_vec;%(ismember(1:n_groups_orig,groupID_index));
        end                
        
        % iterate through groups
        for g = 1:length(groupID_index)
          
            d_ids = find(groupID_vec==groupID_index(g));
            
            if ~isempty(d_ids)
                particle_ids_temp = [];
                init_array = NaN(size(d_ids));
                freq_array = NaN(size(d_ids));
                dur_array = NaN(size(d_ids)); 
                sigma_array = NaN(size(d_ids)); 
                logL_array = NaN(size(d_ids)); 
                outlier_array = NaN(size(d_ids)); 
                R_array = NaN(nStates,nStates,length((d_ids)));
                A_array = NaN(nStates,nStates,length((d_ids)));
                r_array = NaN(nStates,length((d_ids)));
                pi0_array = NaN(nStates,length((d_ids)));
                fluo_array = NaN(size(d_ids)); 
                time_array = NaN(size(d_ids)); 

                for i = 1:length(d_ids)

                    particle_ids_temp = [particle_ids_temp inferenceResults(d_ids(i)).particle_ids];

                    % enforce rank ordering of states by initiation rate
                    [r,ri] = sort(inferenceResults(d_ids(i)).r);
                    if min(r) < 0 && abs(min(r)) >= 0.25*max(r)% && inf > 1
                      warning('negative initiation rate(s) found')
                      continue
                    end
                    A = inferenceResults(d_ids(i)).A_mat(ri,ri);

                    % convert to transition rate matrix
                    R = logm(A) / Tres;
                    if ~isreal(R) || sum(R(:)<0) > nStates
                        if low_order_flag 
                            out = prob_to_rate_fit_sym(A, Tres, 'low_order', .005, 0.05);            
                        else
                            out = prob_to_rate_fit_sym(A, Tres, 'gen', .005, 0.05);            
                        end
                        R = out.R_out;     
                    end
                    [V,D] = eig(R);
                    [~,di] = max(real(diag(D)));
                    ss_vec = V(:,di) / sum(V(:,di));

                    % generate 2 state initiation rates
                    if nStates == 3
                        init_array(i) = (r(2) * ss_vec(2) + r(3) * ss_vec(3)) / (ss_vec(2)+ss_vec(3));   
                        if sum(r.*ss_vec) > 1.1*init_array(i)
                            error('nontrivial "off" state initiation')
                        end            
                    elseif nStates == 2
                        init_array(i) = r(2);
                    else
                        error('Systems with more than 3 states are not currently supported')
                    end
                    
                    % burst freq
                    freq_array(i) = -R(1,1);

                    % burst dur
                    dur_array(i) = -1/R(1,1) * (1/ss_vec(1) - 1);

                    % fluorescence
                    fluo_array(i) = nanmean([inferenceResults(d_ids(i)).fluo_data{:}]);
                    time_array(i) = nanmean([inferenceResults(d_ids(i)).time_data{:}]);

                    % basic inference outputs
                    A_array(:,:,i) = A;
                    R_array(:,:,i) = R;
                    r_array(:,i) = r;
                    sigma_array(i) = sqrt(inferenceResults(d_ids(i)).noise);
                    logL_array(i) = inferenceResults(d_ids(i)).max_logL;
                    pi0_array(:,i) = sqrt(inferenceResults(d_ids(i)).pi0);
                end

                compiledResults.particle_ids{g} = unique(particle_ids_temp);

                % check for outliers
                [freq_array_filt,freq_flags] = rmoutliers(freq_array);
                [dur_array_filt,dur_flags] = rmoutliers(dur_array);
                [init_array_filt,init_flags] = rmoutliers(init_array);

                outlier_filter = freq_flags|dur_flags|init_flags;

                % fraction that were outliers
                compiledResults.outlier_frac(g) = mean(outlier_filter);

                % save full unfiltered results
                compiledResults.freq_results{g} = freq_array;
                compiledResults.dur_results{g} = dur_array;
                compiledResults.init_results{g} = init_array;
                compiledResults.outlier_flags{g} = outlier_filter;

                % calculate average and ste
                compiledResults.init_vec_mean(g) = nanmean(init_array_filt)*60;
                compiledResults.init_vec_ste(g) = nanstd(init_array_filt)*60;

                compiledResults.freq_vec_mean(g) = nanmean(freq_array_filt)*60;
                compiledResults.freq_vec_ste(g) = nanstd(freq_array_filt)*60;

                compiledResults.dur_vec_mean(g) = nanmean(dur_array_filt)/60;
                compiledResults.dur_vec_ste(g) = nanstd(dur_array_filt)/60;

                compiledResults.fluo_mean(g) = nanmean(fluo_array(~outlier_filter));
                compiledResults.fluo_ste(g) = nanstd(fluo_array(~outlier_filter));
                
                compiledResults.time_vec_mean(g) = nanmean(time_array(~outlier_filter));
                compiledResults.time_vec_ste(g) = nanstd(time_array(~outlier_filter));

                compiledResults.sigma_vec_mean(g) = nanmean(sigma_array(~outlier_filter));
                compiledResults.sigma_vec_ste(g) = nanstd(sigma_array(~outlier_filter));

                compiledResults.r_array_mean(:,g) = nanmean(r_array(:,~outlier_filter),2);
                compiledResults.r_array_ste(:,g) = nanstd(r_array(:,~outlier_filter),[],2);

                compiledResults.pi0_array_mean(:,g) = nanmean(pi0_array(:,~outlier_filter),2);
                compiledResults.pi0_array_ste(:,g) = nanstd(pi0_array(:,~outlier_filter),[],2);

                compiledResults.R_array_mean(:,:,g) = nanmean(R_array(:,:,~outlier_filter),3);
                compiledResults.R_array_ste(:,:,g) = nanstd(R_array(:,:,~outlier_filter),[],3);

                A_mean_raw = nanmean(A_array(:,:,~outlier_filter),3);             
                compiledResults.A_array_mean(:,:,g) = A_mean_raw./sum(A_mean_raw);
                compiledResults.A_array_ste(:,:,g) = nanstd(A_array(:,:,~outlier_filter),[],3)./sum(A_mean_raw);

                % basic outputs
                compiledResults.A_results{g} = A_array;
                compiledResults.R_results{g} = R_array;
                compiledResults.r_results{g} = r_array;
                compiledResults.pi0_results{g} = pi0_array;
                compiledResults.sigma_results{g} = sigma_array;
                compiledResults.logL_results{g} = logL_array;
                compiledResults.skippedFlag(g) = false;
            end
        end  

        % save
        save([infDirList(inf).folder filesep 'compiledResults_' infDirList(inf).name '.mat'],'compiledResults')
    end
end

            