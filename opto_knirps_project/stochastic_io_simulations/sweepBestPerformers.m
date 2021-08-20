function sweepInfoBest = sweepBestPerformers(sweepInfoBest,sweepInfo)

    sweepResultsBest = struct;
    
    % initialize score array
    score_array = [];
    for f = 1:length(sweepInfoBest.fit_fields_to_use)         
        score_array = [score_array sweepInfo.(sweepInfoBest.fit_fields_to_use{f})];
    end    
    r2_flag = any(contains(sweepInfoBest.fit_fields_to_use,'R2'));
    if r2_flag
        total_score_vec = sqrt(sum(score_array,2));
        [~,score_indices] = sort(total_score_vec,'ascend');
    else
        total_score_vec = sum(score_array,2);
        [~,score_indices] = sort(total_score_vec,'descend');
    end
    n_keep = sweepInfoBest.n_raw;
    keep_indices = score_indices(1:n_keep);
    for i = 1:length(keep_indices)
        sweepResultsBest(i).bestScore = total_score_vec(keep_indices(i));
        sweepResultsBest(i).bestInd = keep_indices(i);
        sweepResultsBest(i).param_val_vec = sweepInfo.param_val_vec(keep_indices(i),:);    
    end       
    
    % truncate sweep structure                          
    sweepInfoBest.keep_prediction_flag = true;    
    sweepInfoBest.nIterations = n_keep;    
    sweepInfoBest.NumWorkers = min([sweepInfo.NumWorkers ceil(sweepInfoBest.nIterations/4)]);
    sweepInfoBest.bestIndVec = [sweepResultsBest.bestInd]';
    sweepInfoBest.bestScoreVec = [sweepResultsBest.bestScore]';
    
%     sweepResultsBest = sweepResultsBest(ia);  
    [sweepInfoBest, sweepResultsBest] = initializeFitFields(sweepInfoBest,sweepResultsBest);
    % call sweep script
    sweepResultsBest = sweep_par_loop_v3(sweepInfoBest,sweepResultsBest); 
        
    % recombine 
    fnames = fieldnames(sweepResultsBest);
    for f = 1:length(fnames)
        test_var = sweepResultsBest(1).(fnames{f});
        if size(test_var,1)==1
            sweepInfoBest.(fnames{f}) = vertcat(sweepResultsBest.(fnames{f}));
        else
            sweepInfoBest.(fnames{f}) = cat(3,sweepResultsBest.(fnames{f}));
        end
    end
    
    clear sweepResultsBest
    
    %     if any(contains(sweepInfoBest.fit_fields_to_use,'R2'))
%         r2_flag = true;
%         score_array = sqrt(score_array);
% %         score10 = prctile(score_array,10);
% %         score_array_norm = score_array .* (max(score10)./score10);
%         score_array_norm = score_array ./ mean(score_array);
%     else
%         score90 = prctile(score_array,90);
%         score_array_norm = score_array .* (min(score90)./score90);
%     end

    % find the N best networks
%     for i = 1:size(sweepInfoBest.weight_array,1)
%         
%         % find best parameter set
%         if ~r2_flag
%             % first calculate weighted score
%             score_vec = zeros(size(sweepInfo.param_val_vec,1),1);        
%             for f = 1:length(sweepInfoBest.fit_fields_to_use)         
%                 score_vec = score_vec + sweepInfoBest.weight_array(i,f)*score_array_norm(:,f);
%             end
%             score_vec = score_vec ./ sum(sweepInfoBest.weight_array(i,:));
%             [sweepResultsBest(i).bestScore, sweepResultsBest(i).bestInd] = nanmax(score_vec);     
%         else
%             score_vec = zeros(size(sweepInfo.param_val_vec,1),1);        
%             for f = 1:length(sweepInfoBest.fit_fields_to_use)         
%                 score_vec = score_vec + sweepInfoBest.weight_array(i,f)*log(score_array_norm(:,f));
%             end
%             score_vec = exp(score_vec ./ sum(sweepInfoBest.weight_array(i,:)));
%             % find best parameter set
%             [sweepResultsBest(i).bestScore, sweepResultsBest(i).bestInd] = nanmin(score_vec);     
%         end
%         sweepResultsBest(i).param_val_vec = sweepInfo.param_val_vec(sweepResultsBest(i).bestInd,:);        
%     end
    % There will be duplicates, so let's find unique indices of the best
    % performers
%     [unique_indices, ia] = unique([sweepResultsBest.bestInd])