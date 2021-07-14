function simInfo = calculate_cumulative_dist(simInfo,F_min)

gillespie = simInfo.gillespie;

simInfo.p_on_array = nanmean(gillespie.fluo_ms2_array>F_min,2);

% if contains(simInfo.simType,'off')
%     % initialize array
%     on_off_array = ones(simInfo.seq_length,simInfo.n_traces,length(simInfo.F_min));
%     for f = 1:length(simInfo.F_min)
%         for n = 1:simInfo.n_traces
%             last_i = find(gillespie.fluo_ms2_array(:,n)>simInfo.F_min(f),1,'last');
%             on_off_array(last_i+1:end,n,f) = 0;
%         end
%     end
%     % calculate cumulative curves
%     simInfo.p_on_array = reshape(nanmean(on_off_array,2),simInfo.seq_length,length(simInfo.F_min));
% elseif contains(simInfo.simType,'on')    
%     % initialize array
%     on_off_array = zeros(simInfo.seq_length,simInfo.n_traces,length(simInfo.F_min));
%     for f = 1:length(simInfo.F_min)
%         for n = 1:simInfo.n_traces
%             first_i = find(gillespie.fluo_ms2_array(:,n)>simInfo.F_min(f),1);
%             on_off_array(first_i:end,n,f) = 1;
%         end
%     end
%     % calculate cumulative curves
%     simInfo.p_on_array = reshape(nanmean(on_off_array,2),simInfo.seq_length,length(simInfo.F_min));
% end