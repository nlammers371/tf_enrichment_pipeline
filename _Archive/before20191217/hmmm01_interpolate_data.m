% hmmm01_interpolate_data(project)
%
% DESCRIPTION
% Script to generate samples of local protein contentration at active gene
% loci and selected control locations
%
% ARGUMENTS
% project: master ID variable 
%
% OUTPUT: nucleus_struct: 

function nucleus_struct = hmmm01_interpolate_data(project,varargin)

dropboxFolder =  'E:\Nick\LivemRNA\Dropbox\';
dataPath = [dropboxFolder 'ProcessedEnrichmentData\' project '/'];
for i = 1:numel(varargin)
    if strcmpi(varargin{i},'dropboxFolder')
        dataPath = [varargin{i+1} '\ProcessedEnrichmentData\' project '/'];
    end
    if ischar(varargin{i})
        if ismember(varargin{i},{'minDp','TresInterp'})
            eval([varargin{i} '=varargin{i+1};']);
        end
    end
end
load([dataPath '/nucleus_struct.mat'],'nucleus_struct')
% estimate interpolation time
med_time = nanmedian(diff([nucleus_struct.time]));
if exist('TresInterp')~= 1
    TresInterp = round(med_time);
end
interpGrid = 0:TresInterp:60*60;
%%% Cleaning Parameters
% big_jump1 = prctile([nucleus_struct.fluo],99);
% jump_threshold1 = big_jump1/1.5; % this should be a conservative threshold for single time step increase
minDP = nucleus_struct(1).minDP;
for i = 1:length(nucleus_struct) 
    temp = nucleus_struct(i);
    trace = temp.fluo; %Load full trace, including intervening NaN's    
    pt_time = temp.time;      
%     trace1(pt_time < minTime) = NaN;
    quality_flag = temp.qc_flag; % indicates whether trace suitable for inference
    
    if sum(~isnan(trace)) == 0 
        t_start = 0;
        t_stop = -1;
    else
        start_i = find(interpGrid>=min(pt_time(~isnan(trace))),1);
        t_start = interpGrid(start_i);
        stop_i = find(interpGrid<=max(pt_time(~isnan(trace))),1,'last');
        t_stop = interpGrid(stop_i);
    end
    time_interp = t_start:TresInterp:t_stop;
    
    if quality_flag==1        
        %Null assumption is that all clusters of 6 or more NaNs are 0s. Smaller
        %clusters are assumed to have been missed nonzero dps
        trace_nans = isnan(trace);      
        %Look for clusters of 6 or more NaNs
        kernel = [1,1,1,1,1];
        tn_conv = conv(kernel,trace_nans);
        tn_conv = tn_conv(3:end-2);
        z_ids = find(tn_conv==5);
        z_ids = unique([z_ids-1 z_ids z_ids+1]); % get set of z_ids    
        trace(z_ids) = 0; % set clusters to zeros    
        trace(trace<0) = 0; % deal with negative values    
        % find single dp "blips". These will be replaced via interpolation
        % interpolate remaining NaNs    
        query_points1 = pt_time(isnan(trace));
        interp_t1 = pt_time(~isnan(trace));
        interp_f1 = trace(~isnan(trace)); 
        new_f1 = interp1(interp_t1,interp_f1,query_points1);      
        trace(ismember(pt_time,query_points1)) = new_f1;        
        % Interpolate to standardize spacing    
        trace_interp = interp1(pt_time,trace,time_interp); 
    else
        trace_interp = NaN(size(time_interp));
        time_interp = NaN(size(time_interp));           
    end  
    if sum(~isnan(trace_interp)) < minDP && quality_flag == 1
        error('qc_flag error')
    end
    nucleus_struct(i).fluo_interp = trace_interp;    
    nucleus_struct(i).time_interp = time_interp;
%     nucleus_struct(i).inference_flag = quality_flag;  
    nucleus_struct(i).TresInterp = TresInterp;
end
% save
save([dataPath '/nucleus_struct.mat'],'nucleus_struct')