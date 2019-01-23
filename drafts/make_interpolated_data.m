% main04_interpolate_data(project, zeiss_data,)
%
% DESCRIPTION
% Script to clean and interpolate fluorescence trace
%
% ARGUMENTS
% project: master ID variable 
%
% RawPath: Full or relative path to PreProcessed folder (pr
%               equivalent)

%
% OUTPUT: nucleus_Struct_final: compiled data set
function nucleus_struct_final = main04_interpolate_data(project)

%%%Set Cleaning and Summary Parameters
min_dp = 10; % minimum # dp acceptable for an inference trace
Tres_interp = 20; % time resolution
InterpGrid = 0:Tres_interp:60*50;
%------------------------Import Raw Trace Set------------------------%
%ID's of sets to include
% project = 'kr_eve2_reporter';

%---------------------------Set Paths-------------------------------------%
NucleusPath = ['../dat/' project '/nucleus_struct_ctrl.mat'];
OutPath = ['../dat/' project '/'];
FigPath = ['../fig/' project '/preprocessing'];
TraceSavePath = [FigPath '/traces/'];
mkdir(TraceSavePath);
mkdir(OutPath);
% Save Name
CleanNucleusName = ['inference_nuclei_' project '_dT' num2str(Tres_interp) '.mat'];

%----------------Load Traces and Perform First Pass QC--------------------%
%Load raw traces (saved in struct titled "nucleus_struct_ctrl")
load(NucleusPath);
%%% Cleaning Parameters
big_jump1 = prctile([nucleus_struct_ctrl.fluo],99);
jump_threshold1 = big_jump1/1.5; % this should be a conservative threshold for single time step increase
index_vec = 1:length(nucleus_struct_ctrl); % convenience ref vector
field_names = fieldnames(nucleus_struct_ctrl);

nucleus_struct_final = []; % Structure to store traces after first round of cleaning
pt_vec = [nucleus_struct_ctrl.ParticleID];
pt_indices = find(~isnan(pt_vec));
for i = pt_indices 
    temp = nucleus_struct_ctrl(i);
    trace1 = temp.fluo; %Load full trace, including intervening NaN's    
    if sum(~isnan(trace1)) >= min_dp    
        protein = temp.protein; % full protein trace including NaNs
        nuc = nucleus_struct_ctrl([nucleus_struct_ctrl.Nucleus] == temp.Nucleus);
        nuc = nuc([nuc.setID] == temp.setID);    
        time = temp.time;      
        quality_flag = 1;
        if sum(~isnan(trace1)) < min_dp
            quality_flag = 0;
        end
        start = find(~isnan(trace1),1);
        stop = find(~isnan(trace1),1,'last');
        trace1 = trace1(start:stop);
        t_start = InterpGrid(find(InterpGrid>=time(start),1));
        t_stop = InterpGrid(find(InterpGrid<=time(stop),1,'last'));
        %Null assumption is that all clusters of 6 or more NaNs are 0s. Single
        %,double, or triple NaNs are assumed to have been missed nonzero dps
        trace1_nans = isnan(trace1);      
        %Look for clusters of 6 or more NaNs
        kernel = [1,1,1,1,1];
        tn_conv = conv(kernel,trace1_nans);
        tn_conv = tn_conv(3:end-2);
        z_ids = find(tn_conv==5);
        z_ids = unique([z_ids-1 z_ids z_ids+1]); % get set of z_ids    
        trace1(z_ids) = 0; % set clusters to zeros    
        trace1(trace1<0) = 0; % deal with negative values          

        % find single dp "blips"
        tr_dd1 = abs([0 diff(diff(trace1)) 0]);  % 1 slice
        trace1(tr_dd1>2*jump_threshold1) = NaN;    

        % interpolate remaining NaNs    
        query_points1 = time(isnan(trace1));%InterpGrid((InterpGrid>=min(time))&(InterpGrid<=max(time)));
        interp_t1 = time(~isnan(trace1));
        interp_f1 = trace1(~isnan(trace1));
        new_f1 = interp1(interp_t1,interp_f1,query_points1);  
        trace1(ismember(time,query_points1)) = new_f1;            

        %%% flag traces with unreasonably large rises or falls    
        tr_d1 = diff(trace1);    
        if max(abs(tr_d1)) >= jump_threshold1        
            quality_flag = 0;
        end

        % Interpolate to standardize spacing
                
        time_interp = t_start:Tres_interp:t_stop;
        trace1_interp = interp1(time(start:stop),trace1,time_interp);        

        trace1_interp(isnan(trace1_interp)) = 0;
        temp.fluo_interp = trace1_interp;    
        temp.time_interp = time_interp;   
        temp.inference_flag = quality_flag;            
        % generate downsampled time and fluo fields for viterbir fitting            
    else
        temp.fluo_interp = NaN;    
        temp.time_interp = NaN;   
        temp.inference_flag = 0;
    end
    nucleus_struct_final = [nucleus_struct_final temp];    
end

%%% Add a few useful fields
for i = 1:length(nucleus_struct_final)    
    nucleus_struct_final(i).dT = Tres_interp;                
    nucleus_struct_final(i).alpha_frac = 1302/6544;
    nucleus_struct_final(i).InterpGrid = InterpGrid;    
    nucleus_struct_final(i).MeanAPOrig = mean(nucleus_struct_final(i).ap_vector);   
end
%------------------------- Clean Ellipse Set -----------------------------%
save([OutPath CleanNucleusName],'nucleus_struct_final');