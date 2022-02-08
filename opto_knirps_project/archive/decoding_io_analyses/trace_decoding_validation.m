% script to test validity of viterbi decoding method
% clear
close all

addpath(genpath('../utilities'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%% Load experimental data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

resultsRoot = 'S:\Nick\Dropbox\ProcessedEnrichmentData\parameterSweeps\';
if ~isfolder(resultsRoot)
  resultsRoot = 'C:\Users\nlamm\Dropbox (Personal)\ProcessedEnrichmentData\parameterSweeps\';
end

% mak figure directory  
date_string = '22-Sep-2021';
resultsPath = [resultsRoot date_string filesep];

try
    FigurePath = ['S:\Nick\Dropbox (Personal)\LocalEnrichmentFigures\DecodedTraceAnalyses' filesep date_string filesep];
catch
    FigurePath = ['C:\Users\nlamm\Dropbox (Personal)\LocalEnrichmentFigures\PipelineOutput\DecodedTraceAnalyses' filesep date_string filesep];
end
mkdir(FigurePath)

paramIncVec = [25 25 16 16];
simTypeCell = {'koff_only_2','kon_only_2'};
tfDependentParamCell = {'koff','kon'};
master_struct = struct;

% load sweep results
for s = 1:length(simTypeCell)
    simType = simTypeCell{s};
    nInc = paramIncVec(s);
    readPath = [resultsPath 'sweeps_n' num2str(nInc) filesep];
    % load sweep result
    load([readPath 'sweepInfo_' simType '.mat'],'sweepInfo');    
    % store
    master_struct(s).sweepInfo = sweepInfo;
end

% load raw data
load([resultsRoot 'io_ref_ra.mat'],'io_ref_ra')
load([resultsRoot 'io_ref_wt.mat'],'io_ref_wt')
load([resultsRoot 'io_ref_ON.mat'],'io_ref_ON')

% obtain predicted trends for best systems from sim type
n_traces = 250;
n_best = 10;

for s = 1:length(master_struct)
    % extract
    sweepInfo = master_struct(s).sweepInfo;
    simType = sweepInfo.simType;
    
    % calculate aggregate score
    total_score = -sqrt(sweepInfo.fluo_time_ON_fit.^2 + sweepInfo.fluo_time_fit.^2 + sweepInfo.ra_fit.^2 + sweepInfo.still_on_fit.^2);

    % find best overall performers
    [~,score_ids] = sort(total_score,'descend');
    best_i_list = score_ids(1:n_best);
    
    % now find best parameter-specific scores
    [~,best_ft_i] = nanmax(sweepInfo.fluo_time_fit);
    [~,best_ft_ON_i] = nanmax(sweepInfo.fluo_time_ON_fit);
    [~,best_pon_i] = nanmax(sweepInfo.still_on_fit);
    [~,best_ra_i] = nanmax(sweepInfo.ra_fit);
    
    % get corresponding parameters
    param_vec = sweepInfo.param_val_vec(best_i_list(1),:);
    
    % run sweep for selected networks
    sweepInfoBest = io_sweep_wrapper(resultsRoot,2,sweepInfo.simType,param_vec,true,'n_traces',n_traces);
    
    % store
    master_struct(s).sweepInfoBest = sweepInfoBest;
    master_struct(s).best_i_list = best_i_list;
    master_struct(s).best_ft_i = best_ft_i;
    master_struct(s).best_ra_i = best_ra_i;
    master_struct(s).best_pon_i = best_pon_i;
end

%% Now conduct viterbi fits to see if we can recover the true trends
% get parameters to use
fitParameters = struct;
fitParameters = getMarkovSystemInfo(fitParameters);
fitParameters.A = expm(fitParameters.R2_orig*fitParameters.deltaT);
eps = 1e-4; % not actually used
A_log = log(fitParameters.A);
v = double(fitParameters.r2);
sigma = fitParameters.noise;
pi0_log = log(fitParameters.pi0');
nStates = size(A_log,1);
nSteps = fitParameters.memory;
alpha = fitParameters.t_MS2;
nWorkersMax = 24;

% obtain subset of valid traces  
for m = 1:length(master_struct)
                                 
    ms2_array = master_struct(m).sweepInfoBest.ms2_traces_true_wt;
    ms2_array(isnan(ms2_array)) = 0;
    kni_array = master_struct(m).sweepInfoBest.knirps_traces_wt;
    
    kni_values = cell(size(ms2_array,2),1);                
    fluo_values = cell(size(ms2_array,2),1);  
    time_values = cell(size(ms2_array,2),1);  
    for i = 1:size(ms2_array,2)       
        last_i = size(ms2_array,1);%find((ms2_array(:,i))>0,1,'last');
        fluo_values{i} = ms2_array(1:last_i,i)';    
        kni_values{i} = kni_array(1:last_i,i)';   
        time_values{i} = master_struct(m).sweepInfoBest.time_axis_wt(1:last_i)';   
    end    

    % initialize pool 
    p = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(p)
        parpool(nWorkersMax);
    elseif p.NumWorkers~= nWorkersMax
        delete(p);
        parpool(nWorkersMax);
    end      

    disp('conducting viterbi trace fits...')

    % initialize arrays to store results
    fluo_viterbi_array = cell(1,size(ms2_array,2));    
    path_viterbi_array = cell(1,size(ms2_array,2));    
    parfor f = 1:length(fluo_values)
        viterbi_out = viterbi (fluo_values{f}, v, sigma, pi0_log, A_log, nStates, nSteps, alpha);
        fnames = fieldnames(viterbi_out);
        
        fluo_viterbi_array{f} = viterbi_out.fluo_viterbi;
        path_viterbi_array{f} = viterbi_out.z_viterbi;       
    end
    
    master_struct(m).sweepInfoBest.path_viterbi = path_viterbi_array;
    master_struct(m).sweepInfoBest.fluo_viterbi = fluo_viterbi_array;
    
    % Attempt soft decoding
    local_em_outputs = local_em_MS2_reduced_memory (fluo_values, ...
                                v, sigma, pi0_log, A_log, nStates, nSteps, alpha, 1, eps);
    %
    
    master_struct(m).state_transition_cell = local_em_outputs.soft_struct.p_zz_log_soft;
    master_struct(m).state_prob_cell = local_em_outputs.soft_struct.p_z_log_soft;
    
    master_struct(m).fluo_cell = fluo_values;
    master_struct(m).kni_cell = kni_values;
    master_struct(m).time_cell = time_values;
end    

%% Analyze results
knirps_bins = linspace(0,16,30);
min_time = 5*60;

for m = 1:2
    kni_vec = [master_struct(m).kni_cell{:}];        
    time_vec = [master_struct(m).time_cell{:}];        
    
    % extract viterbi state vec and build on/off switch flags to indicate when
    % next observation is a state change
    state_vec = [master_struct(m).sweepInfoBest.path_viterbi{:}]-1;
    off_switch_flags = NaN(size(state_vec));
    on_switch_flags = NaN(size(state_vec));

    off_switch_probs = NaN(size(state_vec));
    on_switch_probs = NaN(size(state_vec));
    
    off_state_probs = NaN(size(state_vec));
    on_state_probs = NaN(size(state_vec));
    
    start_i = 1;
    for s = 1:length(master_struct(m).sweepInfoBest.path_viterbi)
        v_path = master_struct(m).sweepInfoBest.path_viterbi{s}-1;    

        d_path = [diff(v_path) NaN];
        on_flags = NaN(size(v_path));
        on_flags(v_path==0) = 0;
        on_flags(d_path==1&v_path==0) = 1;
        on_switch_flags(start_i:start_i+length(v_path)-1) = on_flags;

        off_flags = NaN(size(v_path));
        off_flags(v_path==1) = 0;
        off_flags(d_path==-1&v_path==1) = 1;
        off_switch_flags(start_i:start_i+length(v_path)-1) = off_flags;
       
        
        %  now do soft-decoded version
        pzz = exp(master_struct(m).state_transition_cell{s});
        pz = exp(master_struct(m).state_prob_cell{s});
        
        p_turn_on = [(reshape(pzz(2,1,:),[],1))' NaN];
        p_turn_off = [(reshape(pzz(1,2,:),[],1))' NaN];
        
        on_switch_probs(start_i:start_i+length(v_path)-1) = p_turn_on;
        off_switch_probs(start_i:start_i+length(v_path)-1) = p_turn_off;
        
        on_state_probs(start_i:start_i+length(v_path)-1) = pz(2,:);
        off_state_probs(start_i:start_i+length(v_path)-1) = pz(1,:);
        
        start_i = start_i+length(v_path);
    end   
            
    kni_bin_vec = discretize(kni_vec,knirps_bins);    

    % initialize arrays to store MF maps
    kni_state_array = NaN(length(knirps_bins)-1,1);
    kni_on_array = NaN(length(knirps_bins)-1,1);
    kni_off_array = NaN(length(knirps_bins)-1,1);
    kni_on_prob_array = NaN(length(knirps_bins)-1,1);
    kni_off_prob_array = NaN(length(knirps_bins)-1,1);

    minTime = 10*60;
    % AP vs TIME and AP vs Knirps                  
    for k = 1:length(knirps_bins)-1          
        kni_filter = kni_bin_vec>=k-1 & kni_bin_vec<=k+1 & time_vec>minTime;
        mak_filter = kni_filter;
        if sum(mak_filter)>=25
            kni_state_array(k) = nanmean(state_vec(mak_filter));
%             kni_on_prob_array(k) = nansum(on_switch_probs(mak_filter).*off_state_probs(mak_filter))/nansum(off_state_probs(mak_filter));
%             kni_off_prob_array(k) = nansum(off_switch_probs(mak_filter).*on_state_probs(mak_filter))/nansum(on_state_probs(mak_filter));
            kni_on_prob_array(k) = nansum(on_switch_probs(mak_filter))/nansum(off_state_probs(mak_filter));
            kni_off_prob_array(k) = nansum(off_switch_probs(mak_filter))/nansum(on_state_probs(mak_filter));
        end
        if sum(~isnan(on_switch_flags(mak_filter)))>=25
            kni_on_array(k) = nanmean(on_switch_flags(mak_filter));
        end
        if sum(~isnan(off_switch_flags(mak_filter)))>=25
            kni_off_array(k) = nanmean(off_switch_flags(mak_filter));
        end
    end
    master_struct(m).kni_off_array = kni_off_array;
    master_struct(m).kni_on_array = kni_on_array;
    master_struct(m).kni_state_array = kni_state_array;
    
    master_struct(m).kni_on_prob_array = kni_on_prob_array;
    master_struct(m).kni_off_prob_array = kni_off_prob_array;
end    