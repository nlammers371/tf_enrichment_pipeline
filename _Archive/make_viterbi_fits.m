% Script to generate Viterbi Fits for Inference traces
close all
clear 
addpath('E:\Nick\projects\hmmm\src\utilities\');
% addpath('D:\Data\Nick\projects\hmmm\src\utilities\');
%------------------------------Set System Params--------------------------%
alpha_ref = 1.4; % MS2 rise time in time steps
fluo_type = 1; % type of spot integration used
clipped = 1; % if 0, traces are taken to be full length of nc14
clipped_ends = 1;
dynamic_bins = 1; % if 1, use time-resolved region classifications
t_window = 50; % size of sliding window used
t_inf = 25;
%-----------------------------ID Variables--------------------------------%
w_ref = 7; %memory assumed for inference
K = 3; %states used for final inference
Tres_ref = 20; %time resolution
Tres_active = 20;

w_active = round(Tres_ref/Tres_active*w_ref);
alpha_active = round(Tres_ref/Tres_active*alpha_ref,1);
aggregate_fits = 1;  % if 1 apply same params to each trace regardless of stripe identity
%%% id variables
datatype = 'weka';
inference_type = 'dp';
ref_project = 'eve7stripes_inf_2018_04_28'; %use hmm params from wt data
active_project = 'kr_eve2_reporter';
%%% Generate filenames and writepath
id_string = [ '/w' num2str(w_ref) '_t' num2str(Tres_ref) '_alpha' num2str(round(alpha_ref*10)) ...
    '_f' num2str(fluo_type) '_cl' num2str(clipped) '_no_ends' num2str(clipped_ends) ...
    '_tbins' num2str(dynamic_bins) '/K' num2str(K) '_t_window' num2str(t_window) ...
     '_t_inf' num2str(round(t_inf)) '_' inference_type '/']; 

InfResultPath = ['../../all_the_stripes_cd/dat/' ref_project id_string];
OutPath = ['../dat/' active_project id_string '/'];
mkdir(OutPath);
%%% Load Data
load([InfResultPath '\hmm_results_t_window' num2str(t_window) '_t_inf' num2str(t_inf) ...
    '.mat'])
load(['../dat/' active_project '/inference_nuclei_' active_project '_dT' num2str(Tres_active) '.mat'])

%%% Viterbi Fits
hmm_regions = [hmm_results.binID];
alpha_ref = hmm_results(1).alpha;
dT = hmm_results(1).dT;
viterbi_fit_struct = struct;
skipped_stripes = [];
i_pass = 1;

parfor i = 1:numel(nucleus_struct_final)
%     MeanAP = round(nucleus_struct_final(i).MeanAP);            
    if aggregate_fits
        hmm_bin = hmm_results(round(hmm_regions,1)==0);
    else
%         hmm_index = find(MeanAP<=hmm_rb&MeanAP>=hmm_lb;        
        hmm_bin = hmm_results(hmm_regions==stripe_id);
    end        
    viterbi_fit_struct(i).skipped = 0;    
    if isempty(hmm_bin)|| length(nucleus_struct_final(i).fluo_interp) < w_active
%         error('wtf')
        viterbi_fit_struct(i).skipped = 1;
        viterbi_fit_struct(i).ParticleID = nucleus_struct_final(i).ParticleID;
        continue
    end    
    v = hmm_bin.initiation_mean/60*Tres_active;
    noise = hmm_bin.noise_mean;
    pi0_log = log(ones(1,K)/3);
    A_mat = reshape(hmm_bin.A_mean,K,K);                
    A_log = log(A_mat^(Tres_active/Tres_ref));
    fluo = nucleus_struct_final(i).fluo_interp;            
    v_fit = viterbi (fluo, v', noise, pi0_log, A_log, K, w_active, alpha_active);            
    
    v_fit.time_exp = nucleus_struct_final(i).time_interp;
    v_fit.fluo_exp = fluo;            
    v_fit.v = v;
    v_fit.w = w_ref;            
    v_fit.alpha = alpha_ref;
    v_fit.noise = noise;
    v_fit.pi0 = exp(pi0_log);
    v_fit.A = exp(A_log);
    v_fit.agg_fit_flag = aggregate_fits;    
%     v_fit.trace_source = datapath;
    viterbi_fit_struct(i).v_fit = v_fit;
    viterbi_fit_struct(i).ParticleID = nucleus_struct_final(i).ParticleID;
%     i_pass = i_pass + 1;
    disp([num2str(i) ' of ' num2str(length(nucleus_struct_final)) ' completed'])
end
save([OutPath '/viterbi_fits_w' num2str(w_ref) '_t_window' num2str(t_window) '_t_inf' num2str(t_inf) ...
    '.mat'],'viterbi_fit_struct') 