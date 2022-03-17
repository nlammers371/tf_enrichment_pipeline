clear
close all

addpath(genpath('../utilities'));

% define sim parameters
seq_length = 120;
memory = 7;
deltaT = 20;
t_MS2 = 1.4;
noise = 1; % noise level in 
n_traces = 10;
% specify network architecture
RateMatrix = [-1 60 0; 
              1 -62 1; 
              0 2 -1]/60;
r_emission = [0 0 1]; % loading rate for each state
pi0 = [0 0.5 0.5];

% dictate which rate is tf-dependent (assume only one possible for now)
tf_dependent_flags = false(size(RateMatrix));
tf_dependent_flags(1,2) = true;

% specify io response characteristics
KD = 10;
hill_coeff = 6;

% calculate time resolution that we need to use
granularity = 1./max(abs(RateMatrix(:))) / 10;

% generate hypothetical temporal TF input profile (assume linear)
m_tf = 0.25*KD / 5 * 60/ deltaT;
p1 = 1e-3;
p0 = 1-1e-3;
r_vec = [KD*(p0/(1-p0))^(1/hill_coeff) KD*(p1/(1-p1))^(1/hill_coeff)];
r0 = min(r_vec);
r1 = max(r_vec);
t_shift = 10*60/deltaT;
center_time = seq_length/2;
tf_profile = NaN(seq_length+1,1);
tf_profile(1:center_time-t_shift/2) = r0;
tf_profile(center_time+2+t_shift/2:end) = r1;
tf_profile(center_time+1-t_shift/2:center_time+1+t_shift/2) = linspace(r0,r1,t_shift+1);
tf_profile(tf_profile<0) = 1e-6;
%%
n_traces = 100;
tic
gillespie = synthetic_rate_gillespie_io(seq_length,r_emission,pi0,...
            t_MS2, memory, tf_profile, KD, hill_coeff, tf_dependent_flags,...
            RateMatrix, deltaT, noise, granularity, n_traces);
toc

% calculate last shut-off times
detection_thresh = 10;
on_off_array = ones(size(gillespie.fluo_ms2_array));
for n = 1:n_traces
    last_i = find(gillespie.fluo_ms2_array(:,n)>detection_thresh,1,'last');
    on_off_array(last_i+1:end,n) = 0;
end