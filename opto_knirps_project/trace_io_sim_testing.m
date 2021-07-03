clear
close all

addpath(genpath('../utilities'));

% define sim parameters
seq_length = 120;
memory = 7;
deltaT = 20;
t_MS2 = 1.4;
noise = 1; % noise level in 

% specify network architecture
RateMatrix = [-1 1 0; 1 -3 1; 0 2 -1]/60;
r_emission = [0 0 1]; % loading rate for each state
pi0 = [0 0.5 0.5];

% dictate which rate is tf-dependent (assume only one possible for now)
tf_dependent_flags = false(size(RateMatrix));
tf_dependent_flags(1,2) = true;

% specify io response characteristics
KD = 5;
hill_coeff = -2;

% calculate time resolution that we need to use
granularity = 1./max(abs(RateMatrix(:))) / 10;

% generate hypothetical temporal TF input profile (assume linear)
m_tf = 0.25*KD / 5 * 60/ deltaT;
r0 = 1;
r1 = 11;
t_shift = 10*60/deltaT;
center_time = seq_length/2;
tf_profile = NaN(1,seq_length+1);
tf_profile(1:center_time-t_shift/2) = r0;
tf_profile(center_time+2+t_shift/2:end) = r1;
tf_profile(center_time+1-t_shift/2:center_time+1+t_shift/2) = linspace(r0,r1,t_shift+1);
% tf_profile = linspace(KD-seq_length/2*m_tf,KD+seq_length/2*m_tf,seq_length+1);
% tf_profile(tf_profile<1e-6) = 1e-6;

%%

gillespie = synthetic_rate_gillespie_io(seq_length,...
            t_MS2, memory, tf_profile, KD, hill_coeff, tf_dependent_flags,...
            RateMatrix, deltaT, noise, pi0, granularity);

