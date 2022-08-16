clear
close all

addpath(genpath('../utilities'));

w = 7; % number time steps to elongate gene
alpha = 1.4; % MS2 rise time
K = 2;
R2 = [-2 1; 2 -1]/60;
r_emission = [0 1]; % loading rate for each 
noise = 10; % noise level in gluoprescece
%noise = 0; % noise level in gluoprescece
pi0 =[0.5 0.5]; % initaial state PDF
seq_length = 100;
deltaT = 20; % time resolution

gillespie = synthetic_rate_gillespie(seq_length, alpha, ...
                                K, w, R2, deltaT, r_emission, noise, pi0);
                              
figure;
plot(gillespie.fluo_MS2)

%figure;
%stairs(gillespie.naive_states)
%ylim([0.5 2.5])

% temp_gillespie = synthetic_rate_gillespie_linear(seq_length,...
%             t_MS2, w, k_on_1,k_on_2,k_off_1,k_off_2, deltaT, r_1, r_2, noise, ...
%             pi0,t_shift,shift_width,granularity,z_flag)

k1 = 1/1.8/60;
RFull = R2;
RFull(:,3) = 0;
RFull(3,1:2) = k1;
RFull(1,1) = RFull(1,1)-k1;
RFull(2,2) = RFull(2,2)-k1;