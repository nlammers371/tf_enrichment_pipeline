% Experiment with sochastic simulations in which rates vary in time
%% initialization

clear
close all

addpath(genpath('../utilities'));

%% set simulation parameters
w = 7; % number time steps to elongate gene
alpha = 1.4; % MS2 rise time
K = 2;
r_emission = [0 1]; % loading rate for each 
noise = 10; % noise level in gluoprescece
%noise = 0; % noise level in gluoprescece
pi0 =[1/3 1/3 1/3]; % initaial state PDF
seq_length = 100;
deltaT = 20; % time resolution
z_flag = 0;
granularity = 1;
t_MS2 = 20;
t_process = deltaT*seq_length; % length of simulation in seconds
% arrival_rate = r_emission(2)./f_per_mRNA;

t_ref = 0:granularity:t_process; % assuming that jump lengths are much shorter than obs time

r_ref = 1*ones(size(t_ref));
k_on_ref = 2/60*ones(size(t_ref));
k_off_ref = 1/60*ones(size(t_ref));

%k_on_1 = 1;
%k_on_2 = 1;
%k_off_1 = 2;
%k_off_2 = 2;
%r_1 = 1;
%r_2 = 2; 
%t_shift = 1000;
%shift_width = 2000;

%shift_start = t_shift - .5*shift_width;
% shift_stop = t_shift + .5*shift_width;
% if shift_start < 0 || shift_stop > t_process
%     error('incompatible shift width and midpoint parameters')
% end
% n_prior = sum(t_ref<shift_start);
% n_shift = sum(t_ref>=shift_start&t_ref<shift_stop);
% n_post = sum(t_ref>=shift_stop);
% k_on_ref = [repelem(k_on_1, n_prior)...
%             linspace(k_on_1,k_on_2,n_shift) repelem(k_on_2, n_post)];
% k_off_ref = [repelem(k_off_1,n_prior) ...
%             linspace(k_off_1,k_off_2,n_shift) repelem(k_off_2,n_post)];
% r_ref = [repelem(r_1,n_prior) ...
%             linspace(r_1,r_2,n_shift) repelem(r_2,n_post)];


%%% Make Initiation rate Array
r_array = [zeros(length(r_ref),1) r_ref' 2*r_ref'];
%%% Make Rate Array (3D rate matrix)
R_array = zeros(3,3,length(t_ref));
% from 0
R_array(1,1,:) = -2*k_on_ref;
R_array(2,1,:) = 2*k_on_ref;
R_array(3,1,:) = zeros(1,length(t_ref));
% from 1
R_array(1,2,:) = k_off_ref;
R_array(2,2,:) = -(k_on_ref+k_off_ref);
R_array(3,2,:) = k_on_ref;
% from 3
R_array(1,3,:) = zeros(1,length(t_ref));
R_array(2,3,:) = 2*k_off_ref;
R_array(3,3,:) = -2*k_off_ref;

%% Simulations
jump_times = [0];
promoter_states = [randsample(1:3,1,true,pi0)];
obs_time_grid = deltaT:deltaT:t_process;
state_vec = 1:3;
for ts = 1:length(t_ref)
    T = t_ref(ts);
    t_step = 0; % process time within discrete step
    R = R_array(:,:,ts); % extract rate matrix for time step    
    cs = promoter_states(end); % get current state
    options = state_vec(state_vec~=cs);
    while t_step < granularity
        tau = -1/R(cs,cs);
        dt = exprnd(tau); % select jump time
        t_step = t_step + dt;
        options = state_vec(state_vec~=cs);        
        % if jump time is within discrete window, record
        if t_step < granularity 
            jump_times = [jump_times T + t_step];
            cs = randsample(options,1,true,R(options(:),cs)');
            promoter_states = [promoter_states cs];            
        end
    end
end



%%% Generate trace from promoter trajectory
fluo_grid = zeros(1,length(obs_time_grid));
fluo_grid_MS2 = zeros(1,length(obs_time_grid));
for window = 1:length(fluo_grid)
    t_end = obs_time_grid(window);  
    t_start = max(0,t_end - w*deltaT);
    s_ind = find(jump_times<=t_start,1,'last'); % index of first relevant promoter state
    if isempty(s_ind)
        error('wtf')
    end
    e_ind = find(jump_times<t_end,1,'last'); % last relevant state
    times = jump_times(s_ind:e_ind);
    times(1) = t_start;
    times = [times t_end];
    states = promoter_states(s_ind:e_ind);
    F = 0;
    F_alpha = 0;
    for t = 1:(length(times)-1)
        r_mean = mean(r_array(t_ref>=times(t)&t_ref<times(t+1),states(t)));
        if isnan(r_mean)
            r_mean = r_array(find(t_ref>=times(t),1));
        end
        t1 = t_end - times(t);
        t2 = t_end - times(t+1);
        ta1 = min(t1,t_MS2);
        ta2 = min(t2,t_MS2);
        F = F + r_mean*(t1 - t2);
        F_alpha = F_alpha + r_mean*(ta1-ta2)*(ta1+ta2)/t_MS2*.5 + ...
            ((t1 - t2)-(ta1-ta2))*r_mean;
        if isnan(F)
            error('afsa')
        end
    end    
    fluo_grid(window) = F;    
    fluo_grid_MS2(window) = F_alpha;
end
noise_vec = normrnd(0,noise,1,length(fluo_grid));
off_ind = fluo_grid==0;
fluo_grid_noise = fluo_grid + noise_vec;
fluo_grid_noise(fluo_grid_noise<0) = 0;
if z_flag % imitate zero tails in real data
    fluo_grid_noise(off_ind) = 0;
end
fluo_grid_MS2_noise = fluo_grid_MS2 + noise_vec;
fluo_grid_MS2_noise(fluo_grid_MS2_noise<0) = 0;
if z_flag % imitate zero tails in real data
    fluo_grid_MS2_noise(off_ind) = 0;
end
% record params
temp_gillespie.fluo_no_noise = fluo_grid;
temp_gillespie.fluo = fluo_grid_noise;
temp_gillespie.fluo_MS2_no_noise = fluo_grid_MS2;
temp_gillespie.fluo_MS2 = fluo_grid_MS2_noise;
temp_gillespie.naive_states = promoter_states;
temp_gillespie.jump_times = jump_times;
temp_gillespie.t_ref = t_ref;
temp_gillespie.k_on_ref = k_on_ref;
temp_gillespie.k_off_ref = k_off_ref;
temp_gillespie.r_ref = r_ref;

%% plot results
fig = figure;
plot(temp_gillespie.fluo_MS2)

%fig = figure;
%plot(k_on_ref)
%hold on
%plot(k_off_ref)
%plot(r_ref)