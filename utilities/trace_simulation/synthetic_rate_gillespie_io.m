% Experiment with sochastic simulations in which rates vary in time
function temp_gillespie = synthetic_rate_gillespie_io(seq_length,r_emission,pi0,...
            t_MS2, memory, tf_profile, KD, hill_coeff, tf_dependent_flags,...
            RateMatrix, deltaT, noise, granularity)

% calculate time vector
t_process = deltaT*seq_length; % length of simulation in seconds
t_ref_in = 0:deltaT:t_process;
t_ref_out = 0:granularity:t_process; % assuming that jump lengths are much shorter than obs time

% calculate I/O curve
io_curve_in = tf_profile.^hill_coeff ./ (KD.^hill_coeff +  tf_profile.^hill_coeff);
io_curve_out = interp1(t_ref_in,io_curve_in,t_ref_out);

%%% Make Rate Array (3D rate matrix)
R_array = repmat(RateMatrix,1,1,length(t_ref_out));

% apply timne trend
trend_array = repmat(tf_dependent_flags,1,1,length(t_ref_out));
stable_flags = trend_array~=1;
trend_array = trend_array .* reshape(io_curve_out,1,1,[]);
trend_array(stable_flags) = 1;
R_array = R_array.*trend_array;

% renormalize
diag_flags = repmat(eye(size(RateMatrix,1)),1,1,length(t_ref_out));
R_array(diag_flags==1) = 0;
R_array(diag_flags==1) = -sum(R_array,1);

%%% Simulations
jump_times = [0];
promoter_states = [randsample(1:3,1,true,pi0)];
obs_time_grid = deltaT:deltaT:t_process;
state_vec = 1:3;
for ts = 1:length(t_ref_out)
    T = t_ref_out(ts);
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
    t_start = max(0,t_end - memory*deltaT);
    s_ind = find(jump_times<=t_start,1,'last'); % index of first relevant promoter state
  
    e_ind = find(jump_times<t_end,1,'last'); % last relevant state
    times = jump_times(s_ind:e_ind);
    times(1) = t_start;
    times = [times t_end];
    states = promoter_states(s_ind:e_ind);
    F = 0;
    F_alpha = 0;
    for t = 1:(length(times)-1)
        r_state = r_emission(states(t));
  
        t1 = t_end - times(t);
        t2 = t_end - times(t+1);
        ta1 = min(t1,t_MS2);
        ta2 = min(t2,t_MS2);
        F = F + r_state*(t1 - t2);
        F_alpha = F_alpha + r_state*(ta1-ta2)*(ta1+ta2)/t_MS2*.5 + ...
            ((t1 - t2)-(ta1-ta2))*r_state;

    end    
    fluo_grid(window) = F;    
    fluo_grid_MS2(window) = F_alpha;
end
noise_vec = normrnd(0,noise,1,length(fluo_grid));
off_ind = fluo_grid==0;
fluo_grid_noise = fluo_grid + noise_vec;
fluo_grid_noise(fluo_grid_noise<0) = 0;

fluo_grid_MS2_noise = fluo_grid_MS2 + noise_vec;
fluo_grid_MS2_noise(fluo_grid_MS2_noise<0) = 0;

% record params
temp_gillespie.fluo_no_noise = fluo_grid;
temp_gillespie.fluo = fluo_grid_noise;
temp_gillespie.fluo_MS2_no_noise = fluo_grid_MS2;
temp_gillespie.fluo_MS2 = fluo_grid_MS2_noise;
temp_gillespie.naive_states = promoter_states;
temp_gillespie.jump_times = jump_times;
temp_gillespie.t_ref = t_ref_out;
temp_gillespie.k_on_ref = k_on_ref;
temp_gillespie.k_off_ref = k_off_ref;
temp_gillespie.r_ref = r_ref;