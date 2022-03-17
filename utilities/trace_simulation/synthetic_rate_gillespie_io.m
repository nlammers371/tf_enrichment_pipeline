% Experiment with sochastic simulations in which rates vary in time
function temp_gillespie = synthetic_rate_gillespie_io(seq_length,r_emission,pi0,...
            t_MS2, memory, tf_profile, KD, hill_coeff, tf_dependent_flags,...
            RateMatrix, deltaT, noise, granularity, n_traces)

nStates = size(RateMatrix,1);          
% calculate time vector
t_process = deltaT*seq_length; % length of simulation in seconds
t_ref_in = (0:deltaT:t_process)';
t_ref_out = (0:granularity:t_process)'; % assuming that jump lengths are much shorter than obs time

% calculate I/O curve
io_curve_in = tf_profile.^hill_coeff ./ (KD.^hill_coeff +  tf_profile.^hill_coeff);
io_curve_out = interp1(t_ref_in,io_curve_in,t_ref_out);

%%% Make Rate Array (3D rate matrix)
R_array = repmat(RateMatrix,1,1,n_traces);
trend_flags = repmat(tf_dependent_flags,1,1,n_traces);
stable_flags = trend_flags~=1;
diag_flags = repmat(eye(nStates),1,1,n_traces);

% generate indexing vectors
trace_ind_vec = (0:n_traces-1)*nStates^2;

% Initialize state array
promoter_state_array = NaN(length(t_ref_out),n_traces);
promoter_state_array(1,:) = randsample(1:nStates,n_traces,true,pi0);

% iterate through time points
for t = 2:length(t_ref_out)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% calculate updated transition rate array
    % apply time trend        
    trend_array_temp = trend_flags .* reshape(io_curve_out(t,:),1,1,[]);
    trend_array_temp(stable_flags) = 1;
    R_array_temp = R_array.*trend_array_temp;

    % renormalize    
    R_array_temp(diag_flags==1) = 0;
    R_array_temp(diag_flags==1) = -sum(R_array_temp,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Choose next states
    
    % get previous state
    prev_state_vec = promoter_state_array(t-1,:);
    
    % get mean jump times 
    prev_lin_index_vec = prev_state_vec + nStates*(prev_state_vec-1) + trace_ind_vec;
    try
      tau_vec = -1./R_array_temp(prev_lin_index_vec);
    catch
      error('wtf')
    end
    % randomly select next jump time for each trace
    next_jump_times = exprnd(tau_vec);
    
    % randomly select states to jump to
    next_lin_index_vec = (1:nStates)' + nStates*(prev_state_vec-1) + trace_ind_vec;
    next_rates = R_array_temp(next_lin_index_vec);
    next_rates(next_rates<0) = 0;
    next_options = cumsum(next_rates ./ sum(next_rates));
    rnd_vec = rand(1,n_traces);
    next_states = sum(next_options<rnd_vec,1) + 1;
    
    % update states for cases when next jump is within sampling resolution
    accepted_jumps =next_jump_times<granularity;
    promoter_state_array(t,accepted_jumps) = next_states(accepted_jumps);
    promoter_state_array(t,~accepted_jumps) = promoter_state_array(t-1,~accepted_jumps);
            
end

% perform convolution to obtain predicted fluorescence
fluo_kernel = ms2_loading_coeff (t_MS2, memory);
fluo_kernel_full = interp1((0:memory-1)*deltaT,fluo_kernel,0:granularity:(memory-1)*deltaT);
initiation_state_array = r_emission(promoter_state_array);
fluo_ms2_array = convn(fluo_kernel_full',initiation_state_array,'full');
fluo_ms2_array = fluo_ms2_array(1:end-length(fluo_kernel_full)+1,:);

% add gaussian noise
fluo_ms2_array_noise = fluo_ms2_array + normrnd(0,noise,size(fluo_ms2_array));

% downsample
fluo_ms2_array_ds = interp1(t_ref_out,fluo_ms2_array,t_ref_in);
fluo_ms2_array_noise_ds = interp1(t_ref_out,fluo_ms2_array_noise,t_ref_in);

% % record params
temp_gillespie.fluo_ms2_array = fluo_ms2_array_ds;
temp_gillespie.fluo_ms2_array_noise = fluo_ms2_array_noise_ds;
temp_gillespie.fluo_ms2_array_full = fluo_ms2_array;
temp_gillespie.fluo_ms2_array_noise_full = fluo_ms2_array_noise;
temp_gillespie.promoter_state_array = promoter_state_array;
temp_gillespie.initiation_rate_array = initiation_state_array;
temp_gillespie.t_ref_full = t_ref_out;
temp_gillespie.t_ref = t_ref_in;
temp_gillespie.io_ref_in = io_curve_in;
temp_gillespie.io_ref_out = io_curve_out;
