% Experiment with sochastic simulations in which rates vary in time
function temp_gillespie = synthetic_rate_gillespie_io_v3(simInfo,tf_profile_array)

% extract key parameters
KD = simInfo.KD;
HC = simInfo.HC;
deltaT = simInfo.deltaT;
seq_length = size(tf_profile_array,1);
n_traces = size(tf_profile_array,3);
tf_dependent_flags = simInfo.tf_dependent_flags;          
granularity = simInfo.granularity;

% trace sim characteristics
noise = simInfo.noise;
t_MS2 = simInfo.t_MS2;
memory = simInfo.memory;
RateMatrix = simInfo.RateMatrix;
r_emission = simInfo.r_emission;
pi0 = simInfo.pi0;
nStates = size(RateMatrix,1);  

% calculate time vector
t_process = deltaT*(seq_length-1); % length of simulation in seconds
t_ref_in = (0:deltaT:t_process)';
t_ref_out = (0:granularity:t_process)'; % assuming that jump lengths are much shorter than obs time

% calculate I/O curve
io_curve_in = tf_profile_array.^HC ./ (KD.^HC +  tf_profile_array.^HC);
rate_curve_in = NaN(size(io_curve_in));
io_curve_out = interp1(t_ref_in,io_curve_in,t_ref_out);

%%% Make Rate Array (3D rate matrix)
R_array = repmat(RateMatrix,1,1,n_traces);
trend_flags = repmat(tf_dependent_flags,1,1,n_traces);
stable_flags = trend_flags~=1;
diag_flags = repmat(eye(nStates),1,1,n_traces);
offset_array = trend_flags*0;%simInfo.k0;

% generate indexing vectors
trace_ind_vec = (0:n_traces-1)*nStates^2;

% Initialize state array
promoter_state_array = NaN(length(t_ref_out),n_traces);
init_vec = randsample(1:nStates,n_traces,true,pi0);

% iterate through time points
for t = 1:length(t_ref_out)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% calculate updated transition rate array
    % apply time trend        
    trend_array_temp = trend_flags .* io_curve_out(t,1,:);
    
    trend_array_temp(stable_flags) = 1;
    R_array_temp = R_array.*trend_array_temp + offset_array;

    if ismember(t_ref_out(t),t_ref_in) && any(tf_dependent_flags(:))
      rate_curve_in(t_ref_in==t_ref_out(t),1,:) = R_array_temp(trend_flags==1);
    end
    % renormalize    
    R_array_temp(diag_flags==1) = 0;
    R_array_temp(diag_flags==1) = -sum(R_array_temp,1);  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Choose next states
    
    % get previous state
    if t > 1
        prev_state_vec = promoter_state_array(t-1,:);
    else
        prev_state_vec = init_vec;
    end
    % get mean jump times 
    prev_lin_index_vec = prev_state_vec + nStates*(prev_state_vec-1) + trace_ind_vec;    
    tau_vec = -1./R_array_temp(prev_lin_index_vec);
    
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
    accepted_jumps = next_jump_times<granularity;
    promoter_state_array(t,accepted_jumps) = next_states(accepted_jumps);
    if t > 1
        promoter_state_array(t,~accepted_jumps) = promoter_state_array(t-1,~accepted_jumps);
    else
        promoter_state_array(t,~accepted_jumps) = init_vec(~accepted_jumps);
    end    
end

% perform convolution to obtain predicted fluorescence
fluo_kernel = ms2_loading_coeff(t_MS2, memory);
fluo_kernel_full = interp1((0:memory-1)*deltaT,fluo_kernel,0:granularity:(memory-1)*deltaT);
initiation_state_array = r_emission(promoter_state_array) * length(fluo_kernel) / length(fluo_kernel_full);
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
temp_gillespie.rate_curve_in = rate_curve_in;
temp_gillespie.tf_ref_in = tf_profile_array;
temp_gillespie.io_ref_out = io_curve_out;
