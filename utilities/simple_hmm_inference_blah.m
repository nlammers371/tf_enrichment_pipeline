% stripped-down inference script

function output = simple_hmm_inference(fluo_data,w,K,Tres,varargin)
% set defaults
MaxWorkers = 12;
alpha_frac = 1302 / 6444;
n_steps_max = 500;
eps = 1e-4;
n_localEM = 20;

for i = 1:numel(varargin)
    if ischar(varargin{i})
        eval([varargin{i} ' = varargin{i+1};'])
    end
end

alpha = alpha_frac*w;
% Random initialization of model parameters
param_init = initialize_random (K, w, fluo_data);
% Approximate inference assuming iid data for param initialization                
local_iid_out = local_em_iid_reduced_memory_truncated (fluo_data, param_init.v, ...
                    param_init.noise, K, w, alpha, n_steps_max, eps);
noise_iid = 1/sqrt(exp(local_iid_out.lambda_log));
v_iid = exp(local_iid_out.v_logs);            
p = gcp('nocreate');
if isempty(p)
    parpool(MaxWorkers); %12 is the number of cores the Garcia lab server can reasonably handle per user.
elseif p.NumWorkers > MaxWorkers
    delete(gcp('nocreate')); % if pool with too many workers, delete and restart
    parpool(MaxWorkers);
end
parfor i_local = 1:n_localEM % Parallel Local EM                
    % Random initialization of model parameters
    param_init = initialize_random_with_priors(K, noise_iid, v_iid);
    % Get Intial Values
    pi0_log_init = log(param_init.pi0);
    A_log_init = log(param_init.A);
    v_init = param_init.v;                        
    noise_init = param_init.noise; 
    %--------------------LocalEM Call-------------------------%
    local_out = local_em_MS2_reduced_memory_truncated(fluo_data, ...
        v_init, noise_init, pi0_log_init', A_log_init, K, w, ...
        alpha, n_steps_max, eps);                    
    %---------------------------------------------------------%                
    % Save Results 
    local_struct(i_local).subset_id = i_local;
    local_struct(i_local).logL = local_out.logL;                
    local_struct(i_local).A = exp(local_out.A_log);
    local_struct(i_local).v = exp(local_out.v_logs).*local_out.v_signs;
    local_struct(i_local).r = exp(local_out.v_logs).*local_out.v_signs / Tres;                                
    local_struct(i_local).noise = 1/exp(local_out.lambda_log);
    local_struct(i_local).pi0 = exp(local_out.pi0_log);
    local_struct(i_local).soft_struct = local_out.soft_struct; 
end
[~, max_index] = max([local_struct.logL]); % Get index of best result                    
% Save parameters from most likely local run
output.pi0 =local_struct(max_index).pi0;                        
output.r = local_struct(max_index).r(:);           
output.noise = local_struct(max_index).noise;
output.A = local_struct(max_index).A(:);
output.A_mat = local_struct(max_index).A;            
% get soft-decoded structure
output.soft_struct = local_struct(max_index).soft_struct;           
output.w = w;
output.K = K;
output.alpha = alpha;
output.deltaT = Tres; 
