% Script to assign likelihood scores to sna and hbP2P traces
function nucleus_struct = classify01_assign_spot_likelihood_scores(project,varargin)
% assumes that two possibilities are snail and hbP2P
sna_project = 'Dl-Ven x snaBAC';
hbP2P_project = 'Dl-Ven x hbP2P';

dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
dataRoot = [dropboxFolder 'ProcessedEnrichmentData\'];
w = 7;
K = 3;

for i = 1:(numel(varargin)-1)  
    if ischar(varargin{i}) && ~strcmpi(varargin{i},'dropboxFolder')
        if i ~= numel(varargin)
            if ~ischar(varargin{i+1})
                eval([varargin{i} '=varargin{i+1};']);        
            end
        end
    elseif strcmpi(varargin{i},'dropboxFolder')
        dataRoot = [varargin{i+1} '\ProcessedEnrichmentData\' project '/'];
    end
end
% load data set
load([dataRoot project '/nucleus_struct.mat'],'nucleus_struct') 
% useful params and vectors
qc_indices = find([nucleus_struct.qc_flag]==1);
Tres = nucleus_struct(1).TresInterp;

% extract fluorescence traces
fluo_data = cell(size(qc_indices));
iter = 1;
for  i = qc_indices
    fvec = nucleus_struct(i).fluo_interp;
    fvec = fvec(~isnan(fvec));
%         fvec = fvec(tvec<max_time);
    fluo_data{iter} = fvec; 
    iter = iter + 1;
end

% load snail reference model first
hmm_suffix =  ['hmm_inference/w' num2str(w) '_K' num2str(K) '/']; 
file_list = dir([dataRoot sna_project '/' hmm_suffix 'hmm_results*.mat']); 
if numel(file_list) > 1
    warning('multiple inference files detected. Ignoring all but first')
end
sna_results = load([dataRoot sna_project '/' hmm_suffix file_list(1).name]);
sna_results = sna_results.output;
disp('performing snail model classification...')
tic
[sna_logL_tot, sna_logL_vec_sub] = ...
        likelihood_reduced_memory_adapted(fluo_data, sna_results.r'*Tres, sqrt(sna_results.noise), ...
                            log(sna_results.pi0),log(sna_results.A_mat), K, w, 1302 / 6444*w);
disp('done.')                        
toc                     
% record results
sna_logL_vec_tot = NaN(size(nucleus_struct));
sna_logL_vec_tot(qc_indices) = sna_logL_vec_sub;
for i = 1:numel(nucleus_struct)
    nucleus_struct(i).sna_logL = sna_logL_vec_tot(i);
    nucleus_struct(i).sna_logL_overall = sna_logL_tot;
end

% load hbP2P reference model 
hmm_suffix =  ['hmm_inference/w' num2str(w) '_K' num2str(K) '/']; 
file_list = dir([dataRoot hbP2P_project '/' hmm_suffix 'hmm_results*.mat']); 
if numel(file_list) > 1
    warning('multiple inference files detected. Ignoring all but first')
end
hbP2P_results = load([dataRoot hbP2P_project '/' hmm_suffix file_list(1).name]);
hbP2P_results = hbP2P_results.output;
disp('performing hbP2P model classification...')
tic
[hbP2P_logL_tot, hbP2P_logL_vec_sub] = ...
        likelihood_reduced_memory_adapted(fluo_data, hbP2P_results.r'*Tres, sqrt(hbP2P_results.noise), ...
                            log(hbP2P_results.pi0),log(hbP2P_results.A_mat), K, w, 1302 / 6444*w);
disp('done.') 
toc

% record results
hbP2P_logL_vec_tot = NaN(size(nucleus_struct));
hbP2P_logL_vec_tot(qc_indices) = hbP2P_logL_vec_sub;
for i = 1:numel(nucleus_struct)
    nucleus_struct(i).hbP2P_logL = hbP2P_logL_vec_tot(i);
    nucleus_struct(i).hbP2P_logL_overall = hbP2P_logL_tot;
end

% save structure
save([dataRoot project '/nucleus_struct.mat'],'nucleus_struct')