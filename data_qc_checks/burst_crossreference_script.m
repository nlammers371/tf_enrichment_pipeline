% Script to generate list of burst frames with apparent surge signal to
% cross-reference with raw data in CheckParticleTracking
clear 
close all
addpath('utilities')
% path to LivemRNA folder 
% livemRNAPath = 'C:\Users\nlamm\projects\LivemRNA';
currentDir - pwd;
livemRNAPath = 'P:\Nick\LivemRNA';
addpath(genpath([livemRNAPath '\mRNADynamics\']));
% path to processed data
% DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\';
DropboxFolder = 'S:\Nick\Dropbox\';
% set ID variables
project_cell = {'2xDl-Ven_snaBAC-mCh'};

% Params
fluo_dim = 2;
protein_dim = 2;
K = 3;
w = 7;

[~, ~, FigureRoot] =   header_function(DropboxFolder, project_cell{1}); 



% load data
master_struct = struct;
for i = 1:length(project_cell)
  % set write paths
  [~, DataPath] =   header_function(DropboxFolder, project_cell{i});   
  
  load([DataPath 'hmm_input_output_results_w' num2str(w) '_K' num2str(K) '_f' num2str(fluo_dim) 'D_p' num2str(protein_dim) 'D.mat'])
  master_struct(i).results_struct = results_struct;
  clear results_struct;
  
  load([DataPath 'nucleus_struct.mat'])
  master_struct(i).nucleus_struct = nucleus_struct;
  clear nucleus_struct;
end


Tres = 20; % seconds
% extract relevant arrays from target project 
for i = 1:length(master_struct)
  master_struct(i).lag_dur_vec = master_struct(i).results_struct.lag_dur_vec;
  master_struct(i).lead_dur_vec = master_struct(i).results_struct.lead_dur_vec;
  master_struct(i).hmm_array = master_struct(i).results_struct.hmm_array;   
  
  master_struct(i).feature_sign_vec = master_struct(i).results_struct.feature_sign_vec;   
end

% set basic analyisis parameters
nBoots = 100; % number of bootstrap samples to use
min_pause_len = 4; % minimum length of preceding OFF period (in time steps)
max_pause_len = 100;
min_burst_len = 2;
max_burst_len = 1000;

% initialize data structure to store key cross-reference info
cross_ref_struct = struct;

for i = 1:length(master_struct)
  % generate basic filter for target locus and computational controls
  master_struct(i).burst_ft = master_struct(i).feature_sign_vec == 1&master_struct(i).lead_dur_vec>=...
    min_pause_len&master_struct(i).lead_dur_vec<=max_pause_len...
      & master_struct(i).lag_dur_vec>=min_burst_len&master_struct(i).lag_dur_vec<=max_burst_len;%&target_swap_qc&target_virtual_qc;; % filter for rise events
    
  event_indices = find(master_struct(i).burst_ft);
  event_particles = master_struct(1).results_struct.particle_id_vec(event_indices);
  
  nucleus_struct = master_struct(i).nucleus_struct;
  particle_vec_nc = [nucleus_struct.ParticleID];
  
  cross_ref_struct(i).prefix_cell = cell(1,length(event_indices));
  cross_ref_struct(i).particles_index_vec = NaN(1,length(event_indices));
  cross_ref_struct(i).frame_vec = NaN(1,length(event_indices));
  
  for e = 1:length(event_indices)
    % get prefix
    source_path = nucleus_struct(particle_vec_nc==event_particles(e)).source_path;
    slashes = regexp(source_path,'[\\|/]');
    
    cross_ref_struct(i).prefix_cell{e} = source_path(slashes(end)+1:end);
    
    % get frame    
    event_time = master_struct(i).results_struct.center_time_vec(event_indices(e));
    [~, cross_ref_struct(i).frame_vec(e)] = min(abs(nucleus_struct(particle_vec_nc==event_particles(e)).time-event_time));
    
    % get corresponding index  in "Particles" structure
    ParticleID = event_particles(e);
    cross_ref_struct(i).particles_index_vec(e) = 1e4*(ParticleID-floor(ParticleID));
  end
end

% save data structure 
save([DataPath 'cross_ref_struct.mat'],'cross_ref_struct');

%% change directory so that we can use livemRNA functions
cd(livemRNAPath)

%% get unique list of Prefixes
for i = 1:length(cross_ref_struct)
  prefix_index = unique(cross_ref_struct(i).prefix_cell);  
  for p = 1:length(prefix_index)
    Prefix = prefix_index{p};
    liveExperiment = LiveExperiment(Prefix);
    % get indices that pertain to this prefix
    PrefixFilter = strcmp(cross_ref_struct(i).prefix_cell,Prefix);
    % make a structure
    surge_cross_ref_structure.frame_vec = cross_ref_struct(i).frame_vec(PrefixFilter);
    surge_cross_ref_structure.particles_index = cross_ref_struct(i).particles_index_vec(PrefixFilter);
    % save a reference structure
    save([liveExperiment.resultsFolder 'surge_cross_ref_structure.mat'],'surge_cross_ref_structure')
  end
end


%% lastly, change back
cd(origDir);
