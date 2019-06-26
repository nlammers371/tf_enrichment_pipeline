% Script to generate training set for event prediction

clear
close all

% define core ID variables
project = 'Dl-Ven_snaBAC-mCh';
% project = 'Dl-Ven_hbP2P-mCh';
dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
dataPath = [dropboxFolder 'ProcessedEnrichmentData\' project '/'];
% load input-output data set
K = 3;
w = 7;
load([dataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'],'hmm_input_output')
load([dataPath 'nucleus_struct_protein.mat'],'nucleus_struct_protein')

%%
% specify indexing vectors
particle_vec_hmm = [hmm_input_output.ParticleID];
% temporal proximity tolerance to burst stop or start event that
dT = 30; % seconds
% generate snip stacks
locus_protein_stack = cat(3,nucleus_struct_protein.spot_protein_snips);
locus_mcp_stack = cat(3,nucleus_struct_protein.spot_mcp_snips);
control_protein_stack = cat(3,nucleus_struct_protein.edge_null_protein_snips);
% initialize class vectors
event_class_vec = NaN(size(locus_protein_stack,3),1);
burst_class_vec = NaN(size(locus_protein_stack,3),1);
mf_protein_vec =  NaN(size(locus_protein_stack,3),1);
% generate indexing vectors
particle_vec_protein = NaN(size(event_class_vec));
time_vec_protein = [nucleus_struct_protein.time];
iter = 1;
for i = 1:numel(nucleus_struct_protein)
    fluo = nucleus_struct_protein(i).fluo;
    pID = nucleus_struct_protein(i).ParticleID;
    particle_vec_protein(iter:iter+numel(fluo)-1) = pID;
    iter = iter + numel(fluo);
end

% iterate through input/output data set and document events
for i = 1:numel(hmm_input_output)
    % extract HMM data
    change_points = hmm_input_output(i).change_points;
    z_diff_vec = hmm_input_output(i).z_diff_vec;
    z_vec = hmm_input_output(i).z_vec;
    time = hmm_input_output(i).time;
    mf_protein = hmm_input_output(i).mf_protein;    
    ptID = hmm_input_output(i).ParticleID;
    % extract raw protein time data
    pt_ids = find(particle_vec_protein==ptID)';
    for id = pt_ids
        pt_time = time_vec_protein(id);
        [minTime, mi] = min(abs(pt_time-time));
        if minTime <= dT
            event_class_vec(id) = z_diff_vec(mi);
            burst_class_vec(id) = z_vec(mi);
            mf_protein_vec(id) = mf_protein(mi);
        end
    end
end
    
% Finally, remove entries for which either class or data are Missing
nan_pt_stack_flag = reshape(sum(sum(isnan(locus_protein_stack),1),2),1,[]) < 100;
nan_class_vec = isnan(event_class_vec);
nan_ft = ~nan_class_vec' & ~nan_pt_stack_flag;
% filter sets
event_class_vec = event_class_vec(nan_ft);
burst_class_vec = burst_class_vec(nan_ft);
mf_protein_vec = mf_protein_vec(nan_ft);
time_vec_protein = time_vec_protein(nan_ft);
particle_vec_protein = particle_vec_protein(nan_ft);
locus_protein_stack = locus_protein_stack(:,:,nan_ft);
control_protein_stack = control_protein_stack(:,:,nan_ft);
locus_mcp_stack = locus_mcp_stack(:,:,nan_ft);
% set remaining NaNs in snips to -1
locus_protein_stack(isnan(locus_protein_stack)) = -1;
control_protein_stack(isnan(control_protein_stack)) = -1;
locus_mcp_stack(isnan(locus_mcp_stack)) = -1;
% store in structure
training_struct = struct;
training_struct.event_class_vec = event_class_vec;
training_struct.burst_class_vec = burst_class_vec;
training_struct.mf_protein_vec = mf_protein_vec;
training_struct.particle_vec_protein = particle_vec_protein;
training_struct.time_vec = time_vec_protein;
training_struct.locus_protein_stack = locus_protein_stack;
training_struct.control_protein_stack = control_protein_stack;
training_struct.locus_mcp_stack = locus_mcp_stack;
% save structure
save([dataPath 'event_training_data.mat'],'training_struct','-v7.3') 

