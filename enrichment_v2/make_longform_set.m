% script to write key enrichment info to longform csv table
clear
close all

projectName = 'Bcd-GFP_hbMS2-mCh_Airy_fast';%'Bcd-GFP-McpMcherry-hbP2P-delta6';
liveProject = LiveEnrichmentProject(projectName);
FigurePath = liveProject.figurePath;
resultsRoot = [liveProject.dataPath filesep];

% load data
protein_flag = exist([resultsRoot 'spot_struct_protein.mat']);

load([resultsRoot 'spot_struct_protein.mat'])
load([resultsRoot 'spot_struct.mat'])
pt_id_vec_trunc = [spot_struct.particleID];

% Build stripped-down data structure for export
keep_fields = {'particleID','nc','xPosParticle','yPosParticle','zPosParticle','fluo',...
                  'edge_null_protein_vec','spot_protein_vec','nuclear_protein_vec'};
out_names = {'particleID','nuclearCycle','xPosParticle','yPosParticle','zPosParticle','spotFluo',...
                  'controlSpotProtein','spotProtein','nucleusProtein'};                
% extend_flags = false(size(keep_fields));
% extend_flags(1) = true;

% get list of fieldnames
fnames = fieldnames(spot_struct_protein);
keepFlags = ismember(fnames,keep_fields);

% generate smaller structure
spot_struct_trunc = rmfield(spot_struct_protein,fnames(~keepFlags));
for i = 1:length(spot_struct_trunc)
    ptID = spot_struct_trunc(i).particleID;
    nc = spot_struct(pt_id_vec_trunc==ptID).nc;
    spot_struct_trunc(i).particleID = repelem(ptID,length(spot_struct_trunc(i).fluo));
    spot_struct_trunc(i).nc = repelem(nc,length(spot_struct_trunc(i).fluo));
end

spot_array = [];
kept_flags = false(size(keep_fields));
kept_flags(1:2) = true;
iter = 1;
for k = 1:length(keep_fields)
    try
        spot_array(:,iter) = [spot_struct_trunc.(keep_fields{k})]';
        kept_flags(iter) = true;
        iter = iter + 1;
    catch
        % do nothing
    end
end
spot_table = array2table(spot_array(:,kept_flags),'VariableNames',out_names(kept_flags));
writetable(spot_table,[resultsRoot 'longform_particle_data.csv'])
