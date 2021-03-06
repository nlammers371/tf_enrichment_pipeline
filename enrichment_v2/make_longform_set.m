% script to write key enrichment info to longform csv table
clear
close all

projectName = 'Bcd-GFP_hbMS2-mCh_AiryscanTest_';
liveProject = LiveEnrichmentProject(projectName);
FigurePath = liveProject.figurePath;
resultsRoot = [liveProject.dataPath filesep];

% load data
protein_flag = exist([resultsRoot 'spot_struct_protein.mat']);
if protein_flag
    load([resultsRoot 'spot_struct_protein.mat'])
else
    load([resultsRoot 'spot_struct.mat'])
end    

% Build stripped-down data structure for export
keep_fields = {'particleID','xPosParticle','yPosParticle','zPosParticle','fluo',...
                  'edge_null_protein_vec','spot_protein_vec','nuclear_protein_vec'};
out_names = {'particleID','xPosParticle','yPosParticle','zPosParticle','spotFluo',...
                  'controlSpotProtein','spotProtein','nucleusProtein'};                
extend_flags = false(size(keep_fields));
extend_flags(1) = true;

% get list of fieldnames
fnames = fieldnames(spot_struct_protein);
keepFlags = ismember(fnames,keep_fields);

% generate smaller structure
spot_struct_trunc = rmfield(spot_struct_protein,fnames(~keepFlags));
for i = 1:length(spot_struct_trunc)
    ptID = spot_struct_trunc(i).particleID;
    spot_struct_trunc(i).particleID = repelem(ptID,length(spot_struct_trunc(i).fluo));
end

spot_array = [];
for k = 1:length(keep_fields)
    spot_array(:,k) = [spot_struct_trunc.(keep_fields{k})]';
end
spot_table = array2table(spot_array,'VariableNames',out_names);
writetable(spot_table,[resultsRoot 'longform_particle_data.csv'])
