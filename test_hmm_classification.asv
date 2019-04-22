% Script to test feasibility of using mHMM to discirminate between gene
% types
clear
close all

% set id variables and load data
dropboxFolder = 'E:\Nick\Dropbox (Garcia Lab)\ProcessedEnrichmentData\';
figPath = 'E:\Nick\Dropbox (Garcia Lab)\LocalEnrichmentFigures\';

project_cell = {'Dl_Venus_snaBAC_MCPmCherry_Leica_Zoom2_7uW14uW','Dl_Venus_hbP2P_MCPmCherry_Zoom2_7uw14uw'};
master_struct = struct;
for i = 1:numel(project_cell)
    master_struct(i).project = project_cell{i};
    load([dropboxFolder project_cell{i} '/hmm_input_output_w6_K3.mat'])
    master_struct(i).hmm_input_output = hmm_input_output;
    
    % get other project info
    project = project_cell{i};
    underscores = strfind(project,'_');
    master_struct(i).protein_name = project(1:underscores(1)-1);
    master_struct(i).protein_fluor = project(underscores(1)+1:underscores(2)-1);
    master_struct(i).gene_name = project(underscores(2)+1:underscores(3)-1);
    if numel(underscores) == 3
        ind = numel(project);
    else
        ind = underscores(4)-1;
    end
    master_struct(i).gene_fluor = project(underscores(3)+1:ind);
end

%% Attempt to use mHMM likelihood to differentiate between gene types
w = 6;
K = 2;
Tres = 20;
max_time = 20*60;
addpath('./utilities')
samp_size = 4000;
for i = 1:numel(master_struct)
    ids = find([master_struct(i).hmm_input_output.mcp_qc_flag]);
    fluo_data_full = cell(size(ids'));
    for j = ids
        fvec = master_struct(i).hmm_input_output(j).fluo_hmm';
        fvec = fvec(~isnan(fvec));
        tvec = master_struct(i).hmm_input_output(j).time;
        tvec = tvec(~isnan(fvec));
%         fvec = fvec(tvec<max_time);
        fluo_data_full{j} = fvec; 
    end
    sample_ids = randsample(1:numel(fluo_data_full),numel(fluo_data_full),false);
    ndp = 0;
    indices = [];
    iter = 1;
    s_size = min([samp_size,numel([fluo_data_full{:}])]);
    while ndp < s_size
        f = fluo_data_full{sample_ids(iter)};        
        if numel(f) > 10
            ndp = ndp + numel(f);
            indices = [indices sample_ids(iter)];
        end
        iter = iter + 1;
    end
    fluo_data = fluo_data_full(indices);
    tic
    output = simple_hmm_inference(fluo_data,w,K,Tres);
    toc
    master_struct(i).hmm_params = output;
    master_struct(i).fluo_data = fluo_data;
    master_struct(i).inference_ids = indices;
end
%%
%%% Now use inferred parameters to score loci...can we use this to discriminate between gene types?
for i = 1:numel(master_struct)
    % get fluo data
    fluo_data = master_struct(i).fluo_data;
    params_one = master_struct(1).hmm_params;
    params_two = master_struct(2).hmm_params;
    tic
    % get scores using first param set
    [master_struct(i).logL_tot_one, master_struct(i).logL_vec_one] = ...
        likelihood_reduced_memory_adapted(fluo_data, params_one.r'*Tres, params_one.noise, ...
                            log(params_one.pi0),log(params_one.A_mat), K, w, 1302 / 6444*w);
    % now second
    [master_struct(i).logL_tot_two, master_struct(i).logL_vec_two] = ...
        likelihood_reduced_memory_adapted(fluo_data, params_two.r'*Tres, params_two.noise, ...
                            log(params_two.pi0),log(params_two.A_mat), K, w, 1302 / 6444*w);                                 
    toc
end   

% make scatter plots of off vs on target likelihoods for each gene
logL_fig = figure;
hold on
gene_name_cell = {};
s = [];
for i = 1:numel(master_struct)
    logL_one_vec = master_struct(i).logL_vec_one;
    logL_two_vec = master_struct(i).logL_vec_two;
    s = [s histogram(logL_one_vec - logL_two_vec)];
    gene_name_cell = [gene_name_cell{:} {master_struct(i).gene_name}];
end
legend(s,gene_name_cell{:})
xlabel([master_struct(1).gene_name ' model score'])
ylabel([master_struct(2).gene_name ' model score'])

saveas(logL_fig,[figPath 'logL_scatter.png'])