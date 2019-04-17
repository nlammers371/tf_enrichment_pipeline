% Script to test feasibility of using mHMM to discirminate between gene
% types
clear
close all

% set id variables and load data
dropboxFolder = 'E:\Nick\Dropbox (Garcia Lab)\ProcessedEnrichmentData\';
figPath = 'E:\Nick\Dropbox (Garcia Lab)\LocalEnrichmentFigures\';

project_cell = {'Dl_Venus_snaBAC_MCPmCherry_Leica_Zoom2_7uW14uW','Bcd_GFP_hbP2P_MCPmCherry_Leica_6uW15uW'};
master_struct = struct;
for i = 1:numel(project_cell)
    master_struct(i).project = project_cell{i};
    load([dropboxFolder project_cell{i} '/nucleus_struct.mat'])
    master_struct(i).nucleus_struct = nucleus_struct;
    
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
%%%

% assign pseudo sister pairings from across two sets
mindp = 20; % only nuclei with 20 time-step temporal overlap will be permitted
% ugh need to interpolate...
interpGrid = 0:20:60*60;
interp_fields = {'fluo'};
for i = 1:numel(master_struct)
    nucleus_struct = master_struct(i).nucleus_struct;
    for j = 1:numel(nucleus_struct)
        nc_time = nucleus_struct(j).time;
        time_interp = interpGrid(interpGrid<=nc_time(end)&interpGrid>=nc_time(1));
        nucleus_struct(j).time_interp = time_interp;
        qc_flag = true;
        on_time = NaN;
        off_time = NaN;
        if numel(nc_time) > 1
            for k = 1:numel(interp_fields)
                vec = nucleus_struct(j).(interp_fields{k});                
                vec_ft = ~isnan(vec);
                qc_flag = sum(vec_ft) > mindp;
                if sum(vec_ft) > 1 
                    nucleus_struct(j).([interp_fields{k} '_interp']) = ...
                        interp1(nc_time(vec_ft),vec(vec_ft),time_interp);
                else
                    nucleus_struct(j).([interp_fields{k} '_interp']) = NaN(size(time_interp));
                end
                on_time = time_interp(find(~isnan(nucleus_struct(j).fluo_interp),1));
                off_time = time_interp(find(~isnan(nucleus_struct(j).fluo_interp),1,'last'));
            end
            nucleus_struct(j).time_interp = time_interp;
        else
            qc_flag = false;
            for k = 1:numel(interp_fields)
                vec = nucleus_struct(j).(interp_fields{k});                
                nucleus_struct(j).([interp_fields{k} '_interp']) = vec;
            end
            time_interp = nc_time;
            nucleus_struct(j).time_interp = time_interp;
        end
        if sum(vec_ft) > mindp && numel(nc_time) > 1 && ~qc_flag
            error('wtf')
        end
        nucleus_struct(j).qc_flag = qc_flag;
        nucleus_struct(j).first_time = on_time;
        nucleus_struct(j).last_time = off_time;
    end
    master_struct(i).nucleus_struct = nucleus_struct;    
end

%% Attempt to use mHMM likelihood to differentiate between gene types
w = 6;
K = 3;
Tres = 20;
max_time = 20*60;
addpath('./utilities')
samp_size = 4000;
for i = 1:numel(master_struct)
    ids = find([master_struct(i).nucleus_struct.qc_flag]);
    fluo_data_full = cell(size(ids'));
    for j = ids
        fvec = master_struct(i).nucleus_struct(j).fluo_interp;
        fvec = fvec(~isnan(fvec));
        tvec = master_struct(i).nucleus_struct(j).time_interp;
        tvec = tvec(~isnan(fvec));
        fvec = fvec(tvec<max_time);
        fluo_data_full{j} = fvec; 
    end
    sample_ids = randsample(1:numel(fluo_data_full),numel(fluo_data_full),false);
    ndp = 0;
    indices = [];
    iter = 1;
    while ndp < samp_size
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
    master_struct(i).inference_ids = indices;
end
%%
%%% Now use inferred parameters to score loci...can we use this to discriminate between gene types?
for i = 1:numel(master_struct)
    % get fluo data
    ids = find([master_struct(i).nucleus_struct.qc_flag]);
    fluo_data_full = cell(size(ids'));   
    for j = ids
        fvec = master_struct(i).nucleus_struct(j).fluo_interp;
        fvec = fvec(~isnan(fvec));
        tvec = master_struct(i).nucleus_struct(j).time_interp;
        tvec = tvec(~isnan(fvec));
        fvec = fvec(tvec<max_time);
        fluo_data_full{j} = fvec; 
    end
    fluo_data = fluo_data_full(master_struct(i).inference_ids);
    params_one = master_struct(1).hmm_params;
    params_two = master_struct(2).hmm_params;
    tic
    % get scores using first param set
    [master_struct(i).logL_tot_one, master_struct(i).logL_cell_one] = ...
        likelihood_reduced_memory_adapted(fluo_data, params_one.r'*Tres, params_one.noise, ...
                            log(params_one.pi0),log(params_one.A_mat), K, w, 1302 / 6444*w);
    % now second
    [master_struct(i).logL_tot_two, master_struct(i).logL_cell_two] = ...
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
    logL_one_vec = NaN(size(master_struct(i).logL_cell_on));
    logL_two_vec = NaN(size(master_struct(i).logL_cell_on));
    for j = 1:numel(logL_one_vec)
        logL_one_vec(j) = mean(master_struct(i).logL_cell_one{j});
        logL_two_vec(j) = mean(master_struct(i).logL_cell_two{j});
    end
    s = [s histogram(logL_one_vec./logL_two_vec)];
    gene_name_cell = [gene_name_cell{:} {master_struct(i).gene_name}];
end
legend(s,gene_name_cell{:})
xlabel([master_struct(1).gene_name ' model score'])
ylabel([master_struct(2).gene_name ' model score'])

saveas(logL_fig,[figPath 'logL_scatter.png'])