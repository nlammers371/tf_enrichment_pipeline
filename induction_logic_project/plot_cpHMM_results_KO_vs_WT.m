% script to compile results from cpHMM inference
clear
close all
addpath(genpath('utilities'))

% projectNameCell = {'EveGtSL','EveGtSL-S1Null','EveWt','EveS1Null'};%};
zeroNameCell = {'20220912_KO_experiments_Oct4','20220912_KO_experiments_Nanog'};
low_order_flag = 0;
iter = 1;
% resultsRoot = 'S:\Nick\Dropbox\InductionLogic\';
zero_struct = struct;
for p = 1:length(zeroNameCell)
    
    % set project to analyze 
    projectName = zeroNameCell{p};

    % get path to results
    try
        liveProject = LiveEnrichmentProject(projectName);
        resultsDir = [liveProject.dataPath 'cpHMM_results' filesep];
    catch
%         resultsRoot = 'S:/Nick/Dropbox/ProcessedEnrichmentData/';
        if isdir('C:\Users\nlamm\Dropbox (Personal)\')
            resultsRoot = 'C:\Users\nlamm\Dropbox (Personal)\ProcessedEnrichmentData';
        elseif isdir([filesep 'Users' filesep 'nick' filesep 'Dropbox (Personal)' filesep 'ProcessedEnrichmentData' filesep])
            resultsRoot = [filesep 'Users' filesep 'nick' filesep 'Dropbox (Personal)' filesep 'ProcessedEnrichmentData' filesep];
        else
           resultsRoot = 'S:\Nick\Dropbox (Personal)\ProcessedEnrichmentData';
        end
        resultsDir = [resultsRoot filesep zeroNameCell{p} filesep 'cpHMM_results' filesep];
    end
    
    % get list of all inference subdirectories. By default, we'll generate
    % summaries for all non-empty inference sub-directory
    infDirList = dir([resultsDir 'w*']);

    % iterate through the directories and compile the results
    for inf = 1:length(infDirList)
        load([infDirList(inf).folder filesep 'compiledResults_' infDirList(inf).name '.mat'],'compiledResults');
        fnames = fieldnames(compiledResults);
        for f = 1:length(fnames)
            zero_struct(iter).(fnames{f}) = compiledResults.(fnames{f});
        end
        iter = iter + 1;
    end
end

% Now load opto results

% projectNameCell = {'EveGtSL','EveGtSL-S1Null','EveWt','EveS1Null'};%};
nzNameCell = {'20220912_KO_experiments_nz_Oct4','20220912_KO_experiments_nz_Nanog'};

iter = 1;
% resultsRoot = 'S:\Nick\Dropbox\InductionLogic\';
nz_struct = struct;
for p = 1:length(nzNameCell)
    
    % set project to analyze 
    projectName = nzNameCell{p};

    % get path to results
    try
        liveProject = LiveEnrichmentProject(projectName);
        resultsDir = [liveProject.dataPath 'cpHMM_results' filesep];
    catch
%         resultsRoot = 'S:/Nick/Dropbox/ProcessedEnrichmentData/';
        if isdir('C:\Users\nlamm\Dropbox (Personal)\')
            resultsRoot = 'C:\Users\nlamm\Dropbox (Personal)\ProcessedEnrichmentData';
        elseif isdir([filesep 'Users' filesep 'nick' filesep 'Dropbox (Personal)' filesep 'ProcessedEnrichmentData' filesep])
            resultsRoot = [filesep 'Users' filesep 'nick' filesep 'Dropbox (Personal)' filesep 'ProcessedEnrichmentData' filesep];
        else
           resultsRoot = 'S:\Nick\Dropbox (Personal)\ProcessedEnrichmentData';
        end
        resultsDir = [resultsRoot filesep nzNameCell{p} filesep 'cpHMM_results' filesep];
    end
    
    % get list of all inference subdirectories. By default, we'll generate
    % summaries for all non-empty inference sub-directory
    infDirList = dir([resultsDir 'w*']);

    % iterate through the directories and compile the results
    for inf = 1:length(infDirList)
        load([infDirList(inf).folder filesep 'compiledResults_' infDirList(inf).name '.mat'],'compiledResults');
        fnames = fieldnames(compiledResults);
        for f = 1:length(fnames)
            nz_struct(iter).(fnames{f}) = compiledResults.(fnames{f});
        end
        iter = iter + 1;
    end
end

%%%%%%%%%%%%%%%
% Make figures
fig_path = [resultsRoot filesep 'KO_vs_WT' filesep];
mkdir(fig_path);

cmap = brewermap(10,'Set2');
close all

kon_vec_mean_z = [zero_struct.freq_vec_mean];
kon_vec_ste_z = [zero_struct.freq_vec_ste];

dur_vec_mean_z = [zero_struct.dur_vec_mean];
dur_vec_ste_z = [zero_struct.dur_vec_ste];

r_vec_mean_z = [zero_struct.init_vec_mean];
r_vec_ste_z = [zero_struct.init_vec_ste];

kon_vec_mean_nz = [nz_struct.freq_vec_mean];
kon_vec_ste_nz = [nz_struct.freq_vec_ste];

dur_vec_mean_nz = [nz_struct.dur_vec_mean];
dur_vec_ste_nz = [nz_struct.dur_vec_ste];

r_vec_mean_nz = [nz_struct.init_vec_mean];
r_vec_ste_nz = [nz_struct.init_vec_ste];

oct4_fig = figure;

hold on
scatter(1,r_vec_mean_z(2)/r_vec_mean_z(1),100,'^','MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k')
scatter(2,dur_vec_mean_z(2)/dur_vec_mean_z(1),100,'s','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k')
scatter(3,kon_vec_mean_z(2)/kon_vec_mean_z(1),100,'o','MarkerFaceColor',cmap(1,:),'MarkerEdgeColor','k')

set(gca,'yscale','log')
grid on
legend('amplitude','duration','frequency')
ylabel('KO/WT (oct4)')
set(gca,'xtick',1:3)
xlim([0.5 3.5])
ylim([1e-1 1e2])
saveas(oct4_fig,[fig_path 'wt_vs_ko_oct4.png'])
saveas(oct4_fig,[fig_path 'wt_vs_ko_oct4.pdf'])


nanog_fig = figure;

hold on
scatter(1,r_vec_mean_z(4)/r_vec_mean_z(3),100,'^','MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k')
scatter(2,dur_vec_mean_z(4)/dur_vec_mean_z(3),100,'s','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k')
scatter(3,kon_vec_mean_z(4)/kon_vec_mean_z(3),100,'o','MarkerFaceColor',cmap(1,:),'MarkerEdgeColor','k')

set(gca,'yscale','log')
grid on
ylabel('KO/WT (nanog)')
set(gca,'xtick',1:3)
xlim([0.5 3.5])
legend('amplitude','duration','frequency')
ylim([5e-1 5e0])
saveas(nanog_fig,[fig_path 'wt_vs_ko_nanog.png'])
saveas(nanog_fig,[fig_path 'wt_vs_ko_nanog.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%
%  calculate p-values

rng(13); % seed random number generator for consistency
n_boots = 1e4;

%%%%%%%%%%%%%%%%%%%%%%%%%
% nanog WT vs. KO
%%%%%%%%%%%%%%%%%%%%%%%%%
p_vec_kon = NaN(1,2);
p_vec_r = NaN(1,2);
p_vec_dur = NaN(1,2);

hyp_vec_kon = NaN(1,2);
hyp_vec_r = NaN(1,2);
hyp_vec_dur = NaN(1,2);

for i = 1:2

    % initialize reference vectors
    wt_outlier_flags = zero_struct(i).outlier_flags{1};
    wt_kon_vec = zero_struct(i).freq_results{1};
    wt_dur_vec = zero_struct(i).dur_results{1};
    wt_r_vec = zero_struct(i).r_results{1}(2,:);
    wt_options = find(~wt_outlier_flags);
    
    
    % extract key vectors
    ko_outlier_flags = zero_struct(i).outlier_flags{2};
    ko_kon_vec = zero_struct(i).freq_results{2};
    ko_dur_vec = zero_struct(i).dur_results{2};
    ko_r_vec = zero_struct(i).r_results{2}(2,:);
    ko_options = find(~ko_outlier_flags);
    
    % calculate value of mean fold qauntities
    fv = zero_struct(i).freq_vec_mean(2)/zero_struct(i).freq_vec_mean(1);
    dv = zero_struct(i).dur_vec_mean(2)/zero_struct(i).dur_vec_mean(1);
    rv = zero_struct(i).init_vec_mean(2)/zero_struct(i).init_vec_mean(1);
    
    % initialize arrays
    kon_boot_vec = NaN(1,n_boots);
    dur_boot_vec = NaN(1,n_boots);
    r_boot_vec = NaN(1,n_boots);
    
    % perform bootsrapping
    n_samp = min([length(ko_options) length(wt_options)]);
    for n = 1:n_boots
        % draw samples
        wt_ids_boot = randsample(wt_options,n_samp,true);
        ko_ids_boot = randsample(ko_options,n_samp,true);
        % calculate bootsrap values
        kon_boot_vec(n) = mean(ko_kon_vec(ko_ids_boot))/mean(wt_kon_vec(wt_ids_boot));
        dur_boot_vec(n) = mean(ko_dur_vec(ko_ids_boot))/mean(wt_dur_vec(wt_ids_boot));
        r_boot_vec(n) = mean(ko_r_vec(ko_ids_boot))/mean(wt_r_vec(wt_ids_boot));
    end
    hyp_vec_kon(i) = fv;
    if fv>1
        p_vec_kon(i) = mean(kon_boot_vec<=1);
    else
        p_vec_kon(i) = mean(kon_boot_vec>1);
    end
    hyp_vec_dur(i) = dv;
    if dv>1
        p_vec_dur(i) = mean(dur_boot_vec<=1);
    else
        p_vec_dur(i) = mean(dur_boot_vec>1);
    end
    hyp_vec_r(i) = rv;
    if rv>1
        p_vec_r(i) = mean(r_boot_vec<=1);
    else
        p_vec_r(i) = mean(r_boot_vec>1);
    end

end

%%  generate summary datasets
gene_cell = {'oct4','oct4','nanog','nanog'};

ko_flags = [zero_struct.additionalGroupVec];
summary_array = cat(2,ko_flags',...
                      [zero_struct.freq_vec_mean]', [zero_struct.freq_vec_ste]',....
                      [zero_struct.dur_vec_mean]', [zero_struct.dur_vec_ste]',...
                      [zero_struct.init_vec_mean]', [zero_struct.init_vec_ste]');
                    
summary_table = array2table(summary_array,'VariableNames',{'KO_flag','freq_mean' ,'freq_ste',...
                                                                         'dur_mean' ,'dur_ste',...
                                                                         'init_mean' ,'init_ste'});

summary_table.gene_name = gene_cell';
summary_table = [summary_table(:,end) summary_table(:,1:end-1)];
                                                                       
writetable(summary_table,[fig_path 'ko_wt_summary_table.csv']);

%% now add p value files
var_cell = {'frequency','duration','initiation'};
p_summary_array = vertcat(cat(2,hyp_vec_kon(1), p_vec_kon(1),hyp_vec_kon(2), p_vec_kon(2)),...
                               cat(2,hyp_vec_dur(1), p_vec_dur(1),hyp_vec_dur(2), p_vec_dur(2)),...
                               cat(2,hyp_vec_r(1), p_vec_r(1),hyp_vec_r(2), p_vec_r(2)));
                    
summary_table = array2table(p_summary_array,'VariableNames',{'oct4_val','oct4_p' ,'nanog_val','nanog_p'});

summary_table.var_name = var_cell';
summary_table = [summary_table(:,end) summary_table(:,1:end-1)];

writetable(summary_table,[fig_path 'pvalue_table.csv']);

