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

%% Make figures
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

dose_fig = figure;

hold on
scatter(3,kon_vec_mean_z(2)/kon_vec_mean_z(1),100,'o','MarkerFaceColor',cmap(1,:),'MarkerEdgeColor','k')
scatter(2,dur_vec_mean_z(2)/dur_vec_mean_z(1),100,'s','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k')
scatter(1,r_vec_mean_z(2)/r_vec_mean_z(1),100,'^','MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k')

set(gca,'yscale','log')
grid on



kon_vec_mean_nz = [nz_struct.freq_vec_mean];
kon_vec_ste_nz = [nz_struct.freq_vec_ste];

dur_vec_mean_nz = [nz_struct.dur_vec_mean];
dur_vec_ste_nz = [nz_struct.dur_vec_ste];

r_vec_mean_nz = [nz_struct.init_vec_mean];
r_vec_ste_nz = [nz_struct.init_vec_ste];

dose_fig = figure;

hold on
scatter(3,kon_vec_mean_nz(2)/kon_vec_mean_nz(1),100,'o','MarkerFaceColor',cmap(1,:),'MarkerEdgeColor','k')
scatter(2,dur_vec_mean_nz(2)/dur_vec_mean_nz(1),100,'s','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k')
scatter(1,r_vec_mean_nz(2)/r_vec_mean_nz(1),100,'^','MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k')

set(gca,'yscale','log')
grid on
% errorbar(1:3,kon_vec_mean,kon_vec_ste,'.','Color','k','LineWidth',1.25)
% ylim([0 1.5*max(kon_vec_mean)])
%%

saveas(dose_fig,[dose_fig_path 'burst_bar_plot.pdf'])
saveas(dose_fig,[dose_fig_path 'burst_bar_plot.png'])

% opto
cmap2 = brewermap(12,'Greens');
cmap3 = brewermap(12,'Greys');

opto_fig_path = [resultsRoot filesep 'opto_figs_v2' filesep];
mkdir(opto_fig_path);

cmap = brewermap(10,'Blues');
% close all


opto_fig = figure('Position',[100 100 768 256]);
tiledlayout(1,3)
nexttile
hold on
bar(1,nz_struct(1).freq_vec_mean(1),'FaceColor',cmap2(2,:))
bar(2,nz_struct(2).freq_vec_mean(1),'FaceColor',cmap3(3,:))

bar(4,nz_struct(1).freq_vec_mean(2),'FaceColor',cmap2(5,:))
bar(5,nz_struct(2).freq_vec_mean(2),'FaceColor',cmap3(5,:))

bar(7,nz_struct(1).freq_vec_mean(3),'FaceColor',cmap2(8,:))
bar(8,nz_struct(2).freq_vec_mean(3),'FaceColor',cmap3(8,:))

errorbar([1 4 7],nz_struct(1).freq_vec_mean,nz_struct(1).freq_vec_ste,'.','Color','k','LineWidth',1.25)
errorbar([2 5 8],nz_struct(2).freq_vec_mean,nz_struct(2).freq_vec_ste,'.','Color','k','LineWidth',1.25)
% ylim([0 1.5*max(kon_vec_mean)])
ylabel('burst frequency (min^{-1})')
set(gca,'xtick',[])
set(gca,'Fontsize',12)
box on

nexttile
hold on
bar(1,nz_struct(1).dur_vec_mean(1),'FaceColor',cmap2(2,:))
bar(2,nz_struct(2).dur_vec_mean(1),'FaceColor',cmap3(3,:))

bar(4,nz_struct(1).dur_vec_mean(2),'FaceColor',cmap2(5,:))
bar(5,nz_struct(2).dur_vec_mean(2),'FaceColor',cmap3(5,:))

bar(7,nz_struct(1).dur_vec_mean(3),'FaceColor',cmap2(8,:))
bar(8,nz_struct(2).dur_vec_mean(3),'FaceColor',cmap3(8,:))

errorbar([1 4 7],nz_struct(1).dur_vec_mean,nz_struct(1).dur_vec_ste,'.','Color','k','LineWidth',1.25)
errorbar([2 5 8],nz_struct(2).dur_vec_mean,nz_struct(2).dur_vec_ste,'.','Color','k','LineWidth',1.25)
%ylim([0 2.6*dur_avg_low])
% ylim([0 1.5*max(dur_vec_mean)])
ylabel('burst duration (min)')
set(gca,'xtick',[])
set(gca,'Fontsize',12)
box on

nexttile
hold on
bar(1,nz_struct(1).init_vec_mean(1),'FaceColor',cmap2(2,:))
bar(2,nz_struct(2).init_vec_mean(1),'FaceColor',cmap3(3,:))

bar(4,nz_struct(1).init_vec_mean(2),'FaceColor',cmap2(5,:))
bar(5,nz_struct(2).init_vec_mean(2),'FaceColor',cmap3(5,:))

bar(7,nz_struct(1).init_vec_mean(3),'FaceColor',cmap2(8,:))
bar(8,nz_struct(2).init_vec_mean(3),'FaceColor',cmap3(8,:))

errorbar([1 4 7],nz_struct(1).init_vec_mean,nz_struct(1).init_vec_ste,'.','Color','k','LineWidth',1.25)
errorbar([2 5 8],nz_struct(2).init_vec_mean,nz_struct(2).init_vec_ste,'.','Color','k','LineWidth',1.25)
ylabel('initiation rate (au/min)')
% errorbar(1:3,r_vec_mean,r_vec_ste,'.','Color','k','LineWidth',1.25)
ylim([0 2e4])
set(gca,'xtick',[])
set(gca,'Fontsize',12)
box on

saveas(opto_fig,[opto_fig_path 'opto_burst_bar_plot.pdf'])
saveas(opto_fig,[opto_fig_path 'opto_burst_bar_plot.png'])
%%  calculate p-values

rng(13); % seed random number generator for consistency
n_boots = 1e4;

%%%%%%%%%%%%%%%%%%%%%%%%%
% dose response
%%%%%%%%%%%%%%%%%%%%%%%%%
p_vec_kon_dose = NaN(1,3);
p_vec_r_dose = NaN(1,3);
p_vec_dur_dose = NaN(1,3);

hyp_vec_kon_dose = NaN(1,3);
hyp_vec_r_dose = NaN(1,3);
hyp_vec_dur_dose = NaN(1,3);

% initialize reference vectors
high_outlier_flags = zero_struct(3).outlier_flags{1};
high_kon_vec = zero_struct(3).freq_results{1};
high_dur_vec = zero_struct(3).dur_results{1};
high_r_vec = zero_struct(3).r_results{1}(2,:);
high_options = find(~high_outlier_flags);

for i = 1:2
    % extract key vectors
    outlier_flags_opto = zero_struct(i).outlier_flags{1};
    kon_vec = zero_struct(i).freq_results{1};
    dur_vec = zero_struct(i).dur_results{1};
    r_vec = zero_struct(i).r_results{1}(2,:);
    temp_options = find(~outlier_flags_opto);
    % calculate value of mean fold qauntities
    fv = zero_struct(i).freq_vec_mean/zero_struct(3).freq_vec_mean;
    dv = zero_struct(i).dur_vec_mean/zero_struct(3).dur_vec_mean;
    rv = zero_struct(i).init_vec_mean/zero_struct(3).init_vec_mean;
    % initialize arrays
    kon_boot_vec = NaN(1,n_boots);
    dur_boot_vec = NaN(1,n_boots);
    r_boot_vec = NaN(1,n_boots);
    % perform bootsrapping
    n_samp = min([length(temp_options) length(high_options)]);
    for n = 1:n_boots
        % draw samples
        high_ids_boot = randsample(high_options,n_samp,true);
        temp_ids_boot_opto = randsample(temp_options,n_samp,true);
        % calculate bootsrap values
        kon_boot_vec(n) = mean(kon_vec(temp_ids_boot_opto))/mean(high_kon_vec(high_ids_boot));
        dur_boot_vec(n) = mean(dur_vec(temp_ids_boot_opto))/mean(high_dur_vec(high_ids_boot));
        r_boot_vec(n) = mean(r_vec(temp_ids_boot_opto))/mean(high_r_vec(high_ids_boot));
    end
    hyp_vec_kon_dose(i) = fv;
    if fv>1
        p_vec_kon_dose(i) = mean(kon_boot_vec<=1);
    else
        p_vec_kon_dose(i) = mean(kon_boot_vec>1);
    end
    hyp_vec_dur_dose(i) = dv;
    if dv>1
        p_vec_dur_dose(i) = mean(dur_boot_vec<=1);
    else
        p_vec_dur_dose(i) = mean(dur_boot_vec>1);
    end
    hyp_vec_r_dose(i) = rv;
    if rv>1
        p_vec_r_dose(i) = mean(r_boot_vec<=1);
    else
        p_vec_r_dose(i) = mean(r_boot_vec>1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% opto response
%%%%%%%%%%%%%%%%%%%%%%%%%
p_vec_kon_opto = NaN(1,3);
p_vec_r_opto = NaN(1,3);
p_vec_dur_opto = NaN(1,3);

hyp_vec_kon_opto = NaN(1,3);
hyp_vec_r_opto = NaN(1,3);
hyp_vec_dur_opto = NaN(1,3);

% extract reference vectors
dark_outlier_flags_opto = nz_struct(1).outlier_flags{1};
dark_kon_vec_opto = nz_struct(1).freq_results{1};
dark_dur_vec_opto = nz_struct(1).dur_results{1};
dark_r_vec_opto = nz_struct(1).r_results{1}(2,:);
dark_options_opto = find(~dark_outlier_flags_opto);

dark_outlier_flags_ctrl = nz_struct(2).outlier_flags{1};
dark_kon_vec_ctrl = nz_struct(2).freq_results{1};
dark_dur_vec_ctrl = nz_struct(2).dur_results{1};
dark_r_vec_ctrl = nz_struct(2).r_results{1}(2,:);
dark_options_ctrl = find(~dark_outlier_flags_ctrl);

for i = 2:3
    % extract key vectors (opto)
    outlier_flags_opto = nz_struct(1).outlier_flags{i};
    kon_vec_opto = nz_struct(1).freq_results{i};
    dur_vec_opto = nz_struct(1).dur_results{i};
    r_vec_opto = nz_struct(1).r_results{i}(2,:);
    temp_options_opto = find(~outlier_flags_opto);

    % extract key vectors (ctrl)
    outlier_flags_ctrl = nz_struct(2).outlier_flags{i};
    kon_vec_ctrl = nz_struct(2).freq_results{i};
    dur_vec_ctrl = nz_struct(2).dur_results{i};
    r_vec_ctrl = nz_struct(2).r_results{i}(2,:);
    temp_options_ctrl = find(~outlier_flags_ctrl);
    % calculate value of mean fold qauntities
    fv = (nz_struct(1).freq_vec_mean(i)/nz_struct(2).freq_vec_mean(i))/(nz_struct(1).freq_vec_mean(1)/nz_struct(2).freq_vec_mean(1));
    dv = (nz_struct(1).dur_vec_mean(i)/nz_struct(2).dur_vec_mean(i))/(nz_struct(1).dur_vec_mean(1)/nz_struct(2).dur_vec_mean(1));
    rv = (nz_struct(1).init_vec_mean(i)/nz_struct(2).init_vec_mean(i))/(nz_struct(1).init_vec_mean(1)/nz_struct(2).init_vec_mean(1));
    % initialize arrays
    kon_boot_vec = NaN(1,n_boots);
    dur_boot_vec = NaN(1,n_boots);
    r_boot_vec = NaN(1,n_boots);
    % perform bootsrapping
    for n = 1:n_boots
        % draw samples
        dark_ids_boot_opto = randsample(dark_options_opto,length(dark_options_opto),true);
        dark_ids_boot_ctrl = randsample(dark_options_ctrl,length(dark_options_ctrl),true);
        temp_ids_boot_opto = randsample(temp_options_opto,length(temp_options_opto),true);
        temp_ids_boot_ctrl = randsample(temp_options_ctrl,length(temp_options_ctrl),true);
        % calculate bootsrap values
        kon_boot_vec(n) = (mean(kon_vec_opto(temp_ids_boot_opto))/mean(kon_vec_ctrl(temp_ids_boot_ctrl)))/...
                          (mean(dark_kon_vec_opto(dark_ids_boot_opto))/mean(dark_kon_vec_ctrl(dark_ids_boot_ctrl)));

        dur_boot_vec(n) = (mean(dur_vec_opto(temp_ids_boot_opto))/mean(dur_vec_ctrl(temp_ids_boot_ctrl)))/...
                          (mean(dark_dur_vec_opto(dark_ids_boot_opto))/mean(dark_dur_vec_ctrl(dark_ids_boot_ctrl)));

        r_boot_vec(n) = (mean(r_vec_opto(temp_ids_boot_opto))/mean(r_vec_ctrl(temp_ids_boot_ctrl)))/...
                          (mean(dark_r_vec_opto(dark_ids_boot_opto))/mean(dark_r_vec_ctrl(dark_ids_boot_ctrl)));
    end
    hyp_vec_kon_opto(i) = fv;
    if fv>1
        p_vec_kon_opto(i) = mean(kon_boot_vec<=1);
    else
        p_vec_kon_opto(i) = mean(kon_boot_vec>1);
    end
    hyp_vec_dur_opto(i) = dv;
    if dv>1
        p_vec_dur_opto(i) = mean(dur_boot_vec<=1);
    else
        p_vec_dur_opto(i) = mean(dur_boot_vec>1);
    end
    hyp_vec_r_opto(i) = rv;
    if rv>1
        p_vec_r_opto(i) = mean(r_boot_vec<=1);
    else
        p_vec_r_opto(i) = mean(r_boot_vec>1);
    end
end

%%  generate summary datasets

opto_summary_array = cat(2,double(nz_struct(1).timeGroupVec'),...
                      nz_struct(1).freq_vec_mean', nz_struct(1).freq_vec_ste', nz_struct(2).freq_vec_mean',nz_struct(2).freq_vec_ste',....
                      nz_struct(1).dur_vec_mean', nz_struct(1).dur_vec_ste', nz_struct(2).dur_vec_mean', nz_struct(2).dur_vec_ste',...
                      nz_struct(1).init_vec_mean', nz_struct(1).init_vec_ste', nz_struct(2).init_vec_mean', nz_struct(2).init_vec_ste');
                    
opto_table = array2table(opto_summary_array,'VariableNames',{'time_group','freq_mean_opto' ,'freq_ste_opto','freq_mean_ctrl' ,'freq_ste_ctrl',...
                                                                         'dur_mean_opto' ,'dur_ste_opto','dur_mean_ctrl' ,'dur_ste_ctrl',...
                                                                         'init_mean_opto' ,'init_ste_opto','init_mean_ctrl' ,'init_ste_ctrl'});
                                                                       
writetable(opto_table,[opto_fig_path 'opto_summary_table.csv']);


dose_summary_array = cat(2,[1:3]',...
                      [zero_struct.freq_vec_mean]', [zero_struct.freq_vec_ste]',....
                      [zero_struct.dur_vec_mean]', [zero_struct.dur_vec_ste]',...
                      [zero_struct.init_vec_mean]', [zero_struct.init_vec_ste]');
                    
dose_table = array2table(dose_summary_array,'VariableNames',{'dose_group','freq_mean' ,'freq_ste',...
                                                                         'dur_mean' ,'dur_ste',...
                                                                         'init_mean' ,'init_ste'});
                                                                       
writetable(dose_table,[dose_fig_path 'dose_summary_table.csv']);

% now add p value files
var_cell = {'frequency','duration','initiation'};
opto_p_summary_array = vertcat(cat(2,hyp_vec_kon_opto(2), p_vec_kon_opto(2),hyp_vec_kon_opto(3), p_vec_kon_opto(3)),...
                               cat(2,hyp_vec_dur_opto(2), p_vec_dur_opto(2),hyp_vec_dur_opto(3), p_vec_dur_opto(3)),...
                               cat(2,hyp_vec_r_opto(2), p_vec_r_opto(2),hyp_vec_r_opto(3), p_vec_r_opto(3)));
                    
opto_p_table = array2table(opto_p_summary_array,'VariableNames',{'opto_mid_val','opto_mid_p' ,'opto_late_val','opto_late_p'});

opto_p_table.var_name = var_cell';
opto_p_table = [opto_p_table(:,end) opto_p_table(:,1:end-1)];

writetable(opto_p_table,[opto_fig_path 'opto_pvalue_table.csv']);

dose_p_summary_array = vertcat(cat(2,hyp_vec_kon_dose(1), p_vec_kon_dose(1),hyp_vec_kon_dose(2), p_vec_kon_dose(2)),...
                               cat(2,hyp_vec_dur_dose(1), p_vec_dur_dose(1),hyp_vec_dur_dose(2), p_vec_dur_dose(2)),...
                               cat(2,hyp_vec_r_dose(1), p_vec_r_dose(1),hyp_vec_r_dose(2), p_vec_r_dose(2)));

dose_p_table = array2table(dose_p_summary_array,'VariableNames',{'dose_low_val','dose_low_p' ,'dose_mid_val','dose_mid_p'});
 
dose_p_table.var_name = var_cell';
dose_p_table = [dose_p_table(:,end) dose_p_table(:,1:end-1)];

writetable(dose_p_table,[dose_fig_path 'dose_pvalue_table.csv']);
