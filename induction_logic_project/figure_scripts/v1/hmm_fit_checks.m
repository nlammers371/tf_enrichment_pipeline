% Script to generate HMM fit check plots
clear 
close all
addpath('../utilities')
% set ID variables
DropboxFolder = 'S:\Nick\Dropbox\';
project_cell = {'Dl-Ven_snaBAC-mCh_v4','Dl-Ven_snaBAC-mCh_v4','Dl-Ven_snaBAC-mCh_F-F-F_v1','Dl-Ven_snaBAC-mCh_F-F-F_v1'};
title_cell = {'OG (old)','OG (new)'};
fluo_dim_vec = [2,3,2,3];
% Params
K = 3;
w = 7;
cmap = brewermap(9,'Set2');
rng(245)
n_traces = 100;
% load data for each project
master_struct = struct;
for p = 1:numel(project_cell)
    project = project_cell{p};
    fluo_dim = fluo_dim_vec(p);
    % set write paths
    [~, DataPath, FigureRoot] =   header_function(DropboxFolder, project); 
    FigPath = [FigureRoot project '\hmm_fit_checks\fluo_' num2str(fluo_dim_vec(p)) 'D\'];
    mkdir(FigPath)
    % intermediate HMM results
    load([DataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '_f' num2str(fluo_dim) 'D_dt.mat'])       
    trace_ids = randsample(1:numel(hmm_input_output),n_traces);
    
    for t = trace_ids
        time = hmm_input_output(t).time;
        fluo = hmm_input_output(t).fluo;
        fluo_hmm = hmm_input_output(t).fluo_hmm;
        
        chk_fig = figure('Visible','off');
        hold on
        plot(time,fluo,'LineWidth',1.5,'Color','black')
        plot(time,fluo_hmm,'LineWidth',1.5,'Color',cmap(2,:))
        
        xlabel('minutes into nc14')
        ylabel('MS2 spot intensity (au)')
        
        legend('raw data','hmm fit')
        set(gca,'FontSize',14)
        
        saveas(chk_fig,[FigPath 'trace_' num2str(t) '.png'])
        close all
    end
end
