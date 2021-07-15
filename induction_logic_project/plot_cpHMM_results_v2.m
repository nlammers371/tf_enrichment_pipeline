% script to plot results from cpHMM inference
clear
close all
addpath(genpath('../utilities'))


% script to compile results from cpHMM inference
clear
close all
addpath(genpath('utilities'))

% projectNameCell = {'EveGtSL','EveGtSL-S1Null','EveWt','EveS1Null'};%};
projectNameCell = {'20210430_Nanog','20210430_Oct4','20210430_Sox2'};%};
% resultsRoot = 'S:\Nick\Dropbox\InductionLogic\';
condition_key = {'SFES (WT)','SFES (KO)','diff (WT)','diff (KO)'};
% useful colors
MarkerSize = 100;
blue = [115 143 193]/256;
purple = [171 133 172]/256;
red = [213 108 85]/256;

% set axes
dur_lims = [0 3.5];
freq_lims = [0 3];
init_lims = [0 18]*1e4;
    
for p = 3%:length(projectNameCell)
    
    % set project to analyze 
    projectName = projectNameCell{p};

    % get path to results
    try
        liveProject = LiveEnrichmentProject(projectName);
        resultsDir = [liveProject.dataPath 'cpHMM_results' filesep];
    catch
        resultsRoot = 'S:/Nick/Dropbox/ProcessedEnrichmentData/';
        resultsDir = [resultsRoot filesep projectNameCell{p} filesep 'cpHMM_results' filesep];
    end
    

    % make figure path
    FigPath = [resultsDir 'inference_figs' filesep];
    mkdir(FigPath);   
        
    
    % get list of projects
    resultList = dir([resultsDir '*result*']);            
    
    for r = 1%:length(resultList) % assume just 1 for now
      
        % load data
        load([resultsDir filesep resultList(r).name]);
                  
        % transfer results
        init_vec_mean = compiledResults.init_vec_mean;
        init_vec_ste = compiledResults.init_vec_ste;
        
        freq_vec_mean = compiledResults.freq_vec_mean;
        freq_vec_ste = compiledResults.freq_vec_ste;
        
        dur_vec_mean = compiledResults.dur_vec_mean;
        dur_vec_ste = compiledResults.dur_vec_ste;
        
        time_vec_mean = compiledResults.time_vec_mean;
        time_vec_ste = compiledResults.time_vec_ste;
        
        fluo_vec_mean = compiledResults.fluo_mean;
        fluo_vec_ste = compiledResults.fluo_ste;
        
    end
    
    x_vec = 1:length(fluo_vec_mean);
    % make figures
    close all
    fluo_trend = figure;
    hm_cm = brewermap(9,'Set2');
    colormap(hm_cm);

    hold on
    
    e = errorbar(x_vec,fluo_vec_mean,fluo_vec_ste,'o','Color','black','LineWidth',1);  
    e.CapSize = 0;
    s = scatter(x_vec,fluo_vec_mean,MarkerSize,'o','MarkerFaceColor',hm_cm(5,:),'MarkerEdgeColor','black');
    
    grid on
    %         xlim(x_lim)
    ylim([0.66*min(fluo_vec_mean) 1.5*max(fluo_vec_mean)])
    % xlabel('time') % NL: need to make this dynamic
    set(gca,'xticklabel',condition_key)
    ylabel('average transcription rate (au/min)')
%     legend(s,'Sox2 (wildtype control)','Sox2 (opto-chronic)','Location','southwest')
    title(['Mean Transcription Rate (r) Across Conditions'])
    set(gca,'Fontsize',14)
    % ylim([2000 3000]);
    xlim([x_vec(1)-0.5 x_vec(end)+0.5])
    set(gca,'yscale','log');
%     ylim([ 3e3 8e3])
    box on
    xtickangle(45)
    saveas(fluo_trend,[FigPath, 'spot_fluo.tif'])
    saveas(fluo_trend,[FigPath, 'spot_fluo.pdf'])
           
    
    r_trend = figure;
    hm_cm = brewermap(9,'Set2');
    colormap(hm_cm);

    hold on
    
    e = errorbar(x_vec,init_vec_mean,init_vec_ste,'o','Color','black','LineWidth',1);  
    e.CapSize = 0;
    s = scatter(x_vec,init_vec_mean,MarkerSize,'o','MarkerFaceColor',hm_cm(1,:),'MarkerEdgeColor','black');
    
    grid on
    %         xlim(x_lim)
    ylim([0.66*min(init_vec_mean) 1.5*max(init_vec_mean)])
    % xlabel('time') % NL: need to make this dynamic
    set(gca,'xticklabel',condition_key)
    ylabel('burst amplitude (au/min)')
%     legend(s,'Sox2 (wildtype control)','Sox2 (opto-chronic)','Location','southwest')
    title(['Burst Amplitude (r) Across Conditions'])
    set(gca,'Fontsize',14)
    % ylim([2000 3000]);
    xlim([x_vec(1)-0.5 x_vec(end)+0.5])
    set(gca,'yscale','log');
%     ylim([ 3e3 8e3])
    box on
    xtickangle(45)
    saveas(r_trend,[FigPath, 'burst_amp.tif'])
    saveas(r_trend,[FigPath, 'burst_amp.pdf'])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% duration
    dur_trend = figure;        

    hold on
    
    e = errorbar(x_vec,dur_vec_mean,dur_vec_ste,'o','Color','black','LineWidth',1);  
    e.CapSize = 0;
    s = scatter(x_vec,dur_vec_mean,MarkerSize,'o','MarkerFaceColor',hm_cm(2,:),'MarkerEdgeColor','black');
    
    grid on
    %         xlim(x_lim)
    ylim([0.66*min(dur_vec_mean) 1.5*max(dur_vec_mean)])
    % xlabel('time') % NL: need to make this dynamic
    set(gca,'xticklabel',condition_key)
    ylabel('burst duration (min)')
%     legend(s,'Sox2 (wildtype control)','Sox2 (opto-chronic)','Location','southwest')
    title(['Burst Duration (1/k_{off}) Across Conditions'])
    set(gca,'Fontsize',14)
    % ylim([2000 3000]);
    xlim([x_vec(1)-0.5 x_vec(end)+0.5])
    set(gca,'yscale','log');
%     ylim([ 3e3 8e3])
    box on
    xtickangle(45)
    saveas(dur_trend,[FigPath, 'burst_dur.tif'])
    saveas(dur_trend,[FigPath, 'burst_dur.pdf'])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    %%%%%% Frequency
    freq_trend = figure;        

    hold on
    
    e = errorbar(x_vec,freq_vec_mean,freq_vec_ste,'o','Color','black','LineWidth',1);  
    e.CapSize = 0;
    s = scatter(x_vec,freq_vec_mean,MarkerSize,'o','MarkerFaceColor',hm_cm(3,:),'MarkerEdgeColor','black');
    
    grid on
    %         xlim(x_lim)
    ylim([0.66*min(freq_vec_mean) 1.5*max(freq_vec_mean)])
    % xlabel('time') % NL: need to make this dynamic
    set(gca,'xticklabel',condition_key)
    ylabel('burst frequency (1/min)')
%     legend(s,'Sox2 (wildtype control)','Sox2 (opto-chronic)','Location','southwest')
    title(['Burst Frequency (k_{on}) Across Conditions'])
    set(gca,'Fontsize',14)
    % ylim([2000 3000]);
    xlim([x_vec(1)-0.5 x_vec(end)+0.5])
    set(gca,'yscale','log');
%     ylim([ 3e3 8e3])
    box on
    xtickangle(45)
    saveas(freq_trend,[FigPath, 'burst_freq.tif'])
    saveas(freq_trend,[FigPath, 'burst_freq.pdf'])
    
    % make fold-change figure
    % duration
    init_fold = init_vec_mean([2 4])./init_vec_mean([1 3]);
    err_1 = (1./init_vec_mean([1 3])).*init_vec_ste([2 4]);
    err_2 = (init_vec_ste([2 4])./ init_vec_mean([1 3]).^2).*init_vec_ste([1 3]);
    init_fold_ste = sqrt(err_1.^2 + err_2.^2);
    
    % duration
    dur_fold = dur_vec_mean([2 4])./dur_vec_mean([1 3]);
    err_1 = (1./dur_vec_mean([1 3])).*dur_vec_ste([2 4]);
    err_2 = (dur_vec_ste([2 4])./ dur_vec_mean([1 3]).^2).*dur_vec_ste([1 3]);
    dur_fold_ste = sqrt(err_1.^2 + err_2.^2);
    
    % frequency
    freq_fold = freq_vec_mean([2 4])./freq_vec_mean([1 3]);
    err_1 = (1./freq_vec_mean([1 3])).*freq_vec_ste([2 4]);
    err_2 = (freq_vec_ste([2 4])./freq_vec_mean([1 3]).^2).*freq_vec_ste([1 3]);
    freq_fold_ste = sqrt(err_1.^2 + err_2.^2);
    
    % fluorescence
    fluo_fold = fluo_vec_mean([2 4])./fluo_vec_mean([1 3]);
    err_1 = (1./fluo_vec_mean([1 3])).*fluo_vec_ste([2 4]);
    err_2 = (fluo_vec_ste([2 4])./fluo_vec_mean([1 3]).^2).*fluo_vec_ste([1 3]);
    fluo_fold_ste = sqrt(err_1.^2 + err_2.^2);
    
    fold_fig = figure;
    hold on
    e = errorbar((1:2)-0.1,init_fold,init_fold_ste,'.','Color','black','LineWidth',1);  
    e.CapSize = 0;
    s1 = scatter((1:2)-0.1,init_fold,MarkerSize,'o','MarkerFaceColor',hm_cm(1,:),'MarkerEdgeColor','black');
    
    e = errorbar((1:2)-0.05,dur_fold,dur_fold_ste,'.','Color','black','LineWidth',1);  
    e.CapSize = 0;
    s2 = scatter((1:2)-0.05,dur_fold,MarkerSize,'s','MarkerFaceColor',hm_cm(2,:),'MarkerEdgeColor','black');
    
    e = errorbar((1:2)+0.05,freq_fold,freq_fold_ste,'.','Color','black','LineWidth',1);  
    e.CapSize = 0;
    s3 = scatter((1:2)+0.05,freq_fold,MarkerSize,'d','MarkerFaceColor',hm_cm(3,:),'MarkerEdgeColor','black');
    
    e = errorbar((1:2)+0.1,fluo_fold,fluo_fold_ste,'.','Color','black','LineWidth',1);  
    e.CapSize = 0;
    s4 = scatter((1:2)+0.1,fluo_fold,MarkerSize,'^','MarkerFaceColor',hm_cm(5,:),'MarkerEdgeColor','black');
    
    set(gca,'Fontsize',14)
    % ylim([2000 3000]);
    xlim([0.5 2.5])
    set(gca,'xtick',1:2,'xticklabels',{'SFES','diff'})
    set(gca,'yscale','log');
    legend([s1 s2 s3 s4],'burst amplitude','burst duration','burst frequency','mean rate','Location','best')
    ylabel('fold change')
    grid on
    
    saveas(fold_fig,[FigPath, 'fold_change.tif'])
    saveas(fold_fig,[FigPath, 'fold_change.pdf'])
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make figures        
close all

