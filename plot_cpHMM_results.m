% script to plot results from cpHMM inference
clear
close all
addpath(genpath('utilities'))

projectNameCell = {'2xDl-Ven_snaBAC-mCh'};

% Set basic plotting parameters
xVar = 'fluo';

MarkerSize = 50;
blue = [115 143 193]/256;
purple = [171 133 172]/256;
red = [213 108 85]/256;

% set axes
dur_lims = [0 3.5];
freq_lims = [0 3];
init_lims = [0 18]*1e4;

for p = 1:length(projectNameCell)
    % set project to analyze 
    projectName = projectNameCell{p};

    % get path to results
    liveProject = LiveEnrichmentProject(projectName);
    resultsDir = [liveProject.dataPath filesep 'cpHMM_results' filesep];
    
    % make figure directory
    figureDir = [resultsDir 'figures' filesep];
    mkdir(figureDir);
    
    % get list of projects
    resultList = dir([resultsDir '*result*']);            
    
    for r = 1:length(resultList)
        % load data
        load([resultsDir filesep resultList(r).name]);
      
        % get index of x axis
        dashes = strfind(resultList(r).name,'_');
        additionalVar = resultList(r).name(dashes(end)+1:strfind(resultList(r).name,'.mat')-1);
        indexVarCols = {'AP','Time','fluo', additionalVar}; %NL: this will eventually be dynamic
        xIndex = find(strcmp(indexVarCols,xVar)); 
        
        % extract index variable array
        indexVarArray = compiledResults.inferenceOptions.indexInfo.indexVarArray;
        % figure out how many unique groups we need to plot
        [newVarArray, mapTo, indexList] = unique(indexVarArray(compiledResults.groupID_index,[1:xIndex-1 xIndex+1:end]),'rows');
        % see which variables actually change
        usedGroupers = [length(compiledResults.inferenceOptions.apBins)>2 length(compiledResults.inferenceOptions.timeBins)>2 ... 
          ~isempty(compiledResults.inferenceOptions.intensityBinVar) ~isempty(compiledResults.inferenceOptions.AdditionalGroupingVariable)];                                        
        
        lgd_flag = 0;
        if sum(usedGroupers)>2
          error('Code currently does not support plots with more than 1 additional grouping variable')
        elseif sum(usedGroupers)==2
          gpIndex = find(~ismember(1:length(indexVarCols),xIndex)&usedGroupers);
          gpVar = indexVarCols{gpIndex};
          lgd_flag = 1;
          lgd_prefix = [gpVar ' '];      
          if gpIndex > xIndex
            gpIndex = gpIndex-1;
          end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make figures
        close all
        
        r_trend = figure;
        hm_cm = flipud(brewermap(size(newVarArray,1)+2,'Set1'));
        colormap(hm_cm);
        lgd_str = {};
        hold on
        for i = 1:size(newVarArray,1)
            ind_list = find(indexList==i);
            e = errorbar(compiledResults.protein_intensity_vec(ind_list),compiledResults.init_vec_mean(ind_list),compiledResults.init_vec_ste(ind_list),'o','Color','black','LineWidth',1);          
            e.CapSize = 0;
        end
        s = [];
        for i = 1:size(newVarArray,1)
            ind_list = find(indexList==i);
            s(i) = scatter(compiledResults.protein_intensity_vec(ind_list),compiledResults.init_vec_mean(ind_list),MarkerSize,'o','MarkerFaceColor',hm_cm(i+1,:),'MarkerEdgeColor','black');
            if lgd_flag
              lgd_str(i) = {[lgd_prefix num2str(newVarArray(i,gpIndex))]};
            end
        end        
        grid on
%         xlim(x_lim)
        % ylim([50 95])
        xlabel('average Dorsal concentration (au) (au)') % NL: need to make this dynamic
        ylabel('burst amplitude (au/min)')
        if lgd_flag
            legend(s,lgd_str{:},'Location','southeast')
        end
        title(['Burst Amplitude (r): ' projectName])
        set(gca,'Fontsize',14)
%         ylim(init_lims);
%         StandardFigure([],gca)
        box on
        saveas(r_trend,[figureDir, 'burst_amp_' additionalVar '.tif'])
        saveas(r_trend,[figureDir, 'burst_amp_' additionalVar '.pdf'])
        %

        dur_trend = figure;        
                
        hold on
        for i = 1:size(newVarArray,1)
            ind_list = find(indexList==i);
            e = errorbar(compiledResults.protein_intensity_vec(ind_list),compiledResults.dur_vec_mean(ind_list),compiledResults.dur_vec_ste(ind_list),'o','Color','black','LineWidth',1);          
            e.CapSize = 0;
        end
        s = [];
        for i = 1:size(newVarArray,1)
            ind_list = find(indexList==i);
            s(i) = scatter(compiledResults.protein_intensity_vec(ind_list),compiledResults.dur_vec_mean(ind_list),MarkerSize,'o','MarkerFaceColor',hm_cm(i+1,:),'MarkerEdgeColor','black');   
        end        
        grid on
        title(['Burst Duration (1/k_{off}): ' projectName])
        xlabel('average Dorsal concentration (au) (au)') % NL: need to make this dynamic
        ylabel('burst duration (min)')
        if lgd_flag
            legend(s,lgd_str{:},'Location','southeast')
        end
        
%         ylim(dur_lims);
        set(gca,'Fontsize',14)
%         StandardFigure([],gca)
        box on
        saveas(dur_trend,[figureDir,'burst_dur_' additionalVar '.tif'])
        saveas(dur_trend,[figureDir,'burst_dur_' additionalVar '.pdf']);

        %   

        freq_trend = figure;
        
        hold on
        for i = 1:size(newVarArray,1)
            ind_list = find(indexList==i);
            e = errorbar(compiledResults.protein_intensity_vec(ind_list),compiledResults.freq_vec_mean(ind_list),compiledResults.freq_vec_ste(ind_list),'o','Color','black','LineWidth',1);          
            e.CapSize = 0;
        end
        s = [];
        for i = 1:size(newVarArray,1)
            ind_list = find(indexList==i);
            s(i) = scatter(compiledResults.protein_intensity_vec(ind_list),compiledResults.freq_vec_mean(ind_list),MarkerSize,'o','MarkerFaceColor',hm_cm(i+1,:),'MarkerEdgeColor','black');   
        end        
        grid on
        title(['Burst Frequency (k_{on}): ' projectName])
        xlabel('average Dorsal concentration (au) (au)') % NL: need to make this dynamic
        ylabel('burst frequency (1/min)')
        if lgd_flag
            legend(s,lgd_str{:},'Location','southeast')
        end        
%         StandardFigure([],gca)
        set(gca,'Fontsize',14)
%         ylim(freq_lims)
        box on
        
        saveas(freq_trend,[figureDir,'burst_freq_' additionalVar '.tif']);
        saveas(freq_trend,[figureDir,'burst_freq_' additionalVar '.pdf']);
    end
end