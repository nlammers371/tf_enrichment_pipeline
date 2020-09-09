clear
close all

typeName = 'eNos-MCP-mCh';

quantiles_in = [15:15:90 98]/100; % specify quantiles to use for analysis
tlims_in = [5 25]; % specify time limits (roughly when the spots should be ON)
nc14_flag_in = true;

% Define basic ID variables
DropboxFolder =  'S:\Meghan\Dropbox\';
figPath = [DropboxFolder 'ModularEnhancersResults' filesep 'eNos-MCP-mCh_saturation' filesep];
mkdir(figPath);

% prefixes = {...
%     '2017-09-14-P2P-MCP-NoNLS-mCherry-doubledosage2',...
%     '2017-09-14-P2P-MCP-NoNLS-mCherry-doubledosage3'};
prefixes = {};

if isempty(prefixes) 
    projects = {'eNosx1-MCP-mCh_VK33_het','eNosx1-MCP-mCh_VK33_homo','eNosx2-MCP-mCh_VK33_homo','eNosx2-MCP-mCh_VK22'};
    projectData(numel(projects)).constructName = '';  %storing in a struct so it's organized by eNos construct
    
    for p = 1:numel(projects)
        tempPrefixes = getProjectPrefixes(projects{p},'onlyApproved');
        projectData(p).constructName = projects{p};
        projectData(p).prefixes = tempPrefixes;
    end
end

% call function to calculate offset and fluo quantiles for each prefix
disp('calculating quantiles...')
for i = 1:numel(projects)
    currNumPrefixes = numel(projectData(i).prefixes);
    
    % intialize data arrays
    spot_mean_array = NaN(numel(quantiles_in),currNumPrefixes);
    spot_se_array = NaN(numel(quantiles_in),currNumPrefixes);

    offset_mean_array = NaN(numel(quantiles_in),currNumPrefixes);
    offset_se_array = NaN(numel(quantiles_in),currNumPrefixes);
    
    for j = 1:currNumPrefixes
        currPrefix = projectData(i).prefixes{j};    %needs to be a string
        disp(['analyzing ' currPrefix '...'])
        [spot_mean_array(:,j), spot_se_array(:,j), offset_mean_array(:,j),...
                offset_se_array(:,j), nc14_flag, quantiles, tlims] = ...
            mcp_saturation_analysis(currPrefix,'quantiles',quantiles_in, ...
                'nc14_flag',nc14_flag_in,'tlims',tlims_in);
                
    end
    
    projectData(i).spotMean = spot_mean_array;
    projectData(i).spotSE = spot_se_array;
    projectData(i).offsetMean = offset_mean_array;
    projectData(i).offsetSE = offset_se_array;
    projectData(i).quantiles = quantiles_in;
    projectData(i).nc14Flag = nc14_flag_in;
    projectData(i).tlims = tlims_in;
end

%% Make spot fluorescence vs MCP offset figure

% PBoC colors
pbocRed = [213,108,85]/255;
pbocBlu = [115,142,193]/255;
pbocYlw = [234, 194,100]/255;
pbocGrn = [122,169,116]/255;
pbocPur = [171,133,172]/255;

pbocColors = [pbocRed ; pbocBlu ; pbocYlw ; pbocGrn ; pbocPur];
ind98pct = 7;

eNosLabels = {'eNosx1 @VK33, het','eNosx1 @VK33, homo','eNosx2 @VK33','eNosx2 @VK22'};

spotVsOffset98pctFig = figure;
hold on
for p = 1:numel(projectData)
    
    nPrefixes = size(projectData(p).offsetMean,2);
    
    offsetMean98pct = projectData(p).offsetMean(ind98pct,:);
    spotMean98pct = projectData(p).spotMean(ind98pct,:);
    
    s = plot(offsetMean98pct, spotMean98pct,'o','MarkerSize',10,'MarkerFaceColor',pbocColors(p,:),'MarkerEdgeColor','black');
end

% ax = gca;
% ax.XLim = [0 5];
% ax.YLim = [0.5 2];
% ax.XTick = [1, 2, 3, 4];
% ax.XTickLabels = {'eNosx1-MCP-mCh @VK33, het','eNosx1-MCP-mCh @VK33 homo','eNosx2-MCP-mCh @VK33','eNosx2-MCP-mCh @VK22'};
legend(eNosLabels(2:end),'Location','northwest');
xlabel('MCP offset (au)')
ylabel('MS2 spot fluo (au)')
box on
hold off

saveas(spotVsOffset98pctFig,[figPath typeName '_spotFluoVsOffset.png'])
saveas(spotVsOffset98pctFig,[figPath typeName '_spotFluoVsOffset.pdf'])

%% Make 98th percentile MCP offset (~nuclear MCP concentration) vs construct figure

mcpOffset98pctFig = figure;
hold on
for p = 1:numel(projectData)
    
    nPrefixes = size(projectData(p).offsetMean,2);
    xData = p * ones(1, nPrefixes);
    
    offsetMean98pct = projectData(p).offsetMean(ind98pct,:);
    offsetSE98pct = projectData(p).offsetSE(ind98pct,:);
    
    s = plot(xData, offsetMean98pct,'o','Color',pbocColors(p,:),'MarkerSize',10,'MarkerFaceColor',pbocColors(p,:),'MarkerEdgeColor','black');
end

ax = gca;
ax.XLim = [0 5];
ax.YLim = [0.5 2];
ax.XTick = [1, 2, 3, 4];
ax.XTickLabels = eNosLabels;
ax.XTickLabelRotation = 45;
ax.FontSize = 12;
ylabel('98th pct MCP offset (au)')
box on
hold off

saveas(mcpOffset98pctFig,[figPath typeName '_mcpOffset98pct.png'])
saveas(mcpOffset98pctFig,[figPath typeName '_mcpOffset98pct.pdf'])


%% Make 98th percentile MS2 spot fluorescence vs construct figure

mcpSpotFluo98pctFig = figure;
hold on
for p = 1:numel(projectData)
    
    nPrefixes = size(projectData(p).spotMean,2);
    xData = p * ones(1, nPrefixes);
    
    spotMean98pct = projectData(p).spotMean(ind98pct,:);
    spotSE98pct = projectData(p).spotSE(ind98pct,:);
    
    e = plot(xData, spotMean98pct,'o','Color',pbocColors(p,:),'MarkerSize',10,'MarkerFaceColor',pbocColors(p,:),'MarkerEdgeColor','black');
end

ax = gca;
ax.XLim = [0 5];
ax.YLim = [300 600];
ax.XTick = [1, 2, 3, 4];
ax.XTickLabels = eNosLabels;
ax.XTickLabelRotation = 45;
ax.FontSize = 12;
ylabel('98th pct MS2 spot fluo. (au)')
box on
hold off

saveas(mcpSpotFluo98pctFig,[figPath typeName '_spotFluo98pct.png'])
saveas(mcpSpotFluo98pctFig,[figPath typeName '_spotFluo98pct.pdf'])


%% save analysis details
save([figPath typeName '_saturationData.mat'],'projectData')
