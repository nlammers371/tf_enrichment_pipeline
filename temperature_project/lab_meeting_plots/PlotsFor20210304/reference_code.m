function PlotLTMSubplotsTrapParamsVsTemp(this, parameter, outdir, varargin)
%%

% PlotTitle, PlottingColors, UseDifferentColors,
% UseDiffProfiles, UsePhysicalAPLength
IncludeFits = true;
x = 1;
while x <= length(varargin)
    if strcmpi(varargin{x}, 'plottitle')
        PlotTitle = varargin{x+1};
        x = x+1;
    elseif strcmpi(varargin{x}, 'plottingcolors')
        PlottingColors = varargin{x+1};
        x = x+1;
    elseif strcmpi(varargin{x}, 'tracetype')
        TraceType = lower(varargin{x+1});
        x = x+1;
    elseif strcmpi(varargin{x}, 'ExcludeFits')
        IncludeFits = false;
    end
    x = x+1;
end


if ~exist('PlottingColors', 'var')
    PlottingColors = 'default';
elseif ~strcmp(lower(PlottingColors), 'default')  & ~strcmp(lower(PlottingColors), 'pboc')
    error('Invalid choice of plotting colors. Can use either "default" or "pboc".') % change to error
end
if ~exist('TraceType', 'var')
    TraceType = 'AnaphaseAligned';
elseif strcmpi(TraceType, 'Fluo3D') | strcmpi(TraceType, 'Unaligned3D')
    TraceType = 'Unaligned3D';
elseif strcmpi(TraceType, 'Fluo')| strcmpi(TraceType, 'Unaligned')
    TraceType = 'Unaligned';
elseif strcmpi(TraceType, 'AnaphaseAligned')
    TraceType = 'AnaphaseAligned';
elseif strcmpi(TraceType, 'AnaphaseAligned3D')
    TraceType = 'AnaphaseAligned3D';
elseif strcmpi(TraceType, 'Tbinned')
    TraceType = 'Tbinned';
elseif strcmpi(TraceType, 'Tbinned3D')
    TraceType = 'Tbinned3D';
else
    error('Invalid choice of trace type. Can use either "fluo", "fluo3d", "anaphasealigned", or "anaphasealigned3d".') % change to error
end

%%
if strcmpi(parameter, 'TimeOns') | strcmpi(parameter, 'TimeOn')
    parameter = 'TimeOns';
elseif strcmpi(parameter, 'TimeOffs') | strcmpi(parameter, 'TimeOff')
    parameter = 'TimeOffs';
elseif strcmpi(parameter, 'TranscriptionWindows')| strcmpi(parameter, 'TranscriptionWindow')
    parameter = 'TranscriptionWindows';
elseif strcmpi(parameter, 'ElongationTimes') | strcmpi(parameter, 'ElongationTime')
    parameter = 'ElongationTimes';
elseif strcmpi(parameter, 'ElongationRates') | strcmpi(parameter, 'ElongationRate')
    parameter = 'ElongationRates';
elseif strcmpi(parameter, 'LoadingRates') | strcmpi(parameter, 'LoadingRate')
    parameter = 'LoadingRates';
elseif strcmpi(parameter, 'PostTranscriptionDurations') | strcmpi(parameter, 'PostTranscriptionDuration')
    parameter = 'PostTranscriptionDurations';
elseif strcmpi(parameter, 'PlateauHeights') | strcmpi(parameter, 'PlateauHeight')
    parameter = 'PlateauHeights';
elseif strcmpi(parameter, 'MaxFluos') | strcmpi(parameter, 'MaxFluo')
    parameter = 'MaxFluos';
end

if ~strcmpi(parameter, 'TimeOns') & ~strcmpi(parameter, 'TranscriptionWindows') & ...
        ~strcmpi(parameter, 'ElongationTimes') & ~strcmpi(parameter, 'ElongationRates') & ...
        ~strcmpi(parameter, 'LoadingRates') & ~strcmpi(parameter, 'PostTranscriptionDurations') & ...
        ~strcmpi(parameter, 'PlateauHeights') & ~strcmpi(parameter, 'MaxFluos') 
    IncludeFits = false;
elseif strcmpi(parameter, 'TimeOns') | strcmpi(parameter, 'TranscriptionWindows') | strcmpi(parameter, 'ElongationTimes') | strcmpi(parameter, 'PostTranscriptionDurations') 
    
    LegSide = 'right';
else
    LegSide = 'left';
    
end


%%
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

Temp_obs = this.Temp_obs;
Temp_sp = this.Temp_sps;

if strcmp(lower(PlottingColors), 'default')
    [~, colors] = getColorPalettes();
elseif strcmp(lower(PlottingColors), 'pboc')
    [colors, ~] = getColorPalettes();
end

%%

R = this.R;
R2bound = this.R2bound;

temperatures = flip(unique(this.Temp_sps(this.ProcessedExperiments)));
NumTemperatures = length(temperatures);


APResolution = this.Experiments{1}.APResolution;
APbins = 0:APResolution:1;
NumAPbins = length(APbins);

NumSets = length(this.ExperimentPrefixes);
UseSet = ismember(1:NumSets, this.ProcessedExperiments);
Nsigfigs = 3;
legend_labels = this.LegendLabels;


MarkerStyles = {'o', 'd', 's', '>', '^','p', 'h', '*', 'x'};

TempMatches = cell(1, NumTemperatures);

for t_index = 1:NumTemperatures
    TempMatches{t_index} = find((this.Temp_sps == temperatures(t_index)) & UseSet);
end


%% Load relevant parameters into memory
[PlottedParams, PlottedParamSEs,R2s, ylab,OutputString,GlobalPlotYmax,GlobalPlotYmin,LogPlotYmin] = ...
    getPlottingVariables(this, parameter,  TraceType, R2bound);

% Calculate Plot Xlims



[BinnedParams, BinnedSEParams, Counts, ParamTemperatures, ParamSETemperatures] = ...
    getBinnedPlottingVariables(this, PlottedParams, PlottedParamSEs,R2s, R2bound);
%%

outdir2 = [outdir,filesep,OutputString];
if ~exist(outdir2, 'dir')
    mkdir(outdir2)
end
outdir3 = [outdir2, filesep, TraceType];
if ~exist(outdir3, 'dir')
    mkdir(outdir3)
end

outdir4 = [outdir3, filesep, datestr(now, 'yyyymmdd')];
if ~exist(outdir4, 'dir')
    mkdir(outdir4)
end





PlottedParamSEs((R2s < R2bound)| (PlottedParams > GlobalPlotYmax) | (PlottedParams < GlobalPlotYmin)) = NaN;
PlottedParams((R2s < R2bound)| (PlottedParams > GlobalPlotYmax) | (PlottedParams < GlobalPlotYmin)) = NaN;

WhereValid = squeeze(sum(~isnan(PlottedParams), 1));
WhereValidAP = sum(WhereValid, 2).';
ValidAPIndices = find(WhereValidAP);
WhereValidNC = sum(WhereValid, 1);
ValidNCIndices = find(WhereValidNC);


SubplotDims = numSubplots(length(ValidAPIndices));
SubplotDims(2) = SubplotDims(2)+1;

SubFigDims = [0.98, 0.9];% 0.9*3072/1920*SubplotDims(1)/SubplotDims(2)*1.2];
SubFigDims(2) = min([SubFigDims(2), 0.95]);
SubFigDims = round(SubFigDims, 2);

LegendSubplots = (1:SubplotDims(1))*SubplotDims(2);


LegendWidth = 0.2;
LegendXPosition = 0.8;
SubplotXBuffer = 0.025;
SubplotYBuffer = 0.075;

SubplotWidth = (LegendXPosition-(0.05-SubplotXBuffer))/(SubplotDims(2)-1);
SubplotHeight = (1-(0.15-SubplotYBuffer))/SubplotDims(1);
SubplotXPositions = 0.025+(0:(SubplotDims(2)-2))*(SubplotWidth);
SubplotYPositions = 0.0+(0:(SubplotDims(1)-1))*(SubplotHeight);



for nc_idx=1:length(ValidNCIndices)
    
    NC = ValidNCIndices(nc_idx)+8;
    % Prepare Traces for plotting
    
    NCMaxParams = NaN(1,NumSets);
    AllNCParams = NaN(NumSets, length(ValidAPIndices));
    AllNCParamSEs = NaN(NumSets, length(ValidAPIndices));
    AllNCTemperatures = NaN(NumSets, length(ValidAPIndices));
    for i=1:NumSets
        SetParams = PlottedParams(i,ValidAPIndices,NC-8);
        SetSEParams = PlottedParamSEs(i,ValidAPIndices,NC-8);
        SetTemps = ones(size(SetParams))*this.Temp_sps(i);
        
        if ~all(isnan(SetSEParams))
            IncludedBins = find(~isnan(SetParams)  & (SetParams./SetSEParams >= 1)) ;
        else
            IncludedBins = find(~isnan(SetParams)) ;
        end
        if ~isempty(IncludedBins)
            TempParams = SetParams;
            TempSEParams = SetSEParams;
            TempSEParams(TempParams > GlobalPlotYmax) = 0;
            TempParams(TempParams > GlobalPlotYmax) = NaN;
            AllNCParams(i,:) = SetParams;
            AllNCParamSEs(i,:) = SetSEParams;
            TempSEParams(isnan(TempSEParams )) = 0;
            SumParams = TempParams +TempSEParams;
            SumParams(SumParams > GlobalPlotYmax)= 0;
            NCMaxParams(i) = max(SumParams);
            AllNCTemperatures(i,:) =SetTemps;
        end
    end
    
    if all(isnan(NCMaxParams))
        continue
    end
    
    NCPlotYmax = max(NCMaxParams)*1.05;
    
    eb = cell(NumSets, length(ValidAPIndices));
    prof = cell(NumSets, length(ValidAPIndices));
    FrameProfAx = cell(1,  length(ValidAPIndices));
    hlegend2 = cell(1,  length(ValidAPIndices));
    PlottedSets = zeros(1, NumSets, 'logical');
    close all
    FrameProfFig = figure(1);
    set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.01, SubFigDims(1), SubFigDims(2)]);
    set(gcf,'color','w');
    SubplotIndex = 0;
    SubplotIndexList = zeros(1, length(ValidAPIndices));
    for l = 1:length(ValidAPIndices)
        APindex = ValidAPIndices(l);
        APbin = APbins(APindex);
        
        SubplotIndex = SubplotIndex + 1;
        if SubplotIndex == 1
            FrameProfAx{l} = subplot(SubplotDims(1), SubplotDims(2), SubplotIndex, gca);
        else
            if mod(SubplotIndex,  SubplotDims(2)) == 0
                SubplotIndex = SubplotIndex+1;
            end
            FrameProfAx{l} = subplot(SubplotDims(1), SubplotDims(2), SubplotIndex);
        end
        SubplotIndexList(l) = SubplotIndex;
        
        
        
        for idx=1:NumSets
            ColIndex = find(temperatures == this.Temp_sps(idx));
            MarkerIndex = find(TempMatches{ColIndex} == idx);
            if isempty(MarkerIndex)
                MarkerIndex = length(MarkerStyles);
            end
            if ~isnan(AllNCParamSEs(idx, l))
                NoDataCondition = isnan(AllNCParams(idx, l)) | (AllNCParams(idx,l)/AllNCParamSEs(idx, l) <  1);
            else
                NoDataCondition = isnan(AllNCParams(idx, l));
            end
            if  NoDataCondition
                
                eb{idx, l} = errorbar(AllNCTemperatures(idx, l), AllNCParams(idx, l),...
                    AllNCParamSEs(idx, l),...
                    'LineStyle', 'none');
                
                hold on
                
                set(eb{idx,l}, 'color', colors(ColIndex,:), 'LineWidth', 1);
                
                set(get(get(eb{idx,l}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                
                
                prof{idx,l} = plot(AllNCTemperatures(idx, l), AllNCParams(idx, l),  MarkerStyles{MarkerIndex},...
                    'MarkerEdgeColor', [0, 0, 0],'MarkerFaceColor', colors(ColIndex,:),...
                    'MarkerSize', 8, 'linestyle', 'none');
                
                set(eb{idx,l},'Visible','off'); %'off' or 'on'
                set(prof{idx,l},'Visible','off'); %'off' or 'on'
                
                
            else
                PlottedSets(idx) = 1;
                eb{idx, l} = errorbar(AllNCTemperatures(idx, l), AllNCParams(idx, l),...
                    AllNCParamSEs(idx, l),...
                    'LineStyle', 'none');
                
                hold on
                
                set(eb{idx,l}, 'color', colors(ColIndex,:), 'LineWidth', 1);
                
                set(get(get(eb{idx,l}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                
                
                prof{idx,l} = plot(AllNCTemperatures(idx, l), AllNCParams(idx, l),MarkerStyles{MarkerIndex},...
                    'MarkerEdgeColor', [0, 0, 0],'MarkerFaceColor', colors(ColIndex,:),...
                    'MarkerSize', 8, 'linestyle', 'none');
                
                
                
                set(eb{idx,l},'Visible','on'); %'off' or 'on'
                set(prof{idx,l},'Visible','on'); %'off' or 'on'
                
            end
        end
        
        
        grid on
        
        SubplotIndex = SubplotIndexList(l);
        ColumnIndex = mod(SubplotIndex, SubplotDims(2));
        if ColumnIndex == 0
            ColumnIndex = SubplotDims(2);
        end
        RowIndex = fix(SubplotIndex/SubplotDims(2))+1;
        %disp([num2str(RowIndex), ', ', num2str(ColumnIndex)]);
        pos = get(FrameProfAx{l}, 'position');
        pos(1) = SubplotXPositions(ColumnIndex)+SubplotXBuffer;
        pos(3) = SubplotWidth-SubplotXBuffer;
        pos(2) = SubplotYPositions(end-(RowIndex-1))+SubplotYBuffer;
        pos(4) = SubplotHeight-SubplotYBuffer;
        set(FrameProfAx{l}, 'position', pos);
        
        
        if IncludeFits
            %             Ea_placeholder = plot(1, 10, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
            [fitx,fity, ci, Ea, se_Ea, LogA, se_LogA, fitR2] = ...
                getActivationEnergyFitTraces(this, parameter, NC, APindex, TraceType);
            if ~isempty(fity) %  & fitR2 > 0.5
                lab1a = MeanSE_num2str(Ea, se_Ea, Nsigfigs);
                Ea_leglabel = ['E_{A} = ', lab1a.m,...
                    ' \pm ', lab1a.se, ' kJ/mol, R^2 = ' num2str(fitR2, 2)];
                
                
                
                if fity(1) < fity(2)
                    SubLegSide = 'left';
                else
                    SubLegSide = 'right';
                end
                dim  = [0 0 0 0];
                if strcmpi(SubLegSide, 'right')
                    dim(1) = FrameProfAx{l}.Position(1) + FrameProfAx{l}.Position(3)- 0.088;
                else
                    dim(1) = FrameProfAx{l}.Position(1) + .005;
                end
                
                dim(3) = 0.085;
                
                dim(2) = FrameProfAx{l}.Position(2) + FrameProfAx{l}.Position(4)-0.0334;
                dim(4) = 0.0286;
                
                hlegend2{l} = annotation('textbox',dim,'String',Ea_leglabel,...
                    'FitBoxToText','on', 'HorizontalAlignment', SubLegSide);
                hlegend2{l}.FontSize = 9;
                
            else
                Ea_leglabel = 'No fit info';
                SubLegSide = LegSide;
                
                
                dim  = [0 0 0 0];
                if strcmpi(SubLegSide, 'right')
                    dim(1) = FrameProfAx{l}.Position(1) + FrameProfAx{l}.Position(3)- 0.048;
                else
                    dim(1) = FrameProfAx{l}.Position(1) + .005;
                end
                
                dim(3) = 0.045;
                
                dim(2) = FrameProfAx{l}.Position(2) + FrameProfAx{l}.Position(4)-0.0334;
                dim(4) = 0.0286;
                
                hlegend2{l} = annotation('textbox',dim,'String',Ea_leglabel,...
                    'FitBoxToText','on', 'HorizontalAlignment', SubLegSide);
                hlegend2{l}.FontSize = 9;
                
                %                 hlegend2{l} = legend(FrameProfAx{l}, Ea_placeholder, {Ea_leglabel}, 'Location', legLoc,...
                %                     'FontSize', 10);
            end
            
            hlegend2{l}.BackgroundColor = [1, 1, 1];
            
        end
        
        hold off
        
        xlim([16, 29])
        
        
   
        ylim([GlobalPlotYmin, NCPlotYmax])
        
        
        FrameProfAx{l}.YAxis.FontSize = 10;
        FrameProfAx{l}.XAxis.FontSize = 10;
        
        title(['AP: ', num2str(round(APbin, 3))])
    end
    
    
    
    
    
    LegendAx = subplot(SubplotDims(1), SubplotDims(2), LegendSubplots);
    legend_profs = {};
    legend_labels3 = {};
    for SetIndex = 1:NumSets
        if ~PlottedSets(SetIndex)
            continue
        end
        ColIndex = find(temperatures == Temp_sp(SetIndex));
        MarkerIndex = find(TempMatches{ColIndex} == SetIndex);
        
        legend_profs{1, length(legend_profs)+1} = plot([0], [0], [MarkerStyles{MarkerIndex}, '-'],...
            'MarkerEdgeColor', [0, 0, 0],'MarkerFaceColor', colors(ColIndex,:), 'Color', colors(ColIndex,:),...
            'MarkerSize', 8);
        hold on
        
        legend_labels3{1, length(legend_labels3)+1} = legend_labels{SetIndex};
        
        
    end
    hold off
    axis off
    
    pos = get(LegendAx, 'Position');
    pos(1) = LegendXPosition;
    pos(3) = LegendWidth;
    set(LegendAx, 'Position', pos);
    
    hlegend = legend(LegendAx,  legend_labels3,...
        'FontSize',10, 'Location', 'east');
    %     hlegend.Position(1) = LegendAx.Position(1)+ LegendAx.Position(3)/2 - hlegend.Position(3)/2;
    %     hlegend.Position(2) = LegendAx.Position(2)+ LegendAx.Position(4)/2 - hlegend.Position(4)/2;
    hlegend.Title.String = 'Data Sets';
    hlegend.Title.FontSize = 10;
    
    for l = 1:length(ValidAPIndices)
        SubplotIndex = SubplotIndexList(l);
        ColumnIndex = mod(SubplotIndex, SubplotDims(2));
        if ColumnIndex == 0
            ColumnIndex = SubplotDims(2);
        end
        RowIndex = fix(SubplotIndex/SubplotDims(2))+1;
        %disp([num2str(RowIndex), ', ', num2str(ColumnIndex)]);
        pos = get(FrameProfAx{l}, 'position');
        pos(1) = SubplotXPositions(ColumnIndex)+SubplotXBuffer;
        pos(3) = SubplotWidth-SubplotXBuffer;
        pos(2) = SubplotYPositions(end-(RowIndex-1))+SubplotYBuffer;
        pos(4) = SubplotHeight-SubplotYBuffer;
        set(FrameProfAx{l}, 'position', pos);
    end
    
    
    
    
    ax = axes(FrameProfFig);
    % Specifiy the active side for the axes ax
    yyaxis(ax, 'left');
    % Specify visibility of the current axis as 'off'
    han = gca;
    han.Visible = 'off';
    han.XLabel.Visible = 'on';
    han.YLabel.Visible = 'on';
    han.Title.Visible = 'on';
    xlabel('Temperature (°C)')
    
    
    ylabel(ylab)
    if exist('PlotTitle', 'var')
        title({PlotTitle, ['Nuclear Cycle ',num2str(NC)]})
    else
        title(['Nuclear Cycle ',num2str(NC)])
    end
    
    
    
    han.XLabel.Color = [0, 0, 0];
    han.XLabel.FontSize = 14;
    xlabel_pos = han.XLabel.Position;
    xlabel_pos(2) = 0.025-han.Position(2);
    xlabel_pos(1) = LegendXPosition/2;
    han.XLabel.Position = xlabel_pos;
    
    han.YLabel.Color = [0, 0, 0];
    han.YLabel.FontSize = 14;
    ylabel_pos = han.YLabel.Position;
    ylabel_pos(1) = -0.025-han.Position(1);
    han.YLabel.Position = ylabel_pos;
    
    
    han.Title.Color = [0, 0, 0];
    han.Title.FontSize = 14;
    title_pos = han.Title.Position;
    title_pos(2) = .96/(han.Position(2)+han.Position(4));
    title_pos(1) = LegendXPosition/2;
    han.Title.Position = title_pos;
    han.Title.FontWeight = 'normal';
    %try
    %
    saveas(FrameProfFig,[outdir4, filesep,...
        'Subplots', OutputString,  '_NC',num2str(NC),'.png']);
end
close all

%%

function PlotLTMLogSubplotsTrapParamsVsTemp(this, parameter, outdir, varargin)
%%

% PlotTitle, PlottingColors, UseDifferentColors,
% UseDiffProfiles, UsePhysicalAPLength
IncludeFits = true;
x = 1;
while x <= length(varargin)
    if strcmpi(varargin{x}, 'plottitle')
        PlotTitle = varargin{x+1};
        x = x+1;
    elseif strcmpi(varargin{x}, 'plottingcolors')
        PlottingColors = varargin{x+1};
        x = x+1;
    elseif strcmpi(varargin{x}, 'tracetype')
        TraceType = lower(varargin{x+1});
        x = x+1;
    elseif strcmpi(varargin{x}, 'ExcludeFits')
        IncludeFits = false;
    end
    x = x+1;
end


if ~exist('PlottingColors', 'var')
    PlottingColors = 'default';
elseif ~strcmp(lower(PlottingColors), 'default')  & ~strcmp(lower(PlottingColors), 'pboc')
    error('Invalid choice of plotting colors. Can use either "default" or "pboc".') % change to error
end
if ~exist('TraceType', 'var')
    TraceType = 'AnaphaseAligned';
elseif strcmpi(TraceType, 'Fluo3D') | strcmpi(TraceType, 'Unaligned3D')
    TraceType = 'Unaligned3D';
elseif strcmpi(TraceType, 'Fluo')| strcmpi(TraceType, 'Unaligned')
    TraceType = 'Unaligned';
elseif strcmpi(TraceType, 'AnaphaseAligned')
    TraceType = 'AnaphaseAligned';
elseif strcmpi(TraceType, 'AnaphaseAligned3D')
    TraceType = 'AnaphaseAligned3D';
elseif strcmpi(TraceType, 'Tbinned')
    TraceType = 'Tbinned';
elseif strcmpi(TraceType, 'Tbinned3D')
    TraceType = 'Tbinned3D';
else
    error('Invalid choice of trace type. Can use either "fluo", "fluo3d", "anaphasealigned", or "anaphasealigned3d".') % change to error
end

%%

if strcmpi(parameter, 'TimeOns') | strcmpi(parameter, 'TimeOn')
    parameter = 'TimeOns';
elseif strcmpi(parameter, 'TimeOffs') | strcmpi(parameter, 'TimeOff')
    parameter = 'TimeOffs';
elseif strcmpi(parameter, 'TranscriptionWindows')| strcmpi(parameter, 'TranscriptionWindow')
    parameter = 'TranscriptionWindows';
elseif strcmpi(parameter, 'ElongationTimes') | strcmpi(parameter, 'ElongationTime')
    parameter = 'ElongationTimes';
elseif strcmpi(parameter, 'ElongationRates') | strcmpi(parameter, 'ElongationRate')
    parameter = 'ElongationRates';
elseif strcmpi(parameter, 'LoadingRates') | strcmpi(parameter, 'LoadingRate')
    parameter = 'LoadingRates';
elseif strcmpi(parameter, 'PostTranscriptionDurations') | strcmpi(parameter, 'PostTranscriptionDuration')
    parameter = 'PostTranscriptionDurations';
elseif strcmpi(parameter, 'PlateauHeights') | strcmpi(parameter, 'PlateauHeight')
    parameter = 'PlateauHeights';
elseif strcmpi(parameter, 'MaxFluos') | strcmpi(parameter, 'MaxFluo')
    parameter = 'MaxFluos';
end

if ~strcmpi(parameter, 'TimeOns') & ~strcmpi(parameter, 'TranscriptionWindows') & ...
        ~strcmpi(parameter, 'ElongationTimes') & ~strcmpi(parameter, 'ElongationRates') & ...
        ~strcmpi(parameter, 'LoadingRates') & ~strcmpi(parameter, 'PostTranscriptionDurations')  & ...
        ~strcmpi(parameter, 'PlateauHeights') & ~strcmpi(parameter, 'MaxFluos') 
    IncludeFits = false;
elseif strcmpi(parameter, 'TimeOns') | strcmpi(parameter, 'TranscriptionWindows') | strcmpi(parameter, 'ElongationTimes') | strcmpi(parameter, 'PostTranscriptionDurations') 

    legLoc = 'northwest';
else
    legLoc = 'northeast';
    
end



%%
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

Temp_obs = this.Temp_obs;
Temp_sp = this.Temp_sps;

if strcmp(lower(PlottingColors), 'default')
    [~, colors] = getColorPalettes();
elseif strcmp(lower(PlottingColors), 'pboc')
    [colors, ~] = getColorPalettes();
end

%%

R = this.R;
R2bound = this.R2bound;

temperatures = flip(unique(this.Temp_sps(this.ProcessedExperiments)));
NumTemperatures = length(temperatures);


APResolution = this.Experiments{1}.APResolution;
APbins = 0:APResolution:1;
NumAPbins = length(APbins);

NumSets = length(this.ExperimentPrefixes);
UseSet = ismember(1:NumSets, this.ProcessedExperiments);
Nsigfigs = 3;
legend_labels = this.LegendLabels;


MarkerStyles = {'o', 'd', 's', '>', '^','p', 'h', '*', 'x'};

TempMatches = cell(1, NumTemperatures);

for t_index = 1:NumTemperatures
    TempMatches{t_index} = find((this.Temp_sps == temperatures(t_index)) & UseSet);
end


%% Load relevant parameters into memory
[PlottedParams, PlottedParamSEs,R2s, ylab,OutputString,GlobalPlotYmax,GlobalPlotYmin,LogPlotYmin] = ...
    getPlottingVariables(this, parameter,  TraceType, R2bound);

% Calculate Plot Xlims

TemperatureVector = 1./(R*(temperatures + 273));
Subplot2DataXmin = min(TemperatureVector);
Subplot2DataXmax = max(TemperatureVector);
Subplot2Xspan = Subplot2DataXmax-Subplot2DataXmin;
Plot2Xmin = Subplot2DataXmin - Subplot2Xspan*.05;
Plot2Xmax = Subplot2DataXmax + Subplot2Xspan*.05;



[BinnedParams, BinnedSEParams, Counts, ParamTemperatures, ParamSETemperatures] = ...
    getBinnedPlottingVariables(this, PlottedParams, PlottedParamSEs,R2s, R2bound);
%%

outdir2 = [outdir,filesep,OutputString];
if ~exist(outdir2, 'dir')
    mkdir(outdir2)
end
outdir3 = [outdir2, filesep, TraceType];
if ~exist(outdir3, 'dir')
    mkdir(outdir3)
end
outdir4 = [outdir3, filesep, datestr(now, 'yyyymmdd')];
if ~exist(outdir4, 'dir')
    mkdir(outdir4)
end





PlottedParamSEs((R2s < R2bound)| (PlottedParams > GlobalPlotYmax) | (PlottedParams < GlobalPlotYmin)) = NaN;
PlottedParams((R2s < R2bound)| (PlottedParams > GlobalPlotYmax) | (PlottedParams < GlobalPlotYmin)) = NaN;

WhereValid = squeeze(sum(~isnan(PlottedParams), 1));
WhereValidAP = sum(WhereValid, 2).';
ValidAPIndices = find(WhereValidAP);
WhereValidNC = sum(WhereValid, 1);
ValidNCIndices = find(WhereValidNC);


SubplotDims = numSubplots(length(ValidAPIndices));
SubplotDims(2) = SubplotDims(2)+1;

SubFigDims = [0.98, 0.9];% 0.9*3072/1920*SubplotDims(1)/SubplotDims(2)*1.2];
SubFigDims(2) = min([SubFigDims(2), 0.95]);
SubFigDims = round(SubFigDims, 2);

LegendSubplots = (1:SubplotDims(1))*SubplotDims(2);


LegendWidth = 0.2;
LegendXPosition = 0.8;
SubplotXBuffer = 0.025;
SubplotYBuffer = 0.075;

SubplotWidth = (LegendXPosition-(0.05-SubplotXBuffer))/(SubplotDims(2)-1);
SubplotHeight = (1-(0.15-SubplotYBuffer))/SubplotDims(1);
SubplotXPositions = 0.025+(0:(SubplotDims(2)-2))*(SubplotWidth);
SubplotYPositions = 0.0+(0:(SubplotDims(1)-1))*(SubplotHeight);



for nc_idx=1:length(ValidNCIndices)
    
    NC = ValidNCIndices(nc_idx)+8;
    % Prepare Traces for plotting
    
    NCMaxParams = NaN(1,NumSets);
    AllNCParams = NaN(NumSets, length(ValidAPIndices));
    AllNCParamSEs = NaN(NumSets, length(ValidAPIndices));
    AllNCTemperatures = NaN(NumSets, length(ValidAPIndices));
    for i=1:NumSets
        SetParams = PlottedParams(i,ValidAPIndices,NC-8);
        SetSEParams = PlottedParamSEs(i,ValidAPIndices,NC-8);
        SetTemps = ones(size(SetParams))*this.Temp_sps(i);
        
        if ~all(isnan(SetSEParams))
            IncludedBins = find(~isnan(SetParams)  & (SetParams./SetSEParams >= 1)) ;
        else
            IncludedBins = find(~isnan(SetParams)) ;
        end
        if ~isempty(IncludedBins)
            TempParams = SetParams;
            TempSEParams = SetSEParams;
            TempSEParams(TempParams > GlobalPlotYmax) = 0;
            TempParams(TempParams > GlobalPlotYmax) = NaN;
            AllNCParams(i,:) = SetParams;
            AllNCParamSEs(i,:) = SetSEParams;
            TempSEParams(isnan(TempSEParams )) = 0;
            SumParams = TempParams +TempSEParams;
            SumParams(SumParams > GlobalPlotYmax)= 0;
            NCMaxParams(i) = max(SumParams);
            AllNCTemperatures(i,:) =SetTemps;
        end
    end
    
    if all(isnan(NCMaxParams))
        continue
    end
    
    NCPlotYmax = max(NCMaxParams)*1.5;
    
    ci_plotline = cell(1,  length(ValidAPIndices));
    fitted_prof = cell(1,  length(ValidAPIndices));
    eb = cell(NumSets, length(ValidAPIndices));
    prof = cell(NumSets, length(ValidAPIndices));
    FrameProfAx = cell(1,  length(ValidAPIndices));
    hlegend2 = cell(1,  length(ValidAPIndices));
    icons = cell(1,  length(ValidAPIndices));
    
    PlottedSets = zeros(1, NumSets, 'logical');
    close all
    FrameProfFig = figure(1);
    set(FrameProfFig,'units', 'normalized', 'position',[0.01, 0.01, SubFigDims(1), SubFigDims(2)]);
    set(gcf,'color','w');
    SubplotIndex = 0;
    SubplotIndexList = zeros(1, length(ValidAPIndices));
    for l = 1:length(ValidAPIndices)
        APindex = ValidAPIndices(l);
        APbin = APbins(APindex);
        
        SubplotIndex = SubplotIndex + 1;
        if SubplotIndex == 1
            FrameProfAx{l} = subplot(SubplotDims(1), SubplotDims(2), SubplotIndex, gca);
        else
            if mod(SubplotIndex,  SubplotDims(2)) == 0
                SubplotIndex = SubplotIndex+1;
            end
            FrameProfAx{l} = subplot(SubplotDims(1), SubplotDims(2), SubplotIndex);
        end
        SubplotIndexList(l) = SubplotIndex;
        
        
        if IncludeFits
            [fitx,fity, ci, Ea, se_Ea, LogA, se_LogA, fitR2] = ...
                getActivationEnergyFitTraces(this, parameter, NC, APindex, TraceType);
            if ~isempty(fity)
                curve1 = ci(:,1).';
                curve2 = ci(:,2).';
                curve1(curve1 > max(GlobalPlotYmax*1.05,1)) = max(GlobalPlotYmax*1.05,1);
                curve1(curve1 < LogPlotYmin) = LogPlotYmin;
                curve2(curve2 > max(GlobalPlotYmax*1.05,1)) = max(GlobalPlotYmax*1.05,1);
                curve2(curve2 < LogPlotYmin) = LogPlotYmin;
                t_vector2 = [fitx, fliplr(fitx)];
                inBetween = [curve1, fliplr(curve2)];
                ci_plotline{l} = fill(t_vector2, inBetween, [0, 0, 0]);
                ci_plotline{l}.FaceAlpha = 0.2;
                set(ci_plotline{l}, 'EdgeColor', 'none');
                hold on
                
                fitted_prof{l} = plot(fitx, fity,  '-', 'Color', [0, 0, 0]);
                
                
            else
                ci_plotline{l} = fill([0, -1, -1, 0], [0, 0, max(GlobalPlotYmax*1.2,1), max(GlobalPlotYmax*1.2,1)], [0, 0, 0] );
                ci_plotline{l}.FaceAlpha = 0.2;
                set(ci_plotline{l}, 'EdgeColor', 'none');
                hold on
                fitted_prof{l} = plot([-1], [0], '-', 'Color', [0, 0, 0]);
                set(fitted_prof{l},'Visible','off'); %'off' or 'on'
                set(ci_plotline{l},'Visible','off'); %'off' or 'on'
                
                
            end
        end
        
        for idx=1:NumSets
            ColIndex = find(temperatures == this.Temp_sps(idx));
            MarkerIndex = find(TempMatches{ColIndex} == idx);
            if isempty(MarkerIndex)
                MarkerIndex = length(MarkerStyles);
            end
            if ~isnan(AllNCParamSEs(idx, l))
                NoDataCondition = isnan(AllNCParams(idx, l)) | (AllNCParams(idx,l)/AllNCParamSEs(idx, l) <  1);
            else
                NoDataCondition = isnan(AllNCParams(idx, l));
            end
            if  NoDataCondition
                
                eb{idx, l} = errorbar(1/(R*(AllNCTemperatures(idx, l)+273)), AllNCParams(idx, l),...
                    AllNCParamSEs(idx, l),...
                    'LineStyle', 'none');
                
                hold on
                
                set(eb{idx,l}, 'color', colors(ColIndex,:), 'LineWidth', 1);
                
                set(get(get(eb{idx,l}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                
                
                prof{idx,l} = plot(1/(R*(AllNCTemperatures(idx, l)+273)), AllNCParams(idx, l),  MarkerStyles{MarkerIndex},...
                    'MarkerEdgeColor', [0, 0, 0],'MarkerFaceColor', colors(ColIndex,:),...
                    'MarkerSize', 8, 'linestyle', 'none');
                
                set(eb{idx,l},'Visible','off'); %'off' or 'on'
                set(prof{idx,l},'Visible','off'); %'off' or 'on'
                
                
            else
                PlottedSets(idx) = 1;
                eb{idx, l} = errorbar(1/(R*(AllNCTemperatures(idx, l)+273)), AllNCParams(idx, l),...
                    AllNCParamSEs(idx, l),...
                    'LineStyle', 'none');
                
                hold on
                
                set(eb{idx,l}, 'color', colors(ColIndex,:), 'LineWidth', 1);
                
                set(get(get(eb{idx,l}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
                
                
                prof{idx,l} = plot(1/(R*(AllNCTemperatures(idx, l)+273)), AllNCParams(idx, l),MarkerStyles{MarkerIndex},...
                    'MarkerEdgeColor', [0, 0, 0],'MarkerFaceColor', colors(ColIndex,:),...
                    'MarkerSize', 8, 'linestyle', 'none');
                
                
                
                set(eb{idx,l},'Visible','on'); %'off' or 'on'
                set(prof{idx,l},'Visible','on'); %'off' or 'on'
                
            end
        end
        
        
        grid on
        
        SubplotIndex = SubplotIndexList(l);
        ColumnIndex = mod(SubplotIndex, SubplotDims(2));
        if ColumnIndex == 0
            ColumnIndex = SubplotDims(2);
        end
        RowIndex = fix(SubplotIndex/SubplotDims(2))+1;
        %disp([num2str(RowIndex), ', ', num2str(ColumnIndex)]);
        pos = get(FrameProfAx{l}, 'position');
        pos(1) = SubplotXPositions(ColumnIndex)+SubplotXBuffer;
        pos(3) = SubplotWidth-SubplotXBuffer;
        pos(2) = SubplotYPositions(end-(RowIndex-1))+SubplotYBuffer;
        pos(4) = SubplotHeight-SubplotYBuffer;
        set(FrameProfAx{l}, 'position', pos);
        
        
        if IncludeFits
            %             Ea_placeholder = plot(1, 10, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
            [fitx,fity, ci, Ea, se_Ea, LogA, se_LogA, fitR2] = ...
                getActivationEnergyFitTraces(this, parameter, NC, APindex, TraceType);
            

                
                
                
                
              
            if ~isempty(fity) 
                lab1a = MeanSE_num2str(Ea, se_Ea, Nsigfigs);
          
                subplot_legend_labels = cell(1, 2);
                
                subplot_legend_labels{1} = ['E_{A} = ', lab1a.m,...
                    ' \pm ', lab1a.se, ' kJ/mol'];
                subplot_legend_labels{2} = ['95% Conf., R^2: ', num2str(fitR2, 2)];
                
                
                if fity(1) < fity(2)
                    legLoc = 'northeast';
                else
                    legLoc = 'northwest';
                end
                
                [hlegend2{l}, icons{l}] = legend(FrameProfAx{l}, [fitted_prof{l}, ci_plotline{l}], subplot_legend_labels, 'Location', legLoc,...
                    'FontSize', 9);
                
                 OldStrPos = icons{l}(1).Position(1);
                 icons{l}(5).FaceAlpha = 0.2;
                 icons{l}(5).Vertices(3,1) = .12;
                 icons{l}(5).Vertices(4,1) = .12;
                 icons{l}(5).Vertices(1,2) = .05;
                 icons{l}(5).Vertices(4,2) = .05;
                 icons{l}(5).Vertices(2,2) = .3;
                 icons{l}(5).Vertices(3,2) = .3;
                 icons{l}(3).XData(2) = .12;
                 icons{l}(4).XData(2) = .12;
                 icons{l}(1).Position(1) = .17;
                 icons{l}(2).Position(1) = .17;
%                  LegPos = get(hlegend2{l}, 'position');
%                  LegPos(3) = LegPos(3) *(1-(OldStrPos-icons{l}(1).Position(1)+.1));
                 set(hlegend2{l}, 'Location', legLoc);
  
                icons{1}(1).FontSize  = 10;
                icons{1}(2).FontSize  = 10;
  
                 
                
%                 h1 = plot(1:10);
%                 [hh,icons,plots,txt] = legend({'Line 1'});
%                 p1 = icons(1).Position;
%                 icons(1).Position = [0.3 p1(2) 0];
%                 icons(2).XData = [0.05 0.2];
                
            else
                subplot_legend_labels = cell(1, 2);
                subplot_legend_labels{2} = '';
                subplot_legend_labels{1} = 'No Fit Info';
                [hlegend2{l}, icons{l}]  = legend(FrameProfAx{l}, [fitted_prof{l}, ci_plotline{l}], subplot_legend_labels, 'Location', legLoc,...
                    'FontSize', 9);
                
                
                OldStrPos = icons{l}(1).Position(1);
                 icons{l}(5).FaceAlpha = 0.2;
                 icons{l}(5).Vertices(3,1) = .12;
                 icons{l}(5).Vertices(4,1) = .12;
                 icons{l}(5).Vertices(1,2) = .05;
                 icons{l}(5).Vertices(4,2) = .05;
                 icons{l}(5).Vertices(2,2) = .3;
                 icons{l}(5).Vertices(3,2) = .3;
                 icons{l}(3).XData(2) = .12;
                 icons{l}(4).XData(2) = .12;
                 icons{l}(1).Position(1) = .17;
                 icons{l}(2).Position(1) = .17;
%                  LegPos = get(hlegend2{l}, 'pos
                set(hlegend2{l}, 'Location', legLoc);
                icons{1}(1).FontSize  = 10;
                icons{1}(2).FontSize  = 10;

            end
            
            
        end
        
        hold off
        
        xlim([Plot2Xmin, Plot2Xmax])
        
        
        
        ylim([LogPlotYmin, 1.5*NCPlotYmax])
        
        
        FrameProfAx{l}.YAxis.FontSize = 10;
        FrameProfAx{l}.XAxis.FontSize = 10;
        
        set(FrameProfAx{l}, 'YScale', 'log')
        
        title(['AP: ', num2str(round(APbin, 3))])
    end
    
    
    
    
    
    LegendAx = subplot(SubplotDims(1), SubplotDims(2), LegendSubplots);
    legend_profs = {};
    legend_labels3 = {};
    for SetIndex = 1:NumSets
        if ~PlottedSets(SetIndex)
            continue
        end
        ColIndex = find(temperatures == Temp_sp(SetIndex));
        MarkerIndex = find(TempMatches{ColIndex} == SetIndex);
        
        legend_profs{1, length(legend_profs)+1} = plot([0], [0], [MarkerStyles{MarkerIndex}, '-'],...
            'MarkerEdgeColor', [0, 0, 0],'MarkerFaceColor', colors(ColIndex,:), 'Color', colors(ColIndex,:),...
            'MarkerSize', 8);
        hold on
        
        legend_labels3{1, length(legend_labels3)+1} = legend_labels{SetIndex};
        
        
    end
    hold off
    axis off
    
    pos = get(LegendAx, 'Position');
    pos(1) = LegendXPosition;
    pos(3) = LegendWidth;
    set(LegendAx, 'Position', pos);
    
    hlegend = legend(LegendAx,  legend_labels3,...
        'FontSize',10, 'Location', 'east');
    %     hlegend.Position(1) = LegendAx.Position(1)+ LegendAx.Position(3)/2 - hlegend.Position(3)/2;
    %     hlegend.Position(2) = LegendAx.Position(2)+ LegendAx.Position(4)/2 - hlegend.Position(4)/2;
    hlegend.Title.String = 'Data Sets';
    hlegend.Title.FontSize = 10;
    
%     for l = 1:length(ValidAPIndices)
%         SubplotIndex = SubplotIndexList(l);
%         ColumnIndex = mod(SubplotIndex, SubplotDims(2));
%         if ColumnIndex == 0
%             ColumnIndex = SubplotDims(2);
%         end
%         RowIndex = fix(SubplotIndex/SubplotDims(2))+1;
%         %disp([num2str(RowIndex), ', ', num2str(ColumnIndex)]);
%         pos = get(FrameProfAx{l}, 'position');
%         pos(1) = SubplotXPositions(ColumnIndex)+SubplotXBuffer;
%         pos(3) = SubplotWidth-SubplotXBuffer;
%         pos(2) = SubplotYPositions(end-(RowIndex-1))+SubplotYBuffer;
%         pos(4) = SubplotHeight-SubplotYBuffer;
%         set(FrameProfAx{l}, 'position', pos);
%     end
    
    
    
    
    ax = axes(FrameProfFig);
    % Specifiy the active side for the axes ax
    yyaxis(ax, 'left');
    % Specify visibility of the current axis as 'off'
    han = gca;
    han.Visible = 'off';
    han.XLabel.Visible = 'on';
    han.YLabel.Visible = 'on';
    han.Title.Visible = 'on';
    xlabel('1/(RT) (mol/kJ)')
    
    
    ylabel(ylab)
    if exist('PlotTitle', 'var')
        title({PlotTitle, ['Nuclear Cycle ',num2str(NC)]})
    else
        title(['Nuclear Cycle ',num2str(NC)])
    end
    
    
    
    han.XLabel.Color = [0, 0, 0];
    han.XLabel.FontSize = 14;
    xlabel_pos = han.XLabel.Position;
    xlabel_pos(2) = 0.025-han.Position(2);
    xlabel_pos(1) = LegendXPosition/2;
    han.XLabel.Position = xlabel_pos;
    
    han.YLabel.Color = [0, 0, 0];
    han.YLabel.FontSize = 14;
    ylabel_pos = han.YLabel.Position;
    ylabel_pos(1) = -0.01-han.Position(1);
    han.YLabel.Position = ylabel_pos;
    
    
    han.Title.Color = [0, 0, 0];
    han.Title.FontSize = 14;
    title_pos = han.Title.Position;
    title_pos(2) = .96/(han.Position(2)+han.Position(4));
    title_pos(1) = LegendXPosition/2;
    han.Title.Position = title_pos;
    han.Title.FontWeight = 'normal';
    %try
    %
    saveas(FrameProfFig,[outdir4, filesep,...
        'LogSubplots', OutputString,  '_NC',num2str(NC),'.png']);
end
close all