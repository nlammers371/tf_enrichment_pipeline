function GenerateProjectSubplots(ResultsPaths, outdir)
close all




if ~exist(outdir, 'dir')
    mkdir(outdir)
end

R = 8.314*10^(-3);
UseSharedYAxis = false;
UseSharedYAxisFrequency = false;
UseSharedYAxisInit = true;

APResolution = 0.025;
APbins = (0:APResolution:1)*100;

CompiledParameters = CompileParameters(ResultsPaths);
NumSets = size(CompiledParameters.InitiationRates, 1);
IncludedAPbins =  find((sum(squeeze(sum(~isnan(CompiledParameters.InitiationRates), 1)), 1)) > 0);
MinAPbin = min(IncludedAPbins);
MaxAPbin = max(IncludedAPbins);

NumTimeBins  =  size(CompiledParameters.InitiationRates, 2);
IncludedTimeBins =  find((sum(squeeze(sum(~isnan(CompiledParameters.InitiationRates), 1)), 2)) > 0).';

InitYmin = nanmin(nanmin(nanmin(CompiledParameters.InitiationRates+CompiledParameters.InitiationRatesStdErr)));
InitYmax = nanmax(nanmax(nanmax(CompiledParameters.InitiationRates+CompiledParameters.InitiationRatesStdErr)));

InitYmin = max(floor(InitYmin/50)*50, 0);
InitYmax = ceil(InitYmax/50)*50;

DurYmin = nanmin(nanmin(nanmin(CompiledParameters.Durations+CompiledParameters.DurationsStdErr)));
DurYmax = nanmax(nanmax(nanmax(CompiledParameters.Durations+CompiledParameters.DurationsStdErr)));

DurYmin = max(floor(DurYmin/2)*2, 0);
DurYmax = min(ceil(DurYmax/2)*2, 20);


FreqYmin = nanmin(nanmin(nanmin(CompiledParameters.Frequencies+CompiledParameters.FrequenciesStdErr)));
FreqYmax = nanmax(nanmax(nanmax(CompiledParameters.Frequencies+CompiledParameters.FrequenciesStdErr)));

FreqYmin = max(floor(FreqYmin/0.1)*0.1, 0);
FreqYmax = min(ceil(FreqYmax/0.1)*0.1, 4);


SubplotDims = [2, 3];
SubplotDims(2) = SubplotDims(2)+1;

SubFigDims = [0.85, 0.8];% 0.9*3072/1920*SubplotDims(1)/SubplotDims(2)*1.2];
SubFigDims(2) = min([SubFigDims(2), 0.95]);
SubFigDims = round(SubFigDims, 2);

LegendSubplots = (1:SubplotDims(1))*SubplotDims(2);


LegendWidth = 0.2;
LegendXPosition = 0.78;
SubplotXBuffer = 0.05;
SubplotYBuffer = 0.1;

SubplotWidth = .23; %(LegendXPosition-(0.05-SubplotXBuffer))/(SubplotDims(2)-1);
SubplotHeight = .4;  %(1-(0.15-SubplotYBuffer))/SubplotDims(1);
SubplotXPositions = 0.025+(0:(SubplotDims(2)-2))*(SubplotWidth);
SubplotYPositions = 0.0+(0:(SubplotDims(1)-1))*(SubplotHeight);

SubplotYPositions = [.050, .5];


%%
DurFigAx = cell(1, NumSets);

dur_fig = figure(1);

set(dur_fig,'units', 'normalized', 'position',[0.01, 0.01, SubFigDims(1), SubFigDims(2)]);
set(gcf,'color','w');
if max(IncludedTimeBins) > 1
cmap = flipud(brewermap(max(IncludedTimeBins),'Spectral'));
else
    cmap = flipud(brewermap(2,'Spectral'));
end
colormap(cmap);
hold on
SubplotIndex = 0;
SubplotIndexList = zeros(1, NumSets);
for SetIndex = 1:NumSets
    SubplotIndex = SubplotIndex + 1;
    if SubplotIndex == 1
        DurFigAx{SetIndex} = subplot(SubplotDims(1), SubplotDims(2), SubplotIndex, gca);
    else
        if mod(SubplotIndex,  SubplotDims(2)) == 0
            SubplotIndex = SubplotIndex+1;
        end
        DurFigAx{SetIndex} = subplot(SubplotDims(1), SubplotDims(2), SubplotIndex);
    end
    SubplotIndexList(SetIndex) = SubplotIndex;
    PlottedTimes = zeros(1, length(CompiledParameters.TimeVector), 'logical');
    
    y_vals = squeeze(CompiledParameters.Durations(SetIndex,:,:));
    y_val_errs = squeeze(CompiledParameters.DurationsStdErr(SetIndex,:,:));
    for t = 1:length(CompiledParameters.TimeVector)
        if all(isnan(y_vals(t,:)))
            continue
        end
        set_ids = ~isnan(y_vals(t,:));
        
        x = APbins(set_ids);
        y = y_vals(t, set_ids);
        y_err = y_val_errs(t,set_ids);
        
        errorbar(DurFigAx{SetIndex}, x,y,y_err,'Color','k','Capsize',0)
        hold on
        scatter(DurFigAx{SetIndex}, x,y,'MarkerFaceColor',cmap(t,:),'MarkerEdgeColor','k')
        
    end
    
    
    grid on
    
    SubplotIndex = SubplotIndexList(SetIndex);
    ColumnIndex = mod(SubplotIndex, SubplotDims(2));
    if ColumnIndex == 0
        ColumnIndex = SubplotDims(2);
    end
    RowIndex = fix(SubplotIndex/SubplotDims(2))+1;
    %disp([num2str(RowIndex), ', ', num2str(ColumnIndex)]);
    pos = get(DurFigAx{SetIndex}, 'position');
    pos(1) = SubplotXPositions(ColumnIndex)+SubplotXBuffer;
    pos(3) = SubplotWidth-SubplotXBuffer;
    pos(2) = SubplotYPositions(end-(RowIndex-1))+SubplotYBuffer;
    pos(4) = SubplotHeight-SubplotYBuffer;
    set(DurFigAx{SetIndex}, 'position', pos);
    
    xlabel('AP position')
    ylabel('burst duration (min)')
    set(DurFigAx{SetIndex} ,'Fontsize',14)
    xlim([(MinAPbin-2)*APResolution*100 (MaxAPbin)*APResolution*100])
    %xlim([(MinAPbin-2)*APResolution*100 (MaxAPbin)*APResolution*100])
     title([CompiledParameters.ReporterLabels{SetIndex}, ' ', num2str(CompiledParameters.SetTemperatures(SetIndex)), '°C w', num2str(CompiledParameters.nSteps(SetIndex))]);
    newylim = get(DurFigAx{SetIndex}, 'ylim');
    newylim(1) = DurYmin;
    if newylim(2) > 15
        newylim(2)  = 15;
    end
    set(DurFigAx{SetIndex}, 'ylim', newylim)
    %     ylim([DurYmin, DurYmax])
    
    
end

LegendAx = subplot(SubplotDims(1), SubplotDims(2), LegendSubplots);
hold on
cmap = flipud(brewermap(max(IncludedTimeBins),'Spectral'));
colormap(cmap);
if max(IncludedTimeBins) > 1
set(LegendAx, 'Clim', [1, max(IncludedTimeBins)])
% % AxesH = axes('CLim', [-12, 12]);
% cbh = colorbar('peer', AxesH, 'h', ...
%                'XTickLabel',{'-12','-9','-6','-3','0','3','6','9','12'}, ...
%                'XTick', -12:3:12)

h = colorbar;
h.Ticks = 1:max(IncludedTimeBins);
h.TickLabels = CompiledParameters.TimeVector(IncludedTimeBins);
ylabel(h,'time cohort (minutes into nc14)')
end
hold off
set(LegendAx,'Fontsize',16)
axis off

pos = get(LegendAx, 'Position');
pos(1) = .6;
pos(3) = .2;
set(LegendAx, 'Position', pos);


for SetIndex = 1:NumSets
    SubplotIndex = SubplotIndexList(SetIndex);
    ColumnIndex = mod(SubplotIndex, SubplotDims(2));
    if ColumnIndex == 0
        ColumnIndex = SubplotDims(2);
    end
    RowIndex = fix(SubplotIndex/SubplotDims(2))+1;
    %disp([num2str(RowIndex), ', ', num2str(ColumnIndex)]);
    pos = get(DurFigAx{SetIndex}, 'position');
    pos(1) = SubplotXPositions(ColumnIndex)+SubplotXBuffer;
    pos(3) = SubplotWidth-SubplotXBuffer;
    pos(2) = SubplotYPositions(end-(RowIndex-1))+SubplotYBuffer;
    pos(4) = SubplotHeight-SubplotYBuffer;
    set(DurFigAx{SetIndex}, 'position', pos);
end

saveas(dur_fig,[outdir, 'Projectsubplots_burst_duration.png'])

%%


FreqFigAx = cell(1, NumSets);

freq_fig = figure(2);

set(freq_fig,'units', 'normalized', 'position',[0.01, 0.01, SubFigDims(1), SubFigDims(2)]);
set(gcf,'color','w');
if max(IncludedTimeBins) > 1
cmap = flipud(brewermap(max(IncludedTimeBins),'Spectral'));
else
    cmap = flipud(brewermap(2,'Spectral'));
end
colormap(cmap);
hold on
SubplotIndex = 0;
SubplotIndexList = zeros(1, NumSets);
for SetIndex = 1:NumSets
    SubplotIndex = SubplotIndex + 1;
    if SubplotIndex == 1
        FreqFigAx{SetIndex} = subplot(SubplotDims(1), SubplotDims(2), SubplotIndex, gca);
    else
        if mod(SubplotIndex,  SubplotDims(2)) == 0
            SubplotIndex = SubplotIndex+1;
        end
        FreqFigAx{SetIndex} = subplot(SubplotDims(1), SubplotDims(2), SubplotIndex);
    end
    SubplotIndexList(SetIndex) = SubplotIndex;
    PlottedTimes = zeros(1, length(CompiledParameters.TimeVector), 'logical');
    
    y_vals = squeeze(CompiledParameters.Frequencies(SetIndex,:,:));
    y_val_errs = squeeze(CompiledParameters.FrequenciesStdErr(SetIndex,:,:));
    for t = 1:length(CompiledParameters.TimeVector)
        if all(isnan(y_vals(t,:)))
            continue
        end
        set_ids = ~isnan(y_vals(t,:));
        
        x = APbins(set_ids);
        y = y_vals(t, set_ids);
        y_err = y_val_errs(t,set_ids);
        
        errorbar(FreqFigAx{SetIndex}, x,y,y_err,'Color','k','Capsize',0)
        hold on
        scatter(FreqFigAx{SetIndex}, x,y,'MarkerFaceColor',cmap(t,:),'MarkerEdgeColor','k')
        
    end
    
    
    grid on
    
    SubplotIndex = SubplotIndexList(SetIndex);
    ColumnIndex = mod(SubplotIndex, SubplotDims(2));
    if ColumnIndex == 0
        ColumnIndex = SubplotDims(2);
    end
    RowIndex = fix(SubplotIndex/SubplotDims(2))+1;
    %disp([num2str(RowIndex), ', ', num2str(ColumnIndex)]);
    pos = get(FreqFigAx{SetIndex}, 'position');
    pos(1) = SubplotXPositions(ColumnIndex)+SubplotXBuffer;
    pos(3) = SubplotWidth-SubplotXBuffer;
    pos(2) = SubplotYPositions(end-(RowIndex-1))+SubplotYBuffer;
    pos(4) = SubplotHeight-SubplotYBuffer;
    set(FreqFigAx{SetIndex}, 'position', pos);
    
    xlabel('AP position')
    ylabel('burst frequency (1/min)')
    set(FreqFigAx{SetIndex} ,'Fontsize',14)
    xlim([(MinAPbin-2)*APResolution*100 (MaxAPbin)*APResolution*100])
    %xlim([(MinAPbin-2)*APResolution*100 (MaxAPbin)*APResolution*100])
    title([CompiledParameters.ReporterLabels{SetIndex}, ' ', num2str(CompiledParameters.SetTemperatures(SetIndex)), '°C w', num2str(CompiledParameters.nSteps(SetIndex))]);
    newylim = get(FreqFigAx{SetIndex}, 'ylim');
    newylim(1) = FreqYmin;
    if newylim(2) > 2
        newylim(2)  = 2;
    end
    set(FreqFigAx{SetIndex}, 'ylim', newylim)
    %     ylim([DurYmin, DurYmax])
    
    
end

LegendAx = subplot(SubplotDims(1), SubplotDims(2), LegendSubplots);
hold on
cmap = flipud(brewermap(max(IncludedTimeBins),'Spectral'));
colormap(cmap);
if max(IncludedTimeBins) > 1
set(LegendAx, 'Clim', [1, max(IncludedTimeBins)])
% % AxesH = axes('CLim', [-12, 12]);
% cbh = colorbar('peer', AxesH, 'h', ...
%                'XTickLabel',{'-12','-9','-6','-3','0','3','6','9','12'}, ...
%                'XTick', -12:3:12)
h = colorbar;
h.Ticks = 1:max(IncludedTimeBins);
h.TickLabels = CompiledParameters.TimeVector(IncludedTimeBins);
ylabel(h,'time cohort (minutes into nc14)')
end
hold off
set(LegendAx,'Fontsize',16)
axis off

pos = get(LegendAx, 'Position');
pos(1) = .6;
pos(3) = .2;
set(LegendAx, 'Position', pos);


for SetIndex = 1:NumSets
    SubplotIndex = SubplotIndexList(SetIndex);
    ColumnIndex = mod(SubplotIndex, SubplotDims(2));
    if ColumnIndex == 0
        ColumnIndex = SubplotDims(2);
    end
    RowIndex = fix(SubplotIndex/SubplotDims(2))+1;
    %disp([num2str(RowIndex), ', ', num2str(ColumnIndex)]);
    pos = get(FreqFigAx{SetIndex}, 'position');
    pos(1) = SubplotXPositions(ColumnIndex)+SubplotXBuffer;
    pos(3) = SubplotWidth-SubplotXBuffer;
    pos(2) = SubplotYPositions(end-(RowIndex-1))+SubplotYBuffer;
    pos(4) = SubplotHeight-SubplotYBuffer;
    set(FreqFigAx{SetIndex}, 'position', pos);
end

saveas(freq_fig,[outdir, 'Projectsubplots_burst_frequency.png'])

%%


InitFigAx = cell(1, NumSets);

init_fig = figure(3);

set(init_fig,'units', 'normalized', 'position',[0.01, 0.01, SubFigDims(1), SubFigDims(2)]);
set(gcf,'color','w');
if max(IncludedTimeBins) > 1
cmap = flipud(brewermap(max(IncludedTimeBins),'Spectral'));
else
    cmap = flipud(brewermap(2,'Spectral'));
end
colormap(cmap);
hold on
SubplotIndex = 0;
SubplotIndexList = zeros(1, NumSets);
for SetIndex = 1:NumSets
    SubplotIndex = SubplotIndex + 1;
    if SubplotIndex == 1
        InitFigAx{SetIndex} = subplot(SubplotDims(1), SubplotDims(2), SubplotIndex, gca);
    else
        if mod(SubplotIndex,  SubplotDims(2)) == 0
            SubplotIndex = SubplotIndex+1;
        end
        InitFigAx{SetIndex} = subplot(SubplotDims(1), SubplotDims(2), SubplotIndex);
    end
    SubplotIndexList(SetIndex) = SubplotIndex;
    PlottedTimes = zeros(1, length(CompiledParameters.TimeVector), 'logical');
    
    y_vals = squeeze(CompiledParameters.InitiationRates(SetIndex,:,:));
    y_val_errs = squeeze(CompiledParameters.InitiationRatesStdErr(SetIndex,:,:));
    for t = 1:length(CompiledParameters.TimeVector)
        if all(isnan(y_vals(t,:)))
            continue
        end
        set_ids = ~isnan(y_vals(t,:));
        
        x = APbins(set_ids);
        y = y_vals(t, set_ids);
        y_err = y_val_errs(t,set_ids);
        
        errorbar(InitFigAx{SetIndex}, x,y,y_err,'Color','k','Capsize',0)
        hold on
        scatter(InitFigAx{SetIndex}, x,y,'MarkerFaceColor',cmap(t,:),'MarkerEdgeColor','k')
        
    end
    
    
    grid on
    
    SubplotIndex = SubplotIndexList(SetIndex);
    ColumnIndex = mod(SubplotIndex, SubplotDims(2));
    if ColumnIndex == 0
        ColumnIndex = SubplotDims(2);
    end
    RowIndex = fix(SubplotIndex/SubplotDims(2))+1;
    %disp([num2str(RowIndex), ', ', num2str(ColumnIndex)]);
    pos = get(InitFigAx{SetIndex}, 'position');
    pos(1) = SubplotXPositions(ColumnIndex)+SubplotXBuffer;
    pos(3) = SubplotWidth-SubplotXBuffer;
    pos(2) = SubplotYPositions(end-(RowIndex-1))+SubplotYBuffer;
    pos(4) = SubplotHeight-SubplotYBuffer;
    set(InitFigAx{SetIndex}, 'position', pos);
    
    xlabel('AP position')
    ylabel('initiation rate (au/min)')
    set(InitFigAx{SetIndex} ,'Fontsize',14)
    xlim([(MinAPbin-2)*APResolution*100 (MaxAPbin)*APResolution*100])
    %xlim([(MinAPbin-2)*APResolution*100 (MaxAPbin)*APResolution*100])
     title([CompiledParameters.ReporterLabels{SetIndex}, ' ', num2str(CompiledParameters.SetTemperatures(SetIndex)), '°C w', num2str(CompiledParameters.nSteps(SetIndex))]);
    newylim = get(InitFigAx{SetIndex}, 'ylim');
    newylim(1) = InitYmin;
   
        newylim(2)  = InitYmax;

    set(InitFigAx{SetIndex}, 'ylim', newylim)
    %     ylim([DurYmin, DurYmax])
    
    
end

LegendAx = subplot(SubplotDims(1), SubplotDims(2), LegendSubplots);
hold on
if max(IncludedTimeBins) > 1
cmap = flipud(brewermap(max(IncludedTimeBins),'Spectral'));
else
   cmap = flipud(brewermap(2,'Spectral'));
end
colormap(cmap);
if max(IncludedTimeBins) > 1
set(LegendAx, 'Clim', [1, max(IncludedTimeBins)])
% % AxesH = axes('CLim', [-12, 12]);
% cbh = colorbar('peer', AxesH, 'h', ...
%                'XTickLabel',{'-12','-9','-6','-3','0','3','6','9','12'}, ...
%                'XTick', -12:3:12)
h = colorbar;
h.Ticks = 1:length(IncludedTimeBins);
h.TickLabels = CompiledParameters.TimeVector(IncludedTimeBins);
ylabel(h,'time cohort (minutes into nc14)')
end
hold off
set(LegendAx,'Fontsize',16)
axis off

pos = get(LegendAx, 'Position');
pos(1) = .6;
pos(3) = .2;
set(LegendAx, 'Position', pos);


for SetIndex = 1:NumSets
    SubplotIndex = SubplotIndexList(SetIndex);
    ColumnIndex = mod(SubplotIndex, SubplotDims(2));
    if ColumnIndex == 0
        ColumnIndex = SubplotDims(2);
    end
    RowIndex = fix(SubplotIndex/SubplotDims(2))+1;
    %disp([num2str(RowIndex), ', ', num2str(ColumnIndex)]);
    pos = get(InitFigAx{SetIndex}, 'position');
    pos(1) = SubplotXPositions(ColumnIndex)+SubplotXBuffer;
    pos(3) = SubplotWidth-SubplotXBuffer;
    pos(2) = SubplotYPositions(end-(RowIndex-1))+SubplotYBuffer;
    pos(4) = SubplotHeight-SubplotYBuffer;
    set(InitFigAx{SetIndex}, 'position', pos);
end

saveas(init_fig,[outdir, 'Projectsubplots_burst_initiation.png'])







