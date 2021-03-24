function GenerateAPSubplots(ResultsPaths, outdir)
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


SubplotDims = [3, 3];
SubplotDims(2) = SubplotDims(2)+1;

SubFigDims = [0.85, 0.9];% 0.9*3072/1920*SubplotDims(1)/SubplotDims(2)*1.2];
SubFigDims(2) = min([SubFigDims(2), 0.95]);
SubFigDims = round(SubFigDims, 2);

LegendSubplots = (1:SubplotDims(1))*SubplotDims(2);


LegendWidth = 0.2;
LegendXPosition = 0.78;
SubplotXBuffer = 0.05;
SubplotYBuffer = 0.1;

SubplotWidth = .23; %(LegendXPosition-(0.05-SubplotXBuffer))/(SubplotDims(2)-1);
SubplotHeight = .3;  %(1-(0.15-SubplotYBuffer))/SubplotDims(1);
SubplotXPositions = 0.025+(0:(SubplotDims(2)-2))*(SubplotWidth);
SubplotYPositions = 0.0+(0:(SubplotDims(1)-1))*(SubplotHeight);

SubplotYPositions = [0, .33, .66];

DurFigAx = cell(1, length(IncludedAPbins));

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
SubplotIndexList = zeros(1,  length(IncludedAPbins));
for i = 1:length(IncludedAPbins)
    APindex = IncludedAPbins(i);
    APbin = APbins(APindex);
    
    SubplotIndex = SubplotIndex + 1;
    if SubplotIndex == 1
        DurFigAx{i} = subplot(SubplotDims(1), SubplotDims(2), SubplotIndex, gca);
    else
        if mod(SubplotIndex,  SubplotDims(2)) == 0
            SubplotIndex = SubplotIndex+1;
        end
        DurFigAx{i} = subplot(SubplotDims(1), SubplotDims(2), SubplotIndex);
    end
    SubplotIndexList(i) = SubplotIndex;
    PlottedTimes = zeros(1, length(CompiledParameters.TimeVector), 'logical');
    
    y_vals = squeeze(CompiledParameters.Durations(:,:,APindex));
    y_val_errs = squeeze(CompiledParameters.DurationsStdErr(:,:,APindex));
    for t = 1:length(CompiledParameters.TimeVector)
        if all(isnan(y_vals(:,t)))
            continue
        end
        set_ids = ~isnan(y_vals(:,t)).';
        
        x = CompiledParameters.SetTemperatures(set_ids);
        y = y_vals(set_ids, t).';
        y_err = y_val_errs(set_ids, t).';
        
        errorbar(DurFigAx{i}, x,y,y_err,'Color','k','Capsize',0)
        hold on
        scatter(DurFigAx{i}, x,y,'MarkerFaceColor',cmap(t,:),'MarkerEdgeColor','k')
        
    end
    
    
    grid on
    
    SubplotIndex = SubplotIndexList(i);
    ColumnIndex = mod(SubplotIndex, SubplotDims(2));
    if ColumnIndex == 0
        ColumnIndex = SubplotDims(2);
    end
    RowIndex = fix(SubplotIndex/SubplotDims(2))+1;
    %disp([num2str(RowIndex), ', ', num2str(ColumnIndex)]);
    pos = get(DurFigAx{i}, 'position');
    pos(1) = SubplotXPositions(ColumnIndex)+SubplotXBuffer;
    pos(3) = SubplotWidth-SubplotXBuffer;
    pos(2) = SubplotYPositions(end-(RowIndex-1))+SubplotYBuffer;
    pos(4) = SubplotHeight-SubplotYBuffer;
    set(DurFigAx{i}, 'position', pos);
    
    xlabel('Temperature (°C)')
    ylabel('duration (min)')
    set(DurFigAx{i} ,'Fontsize',14)
    xlim([16, 29])
    %xlim([(MinAPbin-2)*APResolution*100 (MaxAPbin)*APResolution*100])
    title(['AP Position:  ', num2str(APbin), '%']);
    newylim = get(DurFigAx{i}, 'ylim');
    newylim(1) = DurYmin;
    if newylim(2) > 15
        newylim(2)  = 15;
    end
    set(DurFigAx{i}, 'ylim', newylim)
    %     ylim([DurYmin, DurYmax])
    
    
end

LegendAx = subplot(SubplotDims(1), SubplotDims(2), LegendSubplots);
hold on
cmap = flipud(brewermap(max(IncludedTimeBins),'Spectral'));
colormap(cmap);

% % AxesH = axes('CLim', [-12, 12]);
% cbh = colorbar('peer', AxesH, 'h', ...
%                'XTickLabel',{'-12','-9','-6','-3','0','3','6','9','12'}, ...
%                'XTick', -12:3:12)
if max(IncludedTimeBins) > 1
    set(LegendAx, 'Clim', [1, max(IncludedTimeBins)])
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


for i = 1:length(IncludedAPbins)
    SubplotIndex = SubplotIndexList(i);
    ColumnIndex = mod(SubplotIndex, SubplotDims(2));
    if ColumnIndex == 0
        ColumnIndex = SubplotDims(2);
    end
    RowIndex = fix(SubplotIndex/SubplotDims(2))+1;
    %disp([num2str(RowIndex), ', ', num2str(ColumnIndex)]);
    pos = get(DurFigAx{i}, 'position');
    pos(1) = SubplotXPositions(ColumnIndex)+SubplotXBuffer;
    pos(3) = SubplotWidth-SubplotXBuffer;
    pos(2) = SubplotYPositions(end-(RowIndex-1))+SubplotYBuffer;
    pos(4) = SubplotHeight-SubplotYBuffer;
    set(DurFigAx{i}, 'position', pos);
end

saveas(dur_fig,[outdir, 'APsubplots_burst_duration.png'])

%%
LogDurFigAx = cell(1, length(IncludedAPbins));
log_dur_fig = figure(2);

set(log_dur_fig,'units', 'normalized', 'position',[0.01, 0.01, SubFigDims(1), SubFigDims(2)]);
set(gcf,'color','w');
if max(IncludedTimeBins) > 1
cmap = flipud(brewermap(max(IncludedTimeBins),'Spectral'));
else
    cmap = flipud(brewermap(2,'Spectral'));
end
colormap(cmap);
hold on
SubplotIndex = 0;
SubplotIndexList = zeros(1,  length(IncludedAPbins));
for i = 1:length(IncludedAPbins)
    APindex = IncludedAPbins(i);
    APbin = APbins(APindex);
    
    SubplotIndex = SubplotIndex + 1;
    if SubplotIndex == 1
        LogDurFigAx{i} = subplot(SubplotDims(1), SubplotDims(2), SubplotIndex, gca);
    else
        if mod(SubplotIndex,  SubplotDims(2)) == 0
            SubplotIndex = SubplotIndex+1;
        end
        LogDurFigAx{i} = subplot(SubplotDims(1), SubplotDims(2), SubplotIndex);
    end
    SubplotIndexList(i) = SubplotIndex;
    PlottedTimes = zeros(1, length(CompiledParameters.TimeVector), 'logical');
    
    y_vals = squeeze(CompiledParameters.Durations(:,:,APindex));
    y_val_errs = squeeze(CompiledParameters.DurationsStdErr(:,:,APindex));
    for t = 1:length(CompiledParameters.TimeVector)
        if all(isnan(y_vals(:,t)))
            continue
        end
        set_ids = ~isnan(y_vals(:,t)).';
        
        x = 1./(R*(CompiledParameters.SetTemperatures(set_ids)+273.15));
        y = y_vals(set_ids, t).';
        y_err = y_val_errs(set_ids, t).';
        
        errorbar(LogDurFigAx{i}, x,y,y_err,'Color','k','Capsize',0)
        hold on
        scatter(LogDurFigAx{i}, x,y,'MarkerFaceColor',cmap(t,:),'MarkerEdgeColor','k')
        
    end
    
    
    grid on
    
    SubplotIndex = SubplotIndexList(i);
    ColumnIndex = mod(SubplotIndex, SubplotDims(2));
    if ColumnIndex == 0
        ColumnIndex = SubplotDims(2);
    end
    RowIndex = fix(SubplotIndex/SubplotDims(2))+1;
    %disp([num2str(RowIndex), ', ', num2str(ColumnIndex)]);
    pos = get(LogDurFigAx{i}, 'position');
    pos(1) = SubplotXPositions(ColumnIndex)+SubplotXBuffer;
    pos(3) = SubplotWidth-SubplotXBuffer;
    pos(2) = SubplotYPositions(end-(RowIndex-1))+SubplotYBuffer;
    pos(4) = SubplotHeight-SubplotYBuffer;
    set(LogDurFigAx{i}, 'position', pos);
    
    xlabel('1/(RT) (mol/kJ)')
    ylabel('duration (min)')
    set(LogDurFigAx{i} ,'Fontsize',14)
    xlim([0.398, 0.415])
    %xlim([(MinAPbin-2)*APResolution*100 (MaxAPbin)*APResolution*100])
    title(['AP Position:  ', num2str(APbin), '%']);
    newylim = get(LogDurFigAx{i}, 'ylim');
    newylim(1) = max(DurYmin, 0.5);
    if newylim(2) > 15
        newylim(2)  = 15;
    end
    set(LogDurFigAx{i}, 'ylim', newylim)
    set(LogDurFigAx{i}, 'YScale', 'log')
    
    %     ylim([DurYmin, DurYmax])
    
    
end

LegendAx = subplot(SubplotDims(1), SubplotDims(2), LegendSubplots);
hold on
cmap = flipud(brewermap(max(IncludedTimeBins),'Spectral'));
colormap(cmap);

% % AxesH = axes('CLim', [-12, 12]);
% cbh = colorbar('peer', AxesH, 'h', ...
%                'XTickLabel',{'-12','-9','-6','-3','0','3','6','9','12'}, ...
%                'XTick', -12:3:12)
if max(IncludedTimeBins) > 1
    set(LegendAx, 'Clim', [1, max(IncludedTimeBins)])
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


for i = 1:length(IncludedAPbins)
    SubplotIndex = SubplotIndexList(i);
    ColumnIndex = mod(SubplotIndex, SubplotDims(2));
    if ColumnIndex == 0
        ColumnIndex = SubplotDims(2);
    end
    RowIndex = fix(SubplotIndex/SubplotDims(2))+1;
    %disp([num2str(RowIndex), ', ', num2str(ColumnIndex)]);
    pos = get(LogDurFigAx{i}, 'position');
    pos(1) = SubplotXPositions(ColumnIndex)+SubplotXBuffer;
    pos(3) = SubplotWidth-SubplotXBuffer;
    pos(2) = SubplotYPositions(end-(RowIndex-1))+SubplotYBuffer;
    pos(4) = SubplotHeight-SubplotYBuffer;
    set(LogDurFigAx{i}, 'position', pos);
end

saveas(log_dur_fig,[outdir, 'APLogsubplots_burst_duration.png'])

%%

FreqFigAx = cell(1, length(IncludedAPbins));

freq_fig = figure(3);

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
SubplotIndexList = zeros(1,  length(IncludedAPbins));
for i = 1:length(IncludedAPbins)
    APindex = IncludedAPbins(i);
    APbin = APbins(APindex);
    
    SubplotIndex = SubplotIndex + 1;
    if SubplotIndex == 1
        FreqFigAx{i} = subplot(SubplotDims(1), SubplotDims(2), SubplotIndex, gca);
    else
        if mod(SubplotIndex,  SubplotDims(2)) == 0
            SubplotIndex = SubplotIndex+1;
        end
        FreqFigAx{i} = subplot(SubplotDims(1), SubplotDims(2), SubplotIndex);
    end
    SubplotIndexList(i) = SubplotIndex;
    PlottedTimes = zeros(1, length(CompiledParameters.TimeVector), 'logical');
    
    y_vals = squeeze(CompiledParameters.Frequencies(:,:,APindex));
    y_val_errs = squeeze(CompiledParameters.FrequenciesStdErr(:,:,APindex));
    for t = 1:length(CompiledParameters.TimeVector)
        if all(isnan(y_vals(:,t)))
            continue
        end
        set_ids = ~isnan(y_vals(:,t)).';
        
        x = CompiledParameters.SetTemperatures(set_ids);
        y = y_vals(set_ids, t).';
        y_err = y_val_errs(set_ids, t).';
        
        errorbar(FreqFigAx{i}, x,y,y_err,'Color','k','Capsize',0)
        hold on
        scatter(FreqFigAx{i}, x,y,'MarkerFaceColor',cmap(t,:),'MarkerEdgeColor','k')
        
    end
    
    
    grid on
    
    SubplotIndex = SubplotIndexList(i);
    ColumnIndex = mod(SubplotIndex, SubplotDims(2));
    if ColumnIndex == 0
        ColumnIndex = SubplotDims(2);
    end
    RowIndex = fix(SubplotIndex/SubplotDims(2))+1;
    %disp([num2str(RowIndex), ', ', num2str(ColumnIndex)]);
    pos = get(FreqFigAx{i}, 'position');
    pos(1) = SubplotXPositions(ColumnIndex)+SubplotXBuffer;
    pos(3) = SubplotWidth-SubplotXBuffer;
    pos(2) = SubplotYPositions(end-(RowIndex-1))+SubplotYBuffer;
    pos(4) = SubplotHeight-SubplotYBuffer;
    set(FreqFigAx{i}, 'position', pos);
    
    xlabel('Temperature (°C)')
    ylabel('frequency (1/min)')
    set(FreqFigAx{i} ,'Fontsize',14)
    xlim([16, 29])
    %xlim([(MinAPbin-2)*APResolution*100 (MaxAPbin)*APResolution*100])
    title(['AP Position:  ', num2str(APbin), '%']);
    newylim = get(FreqFigAx{i}, 'ylim');
    newylim(1) = FreqYmin;
    
    newylim(2)  = 2;
    
    set(FreqFigAx{i}, 'ylim', newylim)
    %     ylim([DurYmin, DurYmax])
    
    
end

LegendAx = subplot(SubplotDims(1), SubplotDims(2), LegendSubplots);
hold on
cmap = flipud(brewermap(max(IncludedTimeBins),'Spectral'));
colormap(cmap);
% % AxesH = axes('CLim', [-12, 12]);
% cbh = colorbar('peer', AxesH, 'h', ...
%                'XTickLabel',{'-12','-9','-6','-3','0','3','6','9','12'}, ...
%                'XTick', -12:3:12)
if max(IncludedTimeBins) > 1
    set(LegendAx, 'Clim', [1, max(IncludedTimeBins)])
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


for i = 1:length(IncludedAPbins)
    SubplotIndex = SubplotIndexList(i);
    ColumnIndex = mod(SubplotIndex, SubplotDims(2));
    if ColumnIndex == 0
        ColumnIndex = SubplotDims(2);
    end
    RowIndex = fix(SubplotIndex/SubplotDims(2))+1;
    %disp([num2str(RowIndex), ', ', num2str(ColumnIndex)]);
    pos = get(FreqFigAx{i}, 'position');
    pos(1) = SubplotXPositions(ColumnIndex)+SubplotXBuffer;
    pos(3) = SubplotWidth-SubplotXBuffer;
    pos(2) = SubplotYPositions(end-(RowIndex-1))+SubplotYBuffer;
    pos(4) = SubplotHeight-SubplotYBuffer;
    set(FreqFigAx{i}, 'position', pos);
end

saveas(freq_fig,[outdir, 'APsubplots_burst_frequency.png'])


%%
LogFreqFigAx = cell(1, length(IncludedAPbins));
log_freq_fig = figure(4);

set(log_freq_fig,'units', 'normalized', 'position',[0.01, 0.01, SubFigDims(1), SubFigDims(2)]);
set(gcf,'color','w');
if max(IncludedTimeBins) > 1
cmap = flipud(brewermap(max(IncludedTimeBins),'Spectral'));
else
    cmap = flipud(brewermap(2,'Spectral'));
end
colormap(cmap);
hold on
SubplotIndex = 0;
SubplotIndexList = zeros(1,  length(IncludedAPbins));
for i = 1:length(IncludedAPbins)
    APindex = IncludedAPbins(i);
    APbin = APbins(APindex);
    
    SubplotIndex = SubplotIndex + 1;
    if SubplotIndex == 1
        LogFreqFigAx{i} = subplot(SubplotDims(1), SubplotDims(2), SubplotIndex, gca);
    else
        if mod(SubplotIndex,  SubplotDims(2)) == 0
            SubplotIndex = SubplotIndex+1;
        end
        LogFreqFigAx{i} = subplot(SubplotDims(1), SubplotDims(2), SubplotIndex);
    end
    SubplotIndexList(i) = SubplotIndex;
    PlottedTimes = zeros(1, length(CompiledParameters.TimeVector), 'logical');
    
    y_vals = squeeze(CompiledParameters.Frequencies(:,:,APindex));
    y_val_errs = squeeze(CompiledParameters.FrequenciesStdErr(:,:,APindex));
    for t = 1:length(CompiledParameters.TimeVector)
        if all(isnan(y_vals(:,t)))
            continue
        end
        set_ids = ~isnan(y_vals(:,t)).';
        
        x = 1./(R*(CompiledParameters.SetTemperatures(set_ids)+273.15));
        y = y_vals(set_ids, t).';
        y_err = y_val_errs(set_ids, t).';
        
        errorbar(LogFreqFigAx{i}, x,y,y_err,'Color','k','Capsize',0)
        hold on
        scatter(LogFreqFigAx{i}, x,y,'MarkerFaceColor',cmap(t,:),'MarkerEdgeColor','k')
        
    end
    
    
    grid on
    
    SubplotIndex = SubplotIndexList(i);
    ColumnIndex = mod(SubplotIndex, SubplotDims(2));
    if ColumnIndex == 0
        ColumnIndex = SubplotDims(2);
    end
    RowIndex = fix(SubplotIndex/SubplotDims(2))+1;
    %disp([num2str(RowIndex), ', ', num2str(ColumnIndex)]);
    pos = get(LogFreqFigAx{i}, 'position');
    pos(1) = SubplotXPositions(ColumnIndex)+SubplotXBuffer;
    pos(3) = SubplotWidth-SubplotXBuffer;
    pos(2) = SubplotYPositions(end-(RowIndex-1))+SubplotYBuffer;
    pos(4) = SubplotHeight-SubplotYBuffer;
    set(LogFreqFigAx{i}, 'position', pos);
    
    xlabel('1/(RT) (mol/kJ)')
    ylabel('frequency (1/min)')
    set(LogFreqFigAx{i} ,'Fontsize',14)
    xlim([0.398, 0.415])
    %xlim([(MinAPbin-2)*APResolution*100 (MaxAPbin)*APResolution*100])
    title(['AP Position:  ', num2str(APbin), '%']);
    newylim = get(LogFreqFigAx{i}, 'ylim');
    newylim(1) = max(FreqYmin, 0.01);
    
    newylim(2)  = 2;
    
    set(LogFreqFigAx{i}, 'ylim', newylim)
    set(LogFreqFigAx{i}, 'YScale', 'log')
    
    %     ylim([DurYmin, DurYmax])
    
    
end

LegendAx = subplot(SubplotDims(1), SubplotDims(2), LegendSubplots);
hold on
cmap = flipud(brewermap(max(IncludedTimeBins),'Spectral'));
colormap(cmap);
% % AxesH = axes('CLim', [-12, 12]);
% cbh = colorbar('peer', AxesH, 'h', ...
%                'XTickLabel',{'-12','-9','-6','-3','0','3','6','9','12'}, ...
%                'XTick', -12:3:12)
if max(IncludedTimeBins) > 1
    set(LegendAx, 'Clim', [1, max(IncludedTimeBins)])
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


for i = 1:length(IncludedAPbins)
    SubplotIndex = SubplotIndexList(i);
    ColumnIndex = mod(SubplotIndex, SubplotDims(2));
    if ColumnIndex == 0
        ColumnIndex = SubplotDims(2);
    end
    RowIndex = fix(SubplotIndex/SubplotDims(2))+1;
    %disp([num2str(RowIndex), ', ', num2str(ColumnIndex)]);
    pos = get(LogFreqFigAx{i}, 'position');
    pos(1) = SubplotXPositions(ColumnIndex)+SubplotXBuffer;
    pos(3) = SubplotWidth-SubplotXBuffer;
    pos(2) = SubplotYPositions(end-(RowIndex-1))+SubplotYBuffer;
    pos(4) = SubplotHeight-SubplotYBuffer;
    set(LogFreqFigAx{i}, 'position', pos);
end

saveas(log_freq_fig,[outdir, 'APLogsubplots_burst_frequency.png'])


%%
InitFigAx = cell(1, length(IncludedAPbins));

init_fig = figure(5);

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
SubplotIndexList = zeros(1,  length(IncludedAPbins));
for i = 1:length(IncludedAPbins)
    APindex = IncludedAPbins(i);
    APbin = APbins(APindex);
    
    SubplotIndex = SubplotIndex + 1;
    if SubplotIndex == 1
        InitFigAx{i} = subplot(SubplotDims(1), SubplotDims(2), SubplotIndex, gca);
    else
        if mod(SubplotIndex,  SubplotDims(2)) == 0
            SubplotIndex = SubplotIndex+1;
        end
        InitFigAx{i} = subplot(SubplotDims(1), SubplotDims(2), SubplotIndex);
    end
    SubplotIndexList(i) = SubplotIndex;
    PlottedTimes = zeros(1, length(CompiledParameters.TimeVector), 'logical');
    
    y_vals = squeeze(CompiledParameters.InitiationRates(:,:,APindex));
    y_val_errs = squeeze(CompiledParameters.InitiationRatesStdErr(:,:,APindex));
    for t = 1:length(CompiledParameters.TimeVector)
        if all(isnan(y_vals(:,t)))
            continue
        end
        set_ids = ~isnan(y_vals(:,t)).';
        
        x = CompiledParameters.SetTemperatures(set_ids);
        y = y_vals(set_ids, t).';
        y_err = y_val_errs(set_ids, t).';
        
        errorbar(InitFigAx{i}, x,y,y_err,'Color','k','Capsize',0)
        hold on
        scatter(InitFigAx{i}, x,y,'MarkerFaceColor',cmap(t,:),'MarkerEdgeColor','k')
        
    end
    
    
    grid on
    
    SubplotIndex = SubplotIndexList(i);
    ColumnIndex = mod(SubplotIndex, SubplotDims(2));
    if ColumnIndex == 0
        ColumnIndex = SubplotDims(2);
    end
    RowIndex = fix(SubplotIndex/SubplotDims(2))+1;
    %disp([num2str(RowIndex), ', ', num2str(ColumnIndex)]);
    pos = get(InitFigAx{i}, 'position');
    pos(1) = SubplotXPositions(ColumnIndex)+SubplotXBuffer;
    pos(3) = SubplotWidth-SubplotXBuffer;
    pos(2) = SubplotYPositions(end-(RowIndex-1))+SubplotYBuffer;
    pos(4) = SubplotHeight-SubplotYBuffer;
    set(InitFigAx{i}, 'position', pos);
    
    xlabel('Temperature (°C)')
    ylabel('init. rate (au/min)')
    set(InitFigAx{i} ,'Fontsize',14)
    xlim([16, 29])
    %xlim([(MinAPbin-2)*APResolution*100 (MaxAPbin)*APResolution*100])
    title(['AP Position:  ', num2str(APbin), '%']);
    newylim = get(InitFigAx{i}, 'ylim');
    newylim(1) = InitYmin;
    
    newylim(2)  = InitYmax;
    
    set(InitFigAx{i}, 'ylim', newylim)
    %     ylim([DurYmin, DurYmax])
    
    
end

LegendAx = subplot(SubplotDims(1), SubplotDims(2), LegendSubplots);
hold on
cmap = flipud(brewermap(max(IncludedTimeBins),'Spectral'));
colormap(cmap);
% % AxesH = axes('CLim', [-12, 12]);
% cbh = colorbar('peer', AxesH, 'h', ...
%                'XTickLabel',{'-12','-9','-6','-3','0','3','6','9','12'}, ...
%                'XTick', -12:3:12)
if max(IncludedTimeBins) > 1
    set(LegendAx, 'Clim', [1, max(IncludedTimeBins)])
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


for i = 1:length(IncludedAPbins)
    SubplotIndex = SubplotIndexList(i);
    ColumnIndex = mod(SubplotIndex, SubplotDims(2));
    if ColumnIndex == 0
        ColumnIndex = SubplotDims(2);
    end
    RowIndex = fix(SubplotIndex/SubplotDims(2))+1;
    %disp([num2str(RowIndex), ', ', num2str(ColumnIndex)]);
    pos = get(InitFigAx{i}, 'position');
    pos(1) = SubplotXPositions(ColumnIndex)+SubplotXBuffer;
    pos(3) = SubplotWidth-SubplotXBuffer;
    pos(2) = SubplotYPositions(end-(RowIndex-1))+SubplotYBuffer;
    pos(4) = SubplotHeight-SubplotYBuffer;
    set(InitFigAx{i}, 'position', pos);
end

saveas(init_fig,[outdir, 'APsubplots_burst_initiation.png'])

%%
LogInitFigAx = cell(1, length(IncludedAPbins));
log_init_fig = figure(6);

set(log_init_fig,'units', 'normalized', 'position',[0.01, 0.01, SubFigDims(1), SubFigDims(2)]);
set(gcf,'color','w');
if max(IncludedTimeBins) > 1
cmap = flipud(brewermap(max(IncludedTimeBins),'Spectral'));
else
    cmap = flipud(brewermap(2,'Spectral'));
end
colormap(cmap);
hold on
SubplotIndex = 0;
SubplotIndexList = zeros(1,  length(IncludedAPbins));
for i = 1:length(IncludedAPbins)
    APindex = IncludedAPbins(i);
    APbin = APbins(APindex);
    
    SubplotIndex = SubplotIndex + 1;
    if SubplotIndex == 1
        LogInitFigAx{i} = subplot(SubplotDims(1), SubplotDims(2), SubplotIndex, gca);
    else
        if mod(SubplotIndex,  SubplotDims(2)) == 0
            SubplotIndex = SubplotIndex+1;
        end
        LogInitFigAx{i} = subplot(SubplotDims(1), SubplotDims(2), SubplotIndex);
    end
    SubplotIndexList(i) = SubplotIndex;
    PlottedTimes = zeros(1, length(CompiledParameters.TimeVector), 'logical');
    
    y_vals = squeeze(CompiledParameters.InitiationRates(:,:,APindex));
    y_val_errs = squeeze(CompiledParameters.InitiationRatesStdErr(:,:,APindex));
    for t = 1:length(CompiledParameters.TimeVector)
        if all(isnan(y_vals(:,t)))
            continue
        end
        set_ids = ~isnan(y_vals(:,t)).';
        
        x = 1./(R*(CompiledParameters.SetTemperatures(set_ids)+273.15));
        y = y_vals(set_ids, t).';
        y_err = y_val_errs(set_ids, t).';
        
        errorbar(LogInitFigAx{i}, x,y,y_err,'Color','k','Capsize',0)%, 'Linestyle', 'None')
        hold on
        scatter(LogInitFigAx{i}, x,y,'MarkerFaceColor',cmap(t,:),'MarkerEdgeColor','k')
        
    end
    
    
    grid on
    
    SubplotIndex = SubplotIndexList(i);
    ColumnIndex = mod(SubplotIndex, SubplotDims(2));
    if ColumnIndex == 0
        ColumnIndex = SubplotDims(2);
    end
    RowIndex = fix(SubplotIndex/SubplotDims(2))+1;
    %disp([num2str(RowIndex), ', ', num2str(ColumnIndex)]);
    pos = get(LogInitFigAx{i}, 'position');
    pos(1) = SubplotXPositions(ColumnIndex)+SubplotXBuffer;
    pos(3) = SubplotWidth-SubplotXBuffer;
    pos(2) = SubplotYPositions(end-(RowIndex-1))+SubplotYBuffer;
    pos(4) = SubplotHeight-SubplotYBuffer;
    set(LogInitFigAx{i}, 'position', pos);
    
    xlabel('1/(RT) (mol/kJ)')
    ylabel('init. rate (au/min)')
    set(LogInitFigAx{i} ,'Fontsize',14)
    xlim([0.398, 0.415])
    %xlim([(MinAPbin-2)*APResolution*100 (MaxAPbin)*APResolution*100])
    title(['AP Position:  ', num2str(APbin), '%']);
    newylim = get(LogInitFigAx{i}, 'ylim');
    newylim(1) = max(InitYmin, 10);
    
    newylim(2)  = InitYmax;
    
    set(LogInitFigAx{i}, 'ylim', newylim)
    set(LogInitFigAx{i}, 'YScale', 'log')
    
    %     ylim([DurYmin, DurYmax])
    
    
end

LegendAx = subplot(SubplotDims(1), SubplotDims(2), LegendSubplots);
hold on
cmap = flipud(brewermap(max(IncludedTimeBins),'Spectral'));
colormap(cmap);

% % AxesH = axes('CLim', [-12, 12]);
% cbh = colorbar('peer', AxesH, 'h', ...
%                'XTickLabel',{'-12','-9','-6','-3','0','3','6','9','12'}, ...
%                'XTick', -12:3:12)
if max(IncludedTimeBins) > 1
    set(LegendAx, 'Clim', [1, max(IncludedTimeBins)])
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


for i = 1:length(IncludedAPbins)
    SubplotIndex = SubplotIndexList(i);
    ColumnIndex = mod(SubplotIndex, SubplotDims(2));
    if ColumnIndex == 0
        ColumnIndex = SubplotDims(2);
    end
    RowIndex = fix(SubplotIndex/SubplotDims(2))+1;
    %disp([num2str(RowIndex), ', ', num2str(ColumnIndex)]);
    pos = get(LogInitFigAx{i}, 'position');
    pos(1) = SubplotXPositions(ColumnIndex)+SubplotXBuffer;
    pos(3) = SubplotWidth-SubplotXBuffer;
    pos(2) = SubplotYPositions(end-(RowIndex-1))+SubplotYBuffer;
    pos(4) = SubplotHeight-SubplotYBuffer;
    set(LogInitFigAx{i}, 'position', pos);
end

saveas(log_init_fig,[outdir, 'APLogsubplots_burst_initiation.png'])



