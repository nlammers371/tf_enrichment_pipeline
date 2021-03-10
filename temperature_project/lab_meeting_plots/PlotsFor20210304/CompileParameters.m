function CompiledParameters = CompileParameters(ResultsPaths)
%%


APResolution = 0.025;
APbins = 0:APResolution:1;
NumAPbins = length(APbins);
NumTimeBins = 18;
TimeBinWidth = 15;
TimeBinSep = 5;

NumSets = length(ResultsPaths);
prefix_regex = 'hbBAC-MS2-(?<temperature>[0-9_]+)C\\cpHMM_results\\compiledResults_w(?<nSteps>[0-9]+)_K(?<nStates>[0-9]+)_p0_ap(?<nAPbins>[0-9]+)_t(?<nTimebins>[0-9]+)_f2D.mat';
prefix_regex2 = 'hbBAC-MS2-(?<temperature>[0-9_]+)C\\NC(?<cycle>[0-9]+)\\cpHMM_results\\compiledResults_w(?<nSteps>[0-9]+)_K(?<nStates>[0-9]+)_p0_ap(?<nAPbins>[0-9]+)_t(?<nTimebins>[0-9]+)_f2D.mat';
prefix_regex3 = 'HbMS2JB3\\cpHMM_results\\compiledResults_w(?<nSteps>[0-9]+)_K(?<nStates>[0-9]+)_p0_ap(?<nAPbins>[0-9]+)_t(?<nTimebins>[0-9]+)_f2D.mat';
prefix_regex4 = 'HbMS2JB3\\NC(?<cycle>[0-9]+)\\cpHMM_results\\compiledResults_w(?<nSteps>[0-9]+)_K(?<nStates>[0-9]+)_p0_ap(?<nAPbins>[0-9]+)_t(?<nTimebins>[0-9]+)_f2D.mat';
CompiledParameters = {};
CompiledParameters.ResultsPaths = ResultsPaths;
CompiledParameters.FigurePaths = ResultsPaths;
CompiledParameters.ReporterLabels = cell(1, NumSets);
CompiledParameters.TimeLimits = NaN(NumTimeBins, 2);
CompiledParameters.TimeVector = NaN(1, NumTimeBins);
CompiledParameters.APVector = APbins;
for i = 1:NumTimeBins
    CompiledParameters.TimeLimits(i, 1) = (i-1)*5;
    CompiledParameters.TimeLimits(i, 2) = (i-1)*5+15;
    CompiledParameters.TimeVector(i) = mean(CompiledParameters.TimeLimits(i, :));
end

CompiledParameters.APRange =NaN(1, 2);
CompiledParameters.SetTemperatures = NaN(1,NumSets);
CompiledParameters.NC = NaN(1,NumSets);
CompiledParameters.nSteps = NaN(1,NumSets);
CompiledParameters.nStates = NaN(1,NumSets);
CompiledParameters.nAPbins = NaN(1,NumSets);
CompiledParameters.nTimebins = NaN(1,NumSets);
CompiledParameters.TimeAxes = cell(1, NumSets);
CompiledParameters.Durations = NaN(NumSets, NumTimeBins, NumAPbins);
CompiledParameters.Frequencies = NaN(NumSets, NumTimeBins, NumAPbins);
CompiledParameters.InitiationRates = NaN(NumSets, NumTimeBins, NumAPbins);
CompiledParameters.DurationsStdErr = NaN(NumSets, NumTimeBins, NumAPbins);
CompiledParameters.FrequenciesStdErr = NaN(NumSets, NumTimeBins, NumAPbins);
CompiledParameters.InitiationRatesStdErr = NaN(NumSets, NumTimeBins, NumAPbins);
EnrichmentPath = 'S:\Gabriella\Dropbox\ProcessedEnrichmentData\';
for ResultIndex=1:NumSets
    results_path_split = regexp(ResultsPaths{ResultIndex}, prefix_regex, 'names');
    results_path_split2 = regexp(ResultsPaths{ResultIndex}, prefix_regex2, 'names');
    results_path_split3 = regexp(ResultsPaths{ResultIndex}, prefix_regex3, 'names');
    results_path_split4 = regexp(ResultsPaths{ResultIndex}, prefix_regex4, 'names');
    if ~isempty(results_path_split)
        CompiledParameters.ReporterLabels{ResultIndex} = 'hbBAC-MS2';
        CompiledParameters.NC(ResultIndex) = 14;
        figpath = [EnrichmentPath, 'hbBAC-MS2-', results_path_split.temperature,'C',filesep, 'cpHMM_results',...
            filesep, 'figures', filesep, 'w', results_path_split.nSteps, '_K', results_path_split.nStates, '_p0_ap',...
            results_path_split.nAPbins, '_t', results_path_split.nTimebins, filesep];
        if ~exist(figpath, 'dir')
            mkdir(figpath);
        end
        CompiledParameters.FigurePaths{ResultIndex} = figpath;
        CompiledParameters.SetTemperatures(ResultIndex) = str2num(strrep(results_path_split.temperature, '_', '.'));
        CompiledParameters.nSteps(ResultIndex) = str2num(results_path_split.nSteps);
        CompiledParameters.nStates(ResultIndex) = str2num(results_path_split.nStates);
        CompiledParameters.nAPbins(ResultIndex) = str2num(results_path_split.nAPbins);
        CompiledParameters.nTimebins(ResultIndex) = str2num(results_path_split.nTimebins);
    elseif ~isempty(results_path_split2)
        CompiledParameters.ReporterLabels{ResultIndex} = 'hbBAC-MS2';
        CompiledParameters.NC(ResultIndex) = str2num(results_path_split2.cycle);
        figpath = [EnrichmentPath, 'hbBAC-MS2-', results_path_split2.temperature,'C',filesep,'NC', results_path_split2.cycle, filesep, 'cpHMM_results',...
            filesep, 'figures', filesep, 'w', results_path_split2.nSteps, '_K', results_path_split2.nStates, '_p0_ap',...
            results_path_split2.nAPbins, '_t', results_path_split2.nTimebins, filesep];
        if ~exist(figpath, 'dir')
            mkdir(figpath);
        end
        CompiledParameters.FigurePaths{ResultIndex} = figpath;
        CompiledParameters.SetTemperatures(ResultIndex) = str2num(strrep(results_path_split2.temperature, '_', '.'));
        CompiledParameters.nSteps(ResultIndex) = str2num(results_path_split2.nSteps);
        CompiledParameters.nStates(ResultIndex) = str2num(results_path_split2.nStates);
        CompiledParameters.nAPbins(ResultIndex) = str2num(results_path_split2.nAPbins);
        CompiledParameters.nTimebins(ResultIndex) = str2num(results_path_split2.nTimebins);
    elseif ~isempty(results_path_split3)
        CompiledParameters.ReporterLabels{ResultIndex} = 'Hb-MS2-JB3';
        
        CompiledParameters.SetTemperatures(ResultIndex) = 25;
        figpath = [EnrichmentPath, 'HbMS2JB3', filesep, 'cpHMM_results',...
            filesep, 'figures', filesep, 'w', results_path_split3.nSteps, '_K', results_path_split3.nStates, '_p0_ap',...
            results_path_split3.nAPbins, '_t', results_path_split3.nTimebins, filesep];
        if ~exist(figpath, 'dir')
            mkdir(figpath);
        end
        CompiledParameters.FigurePaths{ResultIndex} = figpath;
        CompiledParameters.nSteps(ResultIndex) = str2num(results_path_split3.nSteps);
        CompiledParameters.nStates(ResultIndex) = str2num(results_path_split3.nStates);
        CompiledParameters.nAPbins(ResultIndex) = str2num(results_path_split3.nAPbins);
        CompiledParameters.nTimebins(ResultIndex) = str2num(results_path_split3.nTimebins);
    else
        CompiledParameters.ReporterLabels{ResultIndex} = 'Hb-MS2-JB3';
        CompiledParameters.NC(ResultIndex) = str2num(results_path_split4.cycle);
        CompiledParameters.SetTemperatures(ResultIndex) = 25;
        figpath = [EnrichmentPath, 'HbMS2JB3', filesep,'NC', results_path_split4.cycle, filesep, 'cpHMM_results',...
            filesep, 'figures', filesep, 'w', results_path_split4.nSteps, '_K', results_path_split4.nStates, '_p0_ap',...
            results_path_split4.nAPbins, '_t', results_path_split4.nTimebins, filesep];
        if ~exist(figpath, 'dir')
            mkdir(figpath);
        end
        CompiledParameters.FigurePaths{ResultIndex} = figpath;
        CompiledParameters.nSteps(ResultIndex) = str2num(results_path_split4.nSteps);
        CompiledParameters.nStates(ResultIndex) = str2num(results_path_split4.nStates);
        CompiledParameters.nAPbins(ResultIndex) = str2num(results_path_split4.nAPbins);
        CompiledParameters.nTimebins(ResultIndex) = str2num(results_path_split4.nTimebins);
    end
    load([EnrichmentPath, ResultsPaths{ResultIndex}]);
    time_group_vec = compiledResults.timeGroupVec;
    time_group_index = unique(time_group_vec);
    ap_group_vec = compiledResults.apGroupVec;
    ap_group_index = unique(ap_group_vec);
    ap_axis = compiledResults.apBins(1:end-1) + diff(compiledResults.apBins);
    time_axis = [];
    for t = 1:length(time_group_index)
        time_axis(t) = mean(compiledResults.timeBins{t})/60;
    end
    CompiledParameters.TimeAxes{ResultIndex} = time_axis;
    
    for t = 1:length(time_group_index)
        time_filter = time_group_vec==t;
        ap_ids = ap_group_vec(time_filter);
        ap_bins = ap_axis(ap_ids)/100;
        [ap_bins, sort_order] = sort(ap_bins);
        AP_index = ismember(round(APbins, 3), round(ap_bins, 3));
        
        dur_vec_mean = compiledResults.dur_vec_mean(time_filter);
        dur_vec_mean = dur_vec_mean(sort_order);
        dur_vec_ste = compiledResults.dur_vec_ste(time_filter);
        dur_vec_ste = dur_vec_ste(sort_order);
        
        CompiledParameters.Durations(ResultIndex, t, AP_index) = dur_vec_mean;
        CompiledParameters.DurationsStdErr(ResultIndex, t, AP_index) = dur_vec_ste;
        
        freq_vec_mean = compiledResults.freq_vec_mean(time_filter);
        freq_vec_mean = freq_vec_mean(sort_order);
        freq_vec_ste = compiledResults.freq_vec_ste(time_filter);
        freq_vec_ste = freq_vec_ste(sort_order);
        
        CompiledParameters.Frequencies(ResultIndex, t, AP_index) = freq_vec_mean;
        CompiledParameters.FrequenciesStdErr(ResultIndex, t, AP_index) = freq_vec_ste;
        
        
        init_vec_mean = compiledResults.init_vec_mean(time_filter);
        init_vec_mean = init_vec_mean(sort_order);
        init_vec_ste = compiledResults.init_vec_ste(time_filter);
        init_vec_ste = init_vec_ste(sort_order);
        
        CompiledParameters.InitiationRates(ResultIndex, t, AP_index) = init_vec_mean;
        CompiledParameters.InitiationRatesStdErr(ResultIndex, t, AP_index) = init_vec_ste;
    end
    
end

%%
