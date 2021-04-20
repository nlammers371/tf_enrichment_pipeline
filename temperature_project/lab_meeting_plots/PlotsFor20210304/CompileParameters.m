function CompiledParameters = CompileParameters(ResultsPaths)
%%


APResolution = 0.025;
APbins = 0:APResolution:1;
NumAPbins = length(APbins);
NumTimeBins = 18;
TimeBinWidth = 15;
TimeBinSep = 5;

NumSets = length(ResultsPaths);


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
CompiledParameters.FluoDims = NaN(1, NumSets);
CompiledParameters.dt = NaN(1,NumSets);
CompiledParameters.TimeAxes = cell(1, NumSets);
CompiledParameters.ElongationTimes = NaN(1,NumSets);
CompiledParameters.Durations = NaN(NumSets, NumTimeBins, NumAPbins);
CompiledParameters.Frequencies = NaN(NumSets, NumTimeBins, NumAPbins);
CompiledParameters.InitiationRates = NaN(NumSets, NumTimeBins, NumAPbins);
CompiledParameters.DurationsStdErr = NaN(NumSets, NumTimeBins, NumAPbins);
CompiledParameters.FrequenciesStdErr = NaN(NumSets, NumTimeBins, NumAPbins);
CompiledParameters.InitiationRatesStdErr = NaN(NumSets, NumTimeBins, NumAPbins);
EnrichmentPath = 'S:\Gabriella\Dropbox\ProcessedEnrichmentData\';

%%

prefix_regex = '(?<projectDirInfo>[A-Za-z0-9\-_\\]+)\\cpHMM_results\\compiledResults_w(?<nSteps>[0-9]+)_K(?<nStates>[0-9]+)_p0_ap(?<nAPbins>[0-9]+)_t(?<nTimebins>[0-9]+)_f(?<FluoDim>[0-9]+)D(?<TimeRes>[a-z0-9._]+)';
cycle_regex = '(?<projectInfo>[A-Za-z\-0-9_]+)\\NC(?<cycle>[0-9]+)';
timeres_regex = '_dt(?<dt>[0-9]+).mat';
projectInfo_regex = 'hbBAC-MS2-(?<temperature>[0-9_]+)C';


%%
for ResultIndex=1:NumSets
    results_path_split = regexp(ResultsPaths{ResultIndex}, prefix_regex, 'names');
    if ~isempty(results_path_split)
        [fp1, fp2, fp3] = fileparts(ResultsPaths{ResultIndex});
        figpath = [EnrichmentPath,  filesep, fp1, filesep, 'figures', filesep, fp2, filesep];
        if ~exist(figpath, 'dir')
            mkdir(figpath);
        end
        CompiledParameters.FigurePaths{ResultIndex} = figpath;
        
        cycle_split = regexp(results_path_split.projectDirInfo, cycle_regex, 'names');
        if ~isempty(cycle_split)
            CompiledParameters.NC(ResultIndex) = str2num(cycle_split.cycle);
            proj_split_v1 = regexp(cycle_split.projectInfo, projectInfo_regex, 'names');
            PotentialProjName = cycle_split.projectInfo;
        else
            CompiledParameters.NC(ResultIndex) = 14;
            proj_split_v1 = regexp(results_path_split.projectDirInfo, projectInfo_regex, 'names');
            PotentialProjName = results_path_split.projectDirInfo;
        end
        
        if ~isempty(proj_split_v1)
            CompiledParameters.ReporterLabels{ResultIndex} = 'hbBAC-MS2';
            CompiledParameters.SetTemperatures(ResultIndex) = str2num(strrep(proj_split_v1.temperature, '_', '.'));
        else
            CompiledParameters.ReporterLabels{ResultIndex} = PotentialProjName;
            CompiledParameters.SetTemperatures(ResultIndex) = 25;
        end
            
        
        
        CompiledParameters.nSteps(ResultIndex) = str2num(results_path_split.nSteps);
        CompiledParameters.nStates(ResultIndex) = str2num(results_path_split.nStates);
        CompiledParameters.nAPbins(ResultIndex) = str2num(results_path_split.nAPbins);
        CompiledParameters.nTimebins(ResultIndex) = str2num(results_path_split.nTimebins);
        CompiledParameters.nTimebins(ResultIndex) = str2num(results_path_split.nTimebins);
        CompiledParameters.FluoDims(ResultIndex) = str2num(results_path_split.FluoDim);
        timeres_split = regexp(results_path_split.TimeRes, timeres_regex, 'names');
        if ~isempty(timeres_split)
            CompiledParameters.dt(ResultIndex) = str2num(timeres_split.dt);
        else
            datadir = strsplit(ResultsPaths{ResultIndex}, '\\cpHMM_results');
            load([EnrichmentPath datadir{1} filesep 'spot_struct.mat']);
            for i = 1:length(spot_struct)
                if length(spot_struct(i).timeInterp) > 1
                   if all(~isnan(spot_struct(i).timeInterp)) 
                       CompiledParameters.dt(ResultIndex) = spot_struct(i).timeInterp(2)-spot_struct(i).timeInterp(1);
                       break
                   end
                end
                       
            end
        end
        
        CompiledParameters.ElongationTimes(ResultIndex) = CompiledParameters.nSteps(ResultIndex)*CompiledParameters.dt(ResultIndex)/60;  
        
        
    else
        error('Inference Info path has unexpected format.')
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