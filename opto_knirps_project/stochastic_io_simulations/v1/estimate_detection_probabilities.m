% Script to scan through parameter space for selected biological variables
function io_ref_wt = estimateDetectionThreshold(io_ref_wt,spot_struct)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%% Load experimental data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% resultsRoot = 'S:\Nick\Dropbox\ProcessedEnrichmentData\parameterSweeps\';
% if ~isfolder(resultsRoot)
%   resultsRoot = 'C:\Users\nlamm\Dropbox (Personal)\ProcessedEnrichmentData\parameterSweeps\';
% end
% Load data
% projectNameWT = 'optokni_eve4+6_WT_FUN'; 

% liveProject = LiveEnrichmentProject(projectNameWT);
% resultsRoot = [liveProject.dataPath filesep];

% load data
% load([resultsRoot 'spot_struct.mat'],'spot_struct')

% load data
% load([resultsRoot 'io_ref_wt.mat'])
% load([resultsRoot 'io_ref_ra.mat'])

% resultsRoot = [resultsRoot 'temp' filesep];
% mkdir(resultsRoot);

% set basic parameters
sweepInfoRaw = struct;
sweepInfoRaw.granularity_wt = 1;
sweepInfoRaw.n_traces_per_ap = 50;
sweepInfoRaw.rate_max = 1; % nothing faster than a second
sweepInfoRaw.keep_prediction_flag = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%% Load cpHMM results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% designate simulation type
simType = 'match_exp';

sweepInfo = sweepInfoRaw;

% save simulation type
sweepInfo.simType = simType;

% load markov system parameter info
sweepInfo = getMarkovSystemInfo(sweepInfo);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Call Sweep function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all force 

% record "true" profile
myCluster = parcluster('local');
max_workers = myCluster.NumWorkers;
sweepInfo.NumWorkers = min([24 max_workers]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate ground truth reference curves       
sweepInfo = addGroundTruthFields(sweepInfo, [], io_ref_wt);

% initialize vectors to store results
sweepResults = struct;
[sweepInfo, sweepResults] = initializeFitFields(sweepInfo,sweepResults);
sweepResults = initializeSweepValues(sweepInfo, sweepResults);              

% call parallel sweep script
tic
sweepResults= sweep_par_loop_v3(sweepInfo,sweepResults,resultsRoot);    
toc

% recombine 
fnames = fieldnames(sweepResults);
for f = 1:length(fnames)
    sweepInfo.(fnames{f}) = vertcat(sweepResults.(fnames{f}));
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Now compare to what we observe %%%%%%%%%%%%%%%%%%%%%%
% close all
% % full distribution
% figure(1);
% hold on
% histogram([spot_struct.fluoInterp],'Normalization','Probability')
% histogram(sweepInfo.ms2_traces_true_wt,'Normalization','Probability')

% look at the lower piece
fluo_thresh = 1e5;
fluo_vec_raw = [spot_struct.fluo];
exp_fluo_vec = fluo_vec_raw(fluo_vec_raw<=fluo_thresh);
pd_fluo_vec = sweepInfo.ms2_traces_true_wt(sweepInfo.ms2_traces_true_wt<=fluo_thresh);


% figure(2);
fluo_bins = linspace(0,1e5,26);
fluo_axis = fluo_bins(1:end-1) + diff(fluo_bins)/2;

% hold on
% histogram(pd_fluo_vec,fluo_bins,'Normalization','Probability')
% histogram(exp_fluo_vec,fluo_bins,'Normalization','Probability')

%% let's try to estimate a curve indicating probability of a missed
% detection as a function of true spot fluorescence
% should be 100% when F=0 and 0% (more or less) when F=1e5

exp_counts = histcounts(exp_fluo_vec,fluo_bins,'Normalization','Probability');
pd_counts = histcounts(pd_fluo_vec,fluo_bins,'Normalization','Probability');

% readjust probabilities
exp_counts_norm = exp_counts .* pd_counts(end) / exp_counts(end);
p_missed_detection = (pd_counts - exp_counts_norm)./pd_counts;

% extend to reinforce need to congerge to zero
p_missed_detection_long = [p_missed_detection zeros(size(p_missed_detection))];
fluo_axis_long = [fluo_axis fluo_axis+fluo_axis(end)];

% let's fit a simple hill curve to this
h_fun = @(x) x(2)^x(1) ./ (x(2)^x(1) + fluo_axis_long.^x(1));
ob_fun = @(x) p_missed_detection_long - h_fun(x);
fit = lsqnonlin(ob_fun,[1 1e4],[-Inf -Inf],[Inf Inf]);

% generate predicted miss probability curve
fluo_ref_curve = linspace(0,2*nanmax(fluo_vec_raw),1e3);
p_miss_ref_vec = fit(2)^fit(1) ./ (fit(2)^fit(1) + fluo_ref_curve.^fit(1));

