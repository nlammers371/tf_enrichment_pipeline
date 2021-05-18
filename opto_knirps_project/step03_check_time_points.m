% Script to incorporate curvature-adjusted AP positions 
% Run this after running main02
clear
close all

projectList = {'optokni_eve4+6_WT'}; % Cell array containing all projects you wish to process
%projectList = {'optokni_eve4+6_MCP-GFP_Homo'}; % Cell array containing all projects you wish to process

%apBins = linspace(-0.12,0.12,8);
apBins = linspace(-0.08,0.12,7);
%apBins = linspace(55,67.5,6);
%timeRange = [7.5*60 37.5*60];
%timeRange = {[0 12.5*60],[12.5*60 25*60]};

% iterate through list of projects and generate required variables
for p = 1:length(projectList)
  
  % get liveProject object
  liveProject = LiveEnrichmentProject(projectList{p});
  
  % load compiled traces dataot
  disp('loading spot structure...')
  
  % load all the data
  %load([liveProject.dataPath filesep 'spot_struct_protein.mat']);
  load([liveProject.dataPath filesep 'spot_struct.mat']);
  analysis_traces = spot_struct;
  
  disp('calculating number of time points for each AP bin')
  % initialize grouping variables and add new position
  
  % apply QC filter calculated in main01
  analysis_traces = analysis_traces([analysis_traces.TraceQCFlag]==1);
  
  trace_length = zeros(size(analysis_traces));
  ap_pos = zeros(size(analysis_traces));
  ap_binNum = zeros(size(analysis_traces));
  bin_len = zeros(length(apBins)-1,1);
  
  for i = 1:length(analysis_traces)
     
    time_interp = analysis_traces(i).timeInterp;
    ap_interp = analysis_traces(i).APPosParticleInterp;
    ap_pos(i) = mean(ap_interp);
    
    ap_interp = ap_interp((time_interp>=timeRange(1)) & (time_interp<=timeRange(2)));
    ap_binNum(i) = discretize(ap_pos(i),apBins);
    trace_length(i) = length(ap_interp);
    
    %update the time points for each ap bin 
    if ~isnan(ap_binNum(i))
        bin_len(ap_binNum(i)) = bin_len(ap_binNum(i)) + trace_length(i);
    end
 
  end
  
  N = histcounts(ap_pos,apBins);

end