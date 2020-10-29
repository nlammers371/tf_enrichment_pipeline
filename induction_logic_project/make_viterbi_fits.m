% Script to conduct 
clear
close all

% set path to utilities
addpath('utilities')

% set project
project = '20200807'; 

% make figure path 
n_plot = 10;
figurePath = ['../fig/' project '/ example_traces/'];
mkdir(figurePath);

% load trace data
dataPath = ['../dat/' project '/']; %Path to inference data
dataname = 'trace_structure.mat'; %name of inference set
load([dataPath dataname]);

% parameter values
K = 2;
w = 2;
alpha = trace_structure(1).alpha_frac*w; 
groups_to_run = 1;

% set path to inference results
infDir =  ['../out/' project '/w' num2str(w) '/K' num2str(K) '/']; 

% compile inference results
fileList = dir([infDir 'raw_cpHMM_output/*.mat']);
inference_struct = struct;
for f = 1:length(fileList)
  load([infDir fileList(f).name])
  fnames = fieldnames(output);
  for n = 1:length(fnames)
    inference_struct(f).(fnames{n}) = output.(fnames{n});
  end
end

A_mean = mean(cat(3,inference_struct.A_mat),3);
A_mean = A_mean ./ sum(A_mean);

r_mean = median([inference_struct.r],2);
v_mean = r_mean*trace_structure(1).Tres;

pi0_mean = mean(vertcat(inference_struct.pi0),1);
pi0_mean = pi0_mean/sum(pi0_mean);

sigma_mean = sqrt(mean([inference_struct.noise]));

fit_indices = find(ismember([trace_structure.setID],groups_to_run));

v_fits = struct;
for f = 1:length(fit_indices)
  viterbi_out = viterbi (trace_structure(fit_indices(f)).fluo, v_mean, sigma_mean, ...
    log(pi0_mean),log(A_mean), K, w, alpha);
  fnames = fieldnames(viterbi_out);
  for j = 1:numel(fnames)
      v_fits(f).(fnames{j}) = viterbi_out.(fnames{j});
  end
end

% make output data table
for f = 1:length(fit_indices)
  if f == 1
    longTable = [f*ones(size(trace_structure(fit_indices(f)).time')),...
      trace_structure(fit_indices(f)).time', trace_structure(fit_indices(f)).fluo',...
        double(v_fits(f).z_viterbi'),v_fits(f).fluo_viterbi'];
  else
    longTable = vertcat(longTable,...
      [f*ones(size(trace_structure(fit_indices(f)).time')),...
      trace_structure(fit_indices(f)).time', trace_structure(fit_indices(f)).fluo',...
        double(v_fits(f).z_viterbi'),v_fits(f).fluo_viterbi']);
  end
end

fitTable = array2table(longTable,'VariableNames',{'column_id','time','raw_fluorescence','promoter_state','predicted_fluorescence'});
writetable(fitTable,[infDir 'single_trace_fits.csv'])

%% plot example traces
close all
time_vec = trace_structure(fit_indices(1)).time;
rng(123);
plot_indices = randsample(1:length(fit_indices),n_plot,false);
for p = plot_indices
  fluo_fig = figure;
  cmap = brewermap([],'Set2');
  hold on
  plot(time_vec/60,trace_structure(fit_indices(p)).fluo,'Color','k','LineWidth',1.5)
  plot(time_vec/60,v_fits(p).fluo_viterbi,'Color',cmap(2,:),'LineWidth',1.5)
%   e = errorbar(time_vec/60,v_fits(p).fluo_viterbi,repelem(sigma_mean,length(time_vec)),'Color',cmap(2,:));
%   e.CapSize = 0;
  grid on
  set(gca,'Fontsize',14);
  xlabel('time (minutes)')
  ylabel('fluorescence intensity (au)')
  saveas(fluo_fig,[figurePath 'trace' sprintf('%03d',p) '_fluo_fit.png'])
  
%   z_fig = figure;
%   cmap = brewermap([],'Set2');
%   hold on
%   plot(time_vec/60,trace_structure(fit_indices(p)).fluo,'Color','k','LineWidth',1.5)
%   ylabel('fluorescence intensity (au)')
%   
%   yyaxis right
%   stairs(time_vec/60,v_fits(p).z_viterbi,'Color',cmap(3,:),'LineWidth',1.5)
%   ylabel('promoter state')
%   ylim([1 2.5])
%   grid on
%   set(gca,'Fontsize',14);
%   xlabel('time (minutes)')
%   
%   saveas(z_fig,[figurePath 'trace' sprintf('%03d',p) '_promoter_state_fit.png'])
end 
  
  
  