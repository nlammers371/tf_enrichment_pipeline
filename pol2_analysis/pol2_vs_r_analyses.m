clear
close all 

projectName = 'Rbp1-GFP_eveBAC-mCh';

% get paths 
liveProject = LiveEnrichmentProject(projectName);
resultsRoot = [liveProject.dataPath filesep];
resultsDir = [resultsRoot 'cpHMM_results' filesep];

% load data
load([resultsDir filesep 'hmm_input_output.mat'])

% make save path for figures
TracePath = [liveProject.figurePath filesep 'pol2_vs_r' filesep 'traces' filesep];
mkdir(TracePath);
TraceDeltaPath = [liveProject.figurePath filesep 'pol2_vs_r' filesep 'traces_delta' filesep];
mkdir(TraceDeltaPath);
  
%% infer mCherry x GFP calibration factor

% generate inference vectors
f_vec_long = vertcat(hmm_input_output(9).fluo_pd_full)';
% pol2_vec_spot_long = [hmm_input_output.spot_protein_vec];

for i = 1:length(hmm_input_output)  
    hmm_input_output(i).ctrl_protein_sm = imgaussfilt(hmm_input_output(i).edge_null_protein_vec,7);
    hmm_input_output(i).spot_protein_dt = hmm_input_output(i).spot_protein_vec - hmm_input_output(i).ctrl_protein_sm;
end

pol2_vec_ctrl_long = [hmm_input_output.ctrl_protein_sm];
time_vec_long = [hmm_input_output.time];

pol2_vec_spot_long = [hmm_input_output(9).spot_protein_dt];
% fit
objective_fun = @(x) x(1)*f_vec_long - pol2_vec_spot_long;
x = lsqnonlin(objective_fun,[0.1],[0],[Inf]);
mdl = fitlm([f_vec_long'] ,pol2_vec_spot_long');
% x = [];
x(1)% = mdl.Coefficients.Estimate(1)
% x(2) = mdl.Coefficients.Estimate(2)
% x(3) = mdl.Coefficients.Estimate(3)

%% generate predicted pol II vector
% pol2_vec_spot_pd_long = x(2)*f_vec_long;

for i = 1:length(hmm_input_output)  
    hmm_input_output(i).pol2_predicted_vec =  1.5*hmm_input_output(i).fluo_pd_full'*x(1);
    hmm_input_output(i).pol2_delta = hmm_input_output(i).spot_protein_dt - hmm_input_output(i).pol2_predicted_vec;
end

%%
plot_id = 9;
close all
trace_fig = figure;
hold on
p1 = plot(hmm_input_output(plot_id).time/60,hmm_input_output(plot_id).pol2_predicted_vec,'-r','LineWidth',1.5);

ylabel('Pol II at the locus (au)')
yyaxis right
% p3 = plot(hmm_input_output(plot_id).time/60,hmm_input_output(plot_id).fluo,'-k','LineWidth',1);
p2 = plot(hmm_input_output(plot_id).time/60,hmm_input_output(plot_id).spot_protein_dt,'-b','LineWidth',1.5);
xlabel('time (minutes)')
ylabel('MS2 spot intensity (au)')
grid on
% legend([p2 p1 p3],'measured pol II at locus','predicted elongating pol II','MS2 trace','location','southeast')
saveas(trace_fig,[TraceDeltaPath 'example_trace_09.png'])

figure;
plot(hmm_input_output(plot_id).time/60,hmm_input_output(plot_id).pol2_delta,'-b','LineWidth',1.5);
yyaxis right
p2 = plot(hmm_input_output(plot_id).time/60,hmm_input_output(plot_id).r_vec_soft,'-k','LineWidth',1.5);
%%
% 9
for plot_id = 1:length(hmm_input_output)
    
    close all
    trace_fig = figure('Visible','off');
    hold on
    plot(hmm_input_output(plot_id).time/60,hmm_input_output(plot_id).spot_protein_vec)
    plot(hmm_input_output(plot_id).time/60,hmm_input_output(plot_id).pol2_predicted_vec)
    
    xlabel('time (minutes)')
    ylabel('Pol 2 intensity (au)')
    
    legend('local Pol 2','predicted Pol 2')
    grid on
    
    saveas(trace_fig,[TracePath 'trace_' sprintf('%03d',plot_id) '.png'])
end