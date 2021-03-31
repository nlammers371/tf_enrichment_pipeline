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
FigPath = [liveProject.figurePath filesep 'pol2_vs_r' filesep];
mkdir(FigPath)
TracePath = [liveProject.figurePath filesep 'pol2_vs_r' filesep 'traces' filesep];
mkdir(TracePath);
TraceDeltaPath = [liveProject.figurePath filesep 'pol2_vs_r' filesep 'traces_delta' filesep];
mkdir(TraceDeltaPath);
  
% Manually screen for "good" traces
good_ids = [1 4 5 7 8 9 10 12 13 15 16 19 21 23 25 26 27 28 32 33 34 36 39  40 43 45 47 48 49  51];

% for ind = 49:length(hmm_input_output)
%     close all
%     figure;   
%     plot(hmm_input_output(ind).time,hmm_input_output(ind).fluo)
%     title(['trace number ' num2str(ind)])
%     pause(3)
% end

hmm_input_output_full = hmm_input_output;
hmm_input_output = hmm_input_output(good_ids);

%% infer mCherry x GFP calibration factor

% generate inference vectors
f_vec_long = vertcat(hmm_input_output.fluo_pd_full)';
pol2_vec_spot_long = [hmm_input_output.spot_protein_vec];

for i = 1:length(hmm_input_output)
    hmm_input_output(i).ctrl_protein_sm = imgaussfilt(hmm_input_output(i).edge_null_protein_vec,7);
    hmm_input_output(i).spot_protein_vec_rel = hmm_input_output(i).spot_protein_vec-hmm_input_output(i).ctrl_protein_sm;
end

% pol2_vec_ctrl_long = [hmm_input_output.ctrl_protein_sm];
pol2_vec_spot_rel = [hmm_input_output.spot_protein_vec_rel];
% time_vec_long = [hmm_input_output(good_ids).time];
nan_filter = ~isnan(pol2_vec_spot_rel);
% fit
objective_fun = @(x) x(1)*f_vec_long(nan_filter) - pol2_vec_spot_rel(nan_filter);
x = lsqnonlin(objective_fun,[2],[0],[Inf]);
% mdl = fitlm([pol2_vec_ctrl_long' f_vec_long'] ,pol2_vec_spot_long');
x(1)% = mdl.Coefficients.Estimate(1);
% x(2) = mdl.Coefficients.Estimate(2);
% x(3) = mdl.Coefficients.Estimate(3)

%% generate predicted pol II vector
% pol2_vec_spot_pd_long = x(2)*f_vec_long;

for i = 1:length(hmm_input_output)
    hmm_input_output(i).pol2_predicted_vec = hmm_input_output(i).fluo_pd_full'*x(1);
    hmm_input_output(i).pol2_delta = hmm_input_output(i).spot_protein_vec - hmm_input_output(i).pol2_predicted_vec;
end
%% Make figures illustrating the process for inferring elongating pol II
ex_id = 3;
close all

trace_fig1 = figure;
cmap = brewermap([],'Set2');

p3 = plot(hmm_input_output(ex_id).time/60,hmm_input_output(ex_id).fluo,'-k','LineWidth',2);
   
xlabel('time (minutes)');
ylabel('MS2 spot intensity')
grid on
set(gca,'FontSize',14)
legend('raw trace');
set(gca,'Color',[228,221,209]/255) 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

trace_fig1 .InvertHardcopy = 'off';
set(gcf,'color','w');
xlim([2 14])
saveas(trace_fig1 ,[FigPath 'trace_fig1.png'])
saveas(trace_fig1 ,[FigPath 'trace_fig1.pdf'])


trace_fig2 = figure;
hold on
plot(hmm_input_output(ex_id).time/60,hmm_input_output(ex_id).fluo,'-k','LineWidth',2);
plot(hmm_input_output(ex_id).time/60,hmm_input_output(ex_id).fluo_pd_MS2,'--','Color',cmap(2,:),'LineWidth',2);
   
xlabel('time (minutes)');
ylabel('MS2 spot intensity')
legend('raw trace','model fit');
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
xlim([2 14])
ax = gca;
% ax.YAxis(1).Color = cmap(5,:);
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

trace_fig2.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(trace_fig2 ,[FigPath 'trace_fig2.png'])
saveas(trace_fig2 ,[FigPath 'trace_fig2.pdf'])

trace_fig3 = figure;

hold on
plot(hmm_input_output(ex_id).time/60,hmm_input_output(ex_id).fluo,'-k','LineWidth',2);
plot(hmm_input_output(ex_id).time/60,hmm_input_output(ex_id).fluo_pd_MS2,'--','Color',cmap(2,:),'LineWidth',2);
plot(hmm_input_output(ex_id).time/60,hmm_input_output(ex_id).fluo_pd_full,'Color',cmap(3,:),'LineWidth',2);
   
xlabel('time (minutes)');
ylabel('MS2 spot intensity')
legend('raw trace','model fit','predicted Pol II');
grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
xlim([2 14])
ax = gca;
% ax.YAxis(1).Color = cmap(5,:);
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

trace_fig3.InvertHardcopy = 'off';
set(gcf,'color','w');
% trace_fig4.InvertHardcopy = 'off';
saveas(trace_fig3 ,[FigPath 'trace_fig3.png'])
saveas(trace_fig3 ,[FigPath 'trace_fig3.pdf'])


trace_fig4 = figure;
hold on
% p1 = plot(hmm_input_output(ex_id).time/60,hmm_input_output(ex_id).pol2_predicted_vec,'Color',cmap(3,:),'LineWidth',2);
p2 = plot(hmm_input_output(ex_id).time/60,hmm_input_output(ex_id).spot_protein_vec_rel,'Color',cmap(5,:),'LineWidth',2);
ylabel('Pol II at the locus (au)')

yyaxis right
plot(hmm_input_output(ex_id).time/60,hmm_input_output(ex_id).fluo,'-k','LineWidth',1);
% plot(hmm_input_output(plot_id).time/60,hmm_input_output(plot_id).fluo_pd_full,'-o','LineWidth',1.5);
legend('predicted Pol II','raw trace');
xlabel('time (minutes)')
ylabel('MS2 spot intensity (au)')

grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
xlim([2 14])
ax = gca;
ax.YAxis(1).Color = cmap(5,:);
ax.YAxis(2).Color = 'k';
ax.XAxis(1).Color = 'k';
trace_fig4.InvertHardcopy = 'off';
set(gcf,'color','w');
% legend([p2 p1 p3],'measured pol II at locus','predicted elongating pol II','MS2 trace','location','southeast')
saveas(trace_fig4,[FigPath 'trace_fig4.png'])

trace_fig5 = figure;
hold on
p1 = plot(hmm_input_output(ex_id).time/60,hmm_input_output(ex_id).pol2_predicted_vec,'Color',cmap(3,:),'LineWidth',2);
p2 = plot(hmm_input_output(ex_id).time/60,hmm_input_output(ex_id).spot_protein_vec_rel,'Color',cmap(5,:),'LineWidth',2);
ylabel('Pol II at the locus (au)')

yyaxis right
plot(hmm_input_output(ex_id).time/60,hmm_input_output(ex_id).fluo,'-k','LineWidth',1);
% plot(hmm_input_output(plot_id).time/60,hmm_input_output(plot_id).fluo_pd_full,'-o','LineWidth',1.5);
legend('predicted Pol II','measured Pol II','raw trace');
xlabel('time (minutes)')
ylabel('MS2 spot intensity (au)')

grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
xlim([2 14])
ax = gca;
ax.YAxis(1).Color = cmap(5,:);
ax.YAxis(2).Color = 'k';
ax.XAxis(1).Color = 'k';

trace_fig5.InvertHardcopy = 'off';
set(gcf,'color','w');
% legend([p2 p1 p3],'measured pol II at locus','predicted elongating pol II','MS2 trace','location','southeast')
saveas(trace_fig5,[FigPath 'trace_fig5.png'])

%%
%6,12: nice bursty trace
%10: large unpredicted pol II surge
%11,14: well-predicted
for ex_id = 1:length(hmm_input_output)
% plot_id = good_ids(16);
  trace_fig5 = figure('Visible','off');
  hold on
  p1 = plot(hmm_input_output(ex_id).time/60,hmm_input_output(ex_id).pol2_predicted_vec,'Color',cmap(3,:),'LineWidth',2);
  p2 = plot(hmm_input_output(ex_id).time/60,hmm_input_output(ex_id).spot_protein_vec_rel,'Color',cmap(5,:),'LineWidth',2);
  ylabel('Pol II at the locus (au)')

  yyaxis right
  plot(hmm_input_output(ex_id).time/60,hmm_input_output(ex_id).fluo,'-k','LineWidth',1);
  % plot(hmm_input_output(plot_id).time/60,hmm_input_output(plot_id).fluo_pd_full,'-o','LineWidth',1.5);
  legend('predicted Pol II','measured Pol II','raw trace');
  xlabel('time (minutes)')
  ylabel('MS2 spot intensity (au)')

  grid on
  set(gca,'FontSize',14)
  set(gca,'Color',[228,221,209]/255) 
%   xlim([2 14])
  ax = gca;
  ax.YAxis(1).Color = cmap(5,:);
  ax.YAxis(2).Color = 'k';
  ax.XAxis(1).Color = 'k';
  trace_fig5.InvertHardcopy = 'off';
  set(gcf,'color','w');
  % legend([p2 p1 p3],'measured pol II at locus','predicted elongating pol II','MS2 trace','location','southeast')
  saveas(trace_fig5,[TracePath 'trace_fig_' sprintf('%03d',ex_id) '.png'])
end
%%
bad_ids = find(~ismember(1:length(hmm_input_output_full),good_ids));
% 9
for ex_id = bad_ids

    bad_trace_fig = figure('Visible','off');
    cmap = brewermap([],'Set2');

    p3 = plot(hmm_input_output_full(ex_id).time/60,hmm_input_output_full(ex_id).fluo,'-k','LineWidth',2);

    xlabel('time (minutes)');
    ylabel('MS2 spot intensity')
    grid on
    set(gca,'FontSize',14)
    legend('raw trace');
    set(gca,'Color',[228,221,209]/255) 

    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.XAxis(1).Color = 'k';

    bad_trace_fig.InvertHardcopy = 'off';
    set(gcf,'color','w');
    % xlim([2 14])
    saveas(bad_trace_fig,[FigPath 'bad_trace_' sprintf('%03d',ex_id) '.png'])
end
% saveas(bad_trace_fig,[FigPath 'bad_trace_fig.pdf'])