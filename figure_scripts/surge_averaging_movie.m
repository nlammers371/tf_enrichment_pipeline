% Script to attempt a systematic dissection of various factors driving
% proteinxtranscription burst coincidence
clear
close all
addpath('../utilities')
% define core ID variables
project = 'Dl-Ven_snaBAC-mCh';
% project = 'Dl-Ven_hbP2P-mCh';
DropboxFolder =  'E:\Nick\LivemRNA\Dropbox (Personal)\';
[~, DataPath, FigRoot] =   header_function(DropboxFolder, project);
FigPath = [FigRoot '\' project '\surge_averaging_movie\snips\'];
mkdir(FigPath)
% load data
load([DataPath 'hmm_input_output_results.mat'])

% define size of window of interest
roi_window = 6; 
window_size = 15;
start = window_size + 2;
% extract roi_vectors from wapo and locus arrays
hmm_array = results_struct.hmm_array;
spot_array_dt = results_struct.spot_array_dt;
virtual_array_dt = results_struct.virtual_array_dt;

feature_sign_vec = results_struct.feature_sign_vec';
lag_size_vec = results_struct.lag_size_vec';
lead_size_vec = results_struct.lead_size_vec';
lag_dur_vec = results_struct.lag_dur_vec';
lead_dur_vec = results_struct.lead_dur_vec';
tr_burst_size_vec = lag_dur_vec.*lag_size_vec;
% make rise filter
min_pause_len = 6;
min_burst_len = 2;
burst_ft_primary = feature_sign_vec == 1&lead_dur_vec>=min_pause_len&lag_dur_vec>min_burst_len; % filter for rise events
min_len_ft = feature_sign_vec == 1&lead_dur_vec==min_pause_len&lag_dur_vec>min_burst_len; % filter for rise events
%%% Get list of indices to use for plot
use_indices = find(burst_ft_primary)';
start_indices = find(min_len_ft)';

start_index = start_indices(7);

% set parametrs for plotting
N = 1000;
increment = 10;
n_iters = N/increment;
Tres = 20;
%  determine snip size
n_col = size(spot_array_dt,2);
window_size = floor(n_col/2);
time_axis = (-window_size:window_size)*Tres/60;

% seed random number generator for consistency
rng(312);
% shuffle indices to show
plot_indices = randsample(use_indices(use_indices~=start_index),N,false);
plot_indices = [start_index plot_indices];
close all
% iterate
for n = 1:n_iters+1
    pt_indices = plot_indices(1:max(1,(n-1)*increment));
    hmm_mean = nanmean(hmm_array(pt_indices,:),1);
    spot_pt_mean = nanmean(spot_array_dt(pt_indices,:),1);
    
    % make figure
    burst_fig = figure('Visible','off');
    cmap1 = brewermap([],'Set2');
    % snail activity
    yyaxis right
%     p1 = area(time_axis,hmm_mean,'FaceColor',cmap1(end,:),'LineWidth',1.5,'FaceAlpha',.4);
    p1 = plot(time_axis,hmm_mean,'--','LineWidth',2,'Color','black');
    ylabel('snail transcription (au)')
    set(gca,'ytick',0:.2:1.4)
    if n == 1
        text(2.5,1.3,'1 sample ','Fontsize',14)
    elseif n < 8
        text(2.5,1.3,[num2str((n-1)*increment) ' samples'],'Fontsize',14)
    else
        text(2.5,1.2-.0714,[num2str((n-1)*increment) ' samples'],'Fontsize',14)
    end
    if n >=8
        ylim([.2 1.2])
    else
        ylim([0 1.4])
    end
    ax = gca;
    ax.YColor = 'black';
    
    % Dorsal activity
    yyaxis left
    hold on
    % fill([time_axis fliplr(time_axis)],[br_spot_ub fliplr(br_spot_lb)],cmap1(2,:),'FaceAlpha',.5,'EdgeAlpha',0)
    p2 = plot(time_axis,spot_pt_mean,'-','Color',[213 108 85]/256,'LineWidth',2);
    p = plot(0,0);
    ylabel('relative Dl enrichment (au)')    
    
%     ax = gca;
    ax.YColor = [213 108 85]/256;   
%     StandardFigurePBoC(p,gca);
    if n < 8
        y_lim = ylim;
        ylim([y_lim(1) -5/4 * y_lim(1)])    
    else
        ylim([-20 25])    
    end
    xlabel('time from burst start (minutes)')
    lgd = legend([p1 p2],'snail transcription','Dl concentration','Location','southeast');
    set(lgd,'color',[228,221,209]/255);
    set(gca,'Fontsize',14,'xtick',-4:2:4)
    set(gca,...
            'Box','off',...
            'Color',[228,221,209]/255,...          
            'TickLength',[0.02,0.05])    
    burst_fig.InvertHardcopy = 'off';        
    % save
    saveas(burst_fig,[FigPath 'surge_averaging_n' sprintf('%04d',1+(n-1)*increment) '.tif'])
end   