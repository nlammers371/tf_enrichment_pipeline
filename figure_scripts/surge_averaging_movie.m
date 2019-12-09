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
single_frame_flag = false;
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
n_nan_ft = sum(~isnan(spot_array_dt),2)>=2*window_size;
%%% Get list of indices to use for plot
use_indices = find(burst_ft_primary)';
single_options = find(burst_ft_primary&n_nan_ft)';
start_indices = find(min_len_ft)';

start_index = start_indices(7);

% set parametrs for plotting
N = 1000;
increment = 10;
n_iters = N/increment;
Tres = 20;
n_start = 10;
%  determine snip size
n_col = size(spot_array_dt,2);
window_size = floor(n_col/2);
time_axis = (-window_size:window_size)*Tres/60;

% seed random number generator for consistency
rng(312);
% shuffle indices to show
plot_indices = randsample(use_indices(use_indices~=start_index),N,false);
plot_indices = [start_index plot_indices];
% single_indices = randsample(single_options,n_start,false);
single_indices = [11518 11656 10935 10651];
close all
if single_frame_flag
    % make series of single snips 
    for n = 1:n_start
        pt_index = single_indices(n);
        hmm_mean = nanmean(hmm_array(pt_index,:),1);
        spot_pt_mean = nanmean(spot_array_dt(pt_index,:),1);

        % make figure
        burst_fig = figure('Visible','off');
        cmap1 = brewermap([],'Set2');
        % snail activity
        yyaxis right
        p1 = plot(time_axis,hmm_mean,'--','LineWidth',2,'Color','black');
        ylabel('snail transcription (au)')    
        ylb = get(gca,'ylabel');    
        set(gca,'ytick',0:.2:1.4)
        ylim([0 1.4])
        ax = gca;
        ax.YColor = 'black';

        % Dorsal activity
        yyaxis left
        hold on
        % fill([time_axis fliplr(time_axis)],[br_spot_ub fliplr(br_spot_lb)],cmap1(2,:),'FaceAlpha',.5,'EdgeAlpha',0)
        p2 = plot(time_axis,spot_pt_mean,'-','Color',[213 108 85]/256,'LineWidth',2);
        p = plot(0,0);
        ylabel('relative Dorsal enrichment (au)');       
        ax.YColor = [213 108 85]/256;      
        ylim([-200 5/4*200])   
        xlabel('time from burst start (minutes)')
        set(gca,'Fontsize',14,'xtick',-4:2:4)
        set(gca,...
                'Box','on',...
                'Color',[228,221,209]/255,...          
                'TickLength',[0.02,0.05])    
        burst_fig.Color = 'white';
        burst_fig.InvertHardcopy = 'off'; 
        % rotate right y label
        ylb.Position(1) = ylb.Position(1)+.7;
        set(ylb,'rotation',-90)

        % save
        saveas(burst_fig,[FigPath 'single_frame_ind' sprintf('%04d',pt_index) '.tif'])
    end   
end

close all
% make averaging movie
for n = 1:n_iters+1
    pt_indices = plot_indices(1:max(1,(n-1)*increment));
    hmm_mean = nanmean(hmm_array(pt_indices,:),1);
    spot_pt_mean = nanmean(spot_array_dt(pt_indices,:),1);
    
    % make figure
    burst_fig = figure('Visible','off');
    cmap1 = brewermap([],'Set2');
    % snail activity
    yyaxis right
    p1 = plot(time_axis,hmm_mean,'--','LineWidth',2,'Color','black');
    ylabel('snail transcription (au)')    
    ylb = get(gca,'ylabel');    
    set(gca,'ytick',0:.2:1.4)
    if n == 1
        text(2.5,1.3,'1 burst','Fontsize',14)
    else
        text(2.5,1.3,[num2str((n-1)*increment) ' bursts'],'Fontsize',14)
    end
      ylim([0 1.4])
    ax = gca;
    ax.YColor = 'black';
    
    % Dorsal activity
    yyaxis left
    hold on
    % fill([time_axis fliplr(time_axis)],[br_spot_ub fliplr(br_spot_lb)],cmap1(2,:),'FaceAlpha',.5,'EdgeAlpha',0)
    p2 = plot(time_axis,spot_pt_mean,'-','Color',[213 108 85]/256,'LineWidth',2);
    p = plot(0,0);
    ylabel('relative Dorsal enrichment (au)');       
    ax.YColor = [213 108 85]/256;   
    if n < 20
        y_lim = ylim;
        if -5/4 * y_lim(1) >= y_lim(2)
            ylim([y_lim(1) -5/4 * y_lim(1)])    
        else
            ylim([-4/5*y_lim(2) y_lim(2)])    
        end
    else
        ylim([-20 25])    
    end
    xlabel('time from burst start (minutes)')
    set(gca,'Fontsize',14,'xtick',-4:2:4)
    set(gca,'Box','on',...
            'Color',[228,221,209]/255,...          
            'TickLength',[0.02,0.05])    
    burst_fig.Color = 'white';
    burst_fig.InvertHardcopy = 'off'; 
    % rotate right y label
    ylb.Position(1) = ylb.Position(1)+.7;
    set(ylb,'rotation',-90)
    
    % save
    saveas(burst_fig,[FigPath 'surge_averaging_n' sprintf('%04d',1+(n-1)*increment) '.tif'])
end   