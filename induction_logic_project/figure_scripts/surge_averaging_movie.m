% Script to attempt a systematic dissection of various factors driving
% proteinxtranscription burst coincidence
clear
close all
addpath(genpath('../utilities'))

% define project
project = '2xDl-Ven_snaBAC-mCh';

% set paths
DropboxFolder =  'S:\Nick\Dropbox\';
[~, DataPath, FigRoot] =   header_function(DropboxFolder, project);
FigPath = [FigRoot '\' project '\surge_averaging_movie\snips\'];

mkdir(FigPath)

% HMM and sampling parameters
w = 7;
K = 3;
fluo_dim = 2;
protein_dim = 2;

% event filtering parameters
min_pause_len = 4; % minimum length of preceding OFF period (in time steps)
max_pause_len = 100;
min_burst_len = 2;
max_burst_len = 1000;

% averaging movie parameters
plot_single_frame_flag = false;
pdf_indices = [1 51 101 501]; 
roi_window = 5; 
window_size = 15;
start = window_size + 2;

% set parametrs for plotting
N = 1000; 
increment = 10;
Tres = 20;

% load data
load([DataPath 'hmm_input_output_results_w' num2str(w) '_K' num2str(K) '_f' ...
                      num2str(fluo_dim) 'D_p' num2str(protein_dim) 'D.mat'])

%% %%%%%%%%%%%%%%%%%%%%%% Build surge data structures %%%%%%%%%%%%%%%%%%%%%

% extract roi_vectors from wapo and locus arrays
hmm_array = results_struct.hmm_array;
spot_array_dt = results_struct.spot_array_dt;
% virtual_array_dt = results_struct.virtual_array_dt;
feature_sign_vec = results_struct.feature_sign_vec';
lag_dur_vec = results_struct.lag_dur_vec';
lead_dur_vec = results_struct.lead_dur_vec';

% make rise filter
burst_filter = feature_sign_vec == 1&lead_dur_vec>=min_pause_len&lag_dur_vec>min_burst_len; % filter for rise events
n_nan_ft = sum(~isnan(spot_array_dt),2)>=2*window_size;

N = min([N sum(burst_filter)]);
n_iters = round(N/increment);

%%% Get list of indices to use for plot
use_indices = find(burst_filter)';
single_options = find(burst_filter&n_nan_ft)';
start_indices = find(burst_filter)';
start_index = start_indices(7);


%  determine snip size
n_col = size(spot_array_dt,2);
window_size = floor(n_col/2);
time_axis = (-window_size:window_size)*Tres/60;


%% %%%%%%%%%%%%%%%%%%% Perform plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% seed random number generator for consistency
rng(312);

% shuffle indices to show
plot_indices = randsample(use_indices(use_indices~=start_index),N-1,false);
plot_indices = [start_indices(10) plot_indices];
% single_indices = randsample(single_options,n_start,false);
single_indices = [11518 11656 10935 10651];

close all
if plot_single_frame_flag
    % make series of single snips 
    for n = 1:numel(single_indices)
        pt_index = single_indices(n);
        hmm_mean = nanmean(hmm_array(pt_index,:),1);
        spot_protein_mean = nanmean(spot_array_dt(pt_index,:),1);

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
        p2 = plot(time_axis,spot_protein_mean,'-','Color',[213 108 85]/256,'LineWidth',2);
        p = plot(0,0);
        ylabel('relative Dorsal enrichment (au)');       
        ax.YColor = [213 108 85]/256;      
        ylim([-200 5/4*200])   
        xlabel('time from burst start (minutes)')
        set(gca,'Fontsize',14,'xtick',-4:2:4)
        set(gca,...
                'Box','on',...
                'Color',[228 220 209]/255,...          
                'TickLength',[0.02,0.05])    
        burst_fig.Color = 'white';
        burst_fig.InvertHardcopy = 'off'; 
        % rotate right y label
        ylb.Position(1) = ylb.Position(1)+.7;
        set(ylb,'rotation',-90)

        % save
        saveas(burst_fig,[FigPath 'single_frame_ind' sprintf('%04d',pt_index) '.pdf'])        
    end   
end

% generate dynamic y limits
y1_vec = -.05 - (-.05+.3)/(n_iters-1).^2 .* (-n_iters+1:0).^2;
y2_vec = .07 - (.07-.3)/(n_iters-1).^2 .* (-n_iters+1:0).^2;

close all
iter = 1;
% make averaging movie
for n = 1:increment:n_iters*increment     

    plot_indices_curr = plot_indices(1:min([max(1,n), length(plot_indices)]));

    hmm_mean = nanmean(hmm_array(plot_indices_curr,:),1);
    spot_protein_mean = nanmean(spot_array_dt(plot_indices_curr,:),1);

    % make figure
    burst_fig = figure('Visible','off');
    cmap1 = brewermap([],'Set2');    
    
    % snail activity
    yyaxis right   
               
    if n == 1
        p1 = stairs(time_axis,hmm_mean,'--','LineWidth',2,'Color','black');
        text(2.5,1.3,'1 burst','Fontsize',14)        
    else
        p1 = plot(time_axis,hmm_mean,'--','LineWidth',2,'Color','black');
        text(2.5,1.3,[num2str((n-1)*increment) ' bursts'],'Fontsize',14)
    end
    
    ylabel('snail transcription (au)')    
    ylb = get(gca,'ylabel'); 
    ylim([-.05 1.05])
    
    ax = gca;
    ax.YColor = 'black';

    % Dorsal activity
    yyaxis left
    hold on
    % fill([time_axis fliplr(time_axis)],[br_spot_ub fliplr(br_spot_lb)],cmap1(2,:),'FaceAlpha',.5,'EdgeAlpha',0)
    p2 = plot(time_axis,spot_protein_mean,'-','Color',[213 108 85]/256,'LineWidth',2);
    p = plot(0,0);
    ylabel('relative Dorsal enrichment (au)');       
    ax.YColor = [213 108 85]/256;   
    if n == 1
        ylim([-1.5 1.5])  
    else
        ylim([y1_vec(iter) y2_vec(iter)])  
        iter = iter + 1;
    end
%     elseif n < 20*increment
%         ylim([-.25 .25])  
%     elseif n < 38*increment
%         ylim([-.05 .1])  
%     else
%         ylim([-.05 .07])
%     end
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

    xlim([time_axis(1) time_axis(end)])
    % save
    saveas(burst_fig,[FigPath 'surge_averaging_n' sprintf('%04d',1+(n-1)*increment) '.tif'])
    if ismember(n,pdf_indices)
        saveas(burst_fig,[FigPath 'surge_averaging_n' sprintf('%04d',1+(n-1)*increment) '.pdf'])
    end    
end   