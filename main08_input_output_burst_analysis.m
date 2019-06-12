% Script to analyze protein feature variation as a function of size and
% duration of transcription bursts
clear
close all
addpath('utilities')
% define core ID variables
project = 'Dl-Ven_snaBAC-mCh';
% project = 'Dl-Ven_hbP2P-mCh';
dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
dataPath = [dropboxFolder 'ProcessedEnrichmentData\' project '\'];
figPath = [dropboxFolder 'LocalEnrichmentFigures\' project '\burst_analyses\'];
mkdir(figPath)
% load data
load([dataPath 'hmm_input_output_results.mat'])
w = 7;
K = 3;
load([dataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'],'hmm_input_output')

%%% Make Illustrative burst calling figure
hm_cm = flipud(brewermap([],'RdYlBu'));
index = 129;
burst_fig = figure;
z_vec = hmm_input_output(129).z_vec > 1;
spot_vec = hmm_input_output(129).spot_protein;
time_vec = hmm_input_output(129).time;

yyaxis left
pp = plot(time_vec/60,spot_vec,'LineWidth',1.5,'Color',hm_cm(end-10,:));
ax = gca;
ax.YColor = 'black';
ylabel('Dorsal concentration (au)')

yyaxis right 
hold on
x = [time_vec/60 time_vec/60];
y = [z_vec' z_vec'];
s = area(x,y,'FaceColor',[0 0 0],'EdgeAlpha',0,'FaceAlpha',.3);
ylim([0 1.35])
ax = gca;
ax.YColor = 'black';
ylabel('activity state')
p = plot(0,0);
xlabel('minutes')
xlim([10 40])
legend([pp s],'protein level','transcriptional activity');
StandardFigure(p,gca)
saveas(burst_fig, [figPath 'burst_example.tif'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% uncontrolled burst duration figures

% burst rises rise
if strcmp(project,'Dl-Ven_hbP2P-mCh')
    burst_range = 5:10;
    burst_window = 4;
else
    burst_range = 3:12;
    burst_window = 3;
end
min_buffer_len = 6;
% extract relevant arrays
lag_dur_vec = results_struct.lag_dur_vec;
lead_dur_vec = results_struct.lead_dur_vec;
lead_sz_vec = results_struct.lead_size_vec;
feature_sign_vec = results_struct.feature_sign_vec;
hmm_array = results_struct.hmm_array;
spot_array = results_struct.spot_array;
swap_array = results_struct.swap_array;
window_size = size(swap_array,2);
% initialize data arrays
burst_rise_dur_hmm_mean = NaN(numel(burst_range),window_size);
burst_rise_dur_spot_mean = NaN(numel(burst_range),window_size);
burst_rise_dur_swap_mean = NaN(numel(burst_range),window_size);
for i = 1:numel(burst_range)
    burst_vec = burst_range(i)-burst_window:burst_range(i)+burst_window;
    burst_ft = feature_sign_vec == 1 & ismember(lag_dur_vec,burst_vec) & lead_dur_vec >= min_buffer_len;
    % calculate averages
    burst_rise_dur_hmm_mean(i,:) = nanmean(hmm_array(burst_ft,:));
    burst_rise_dur_spot_mean(i,:) = nanmean(spot_array(burst_ft,:));
    burst_rise_dur_swap_mean(i,:) = nanmean(swap_array(burst_ft,:));
end

% window_vec = ((1:window_size) - ceil(window_size/2))/3;
cm_rise = brewermap(128,'Reds');
inc = floor(128/numel(burst_range));

burst_rise_dur_hm = figure;
colormap(hm_cm)
% colormap(jet(128)/1.1);
pcolor(flipud(burst_rise_dur_spot_mean))
xlabel('offset (min)')
set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)
ylabel('burst duration (min)')
set(gca,'yticklabels',fliplr(round(burst_range/3,1)))
h = colorbar;
ylabel(h, 'Dorsal levels (au)')
set(gca,'FontSize', 12);
saveas(burst_rise_dur_hm, [figPath 'burst_dur_rise_hm_target.tif'])

swap_rise_dur_hm = figure;
colormap(hm_cm)
% colormap(jet(128)/1.1);
pcolor(flipud(burst_rise_dur_swap_mean))
xlabel('offset (min)')
set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)
ylabel('burst duration (min)')
set(gca,'yticklabels',fliplr(round(burst_range/3,1)))
h = colorbar;
ylabel(h, 'Dorsal levels (au)')
set(gca,'FontSize', 12);
saveas(swap_rise_dur_hm, [figPath 'burst_dur_rise_hm_swap.tif'])

hmm_rise_dur_hm = figure;
colormap(hm_cm)
pcolor(flipud(burst_rise_dur_hmm_mean))
xlabel('offset (min)')
set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)
ylabel('burst duration (min)')
set(gca,'yticklabels',fliplr(round(burst_range/3,1)))
h = colorbar;
ylabel(h, 'sna activity (au)')
set(gca,'FontSize', 12);
saveas(hmm_rise_dur_hm, [figPath 'burst_dur_rise_hm_hmm.tif'])

burst_rise_dur_wt = figure;
index_vec = 1:numel(burst_range);
hold on
for i = 1:numel(burst_range)
    temp = burst_rise_dur_spot_mean;
    temp(index_vec~=i,:) = NaN;
    w = waterfall(temp-nanmin(temp(:)),repmat((1:numel(burst_range))',1,window_size));
    w.FaceColor = cm_rise(1+(i-1)*inc,:);
    w.FaceAlpha = .8;
    w.EdgeColor = 'black';
end
xlabel('offset (min)')
set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)    
zlabel('Dorsal levels (au)')    
view(-15,20)
grid on
saveas(burst_rise_dur_wt, [figPath 'burst_dur_rise_waterfall_target.tif'])

swap_rise_dur_wt = figure;
index_vec = 1:numel(burst_range);
hold on
for i = 1:numel(burst_range)
    temp = burst_rise_dur_swap_mean;
    temp(index_vec~=i,:) = NaN;
    w = waterfall(temp-nanmin(temp(:)),repmat((1:numel(burst_range))',1,window_size));
    w.FaceColor = cm_rise(1+(i-1)*inc,:);
    w.FaceAlpha = .8;
    w.EdgeColor = 'black';
end
xlabel('offset (min)')
set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)    
zlabel('Dorsal levels (au)')    
view(-15,20)
grid on
saveas(swap_rise_dur_wt, [figPath 'burst_dur_rise_waterfall_swap.tif'])

hmm_rise_dur_wt = figure;
index_vec = 1:numel(burst_range);
hold on
for i = 1:numel(burst_range)
    temp = burst_rise_dur_hmm_mean;
    temp(index_vec~=i,:) = NaN;
    w = waterfall(temp-nanmin(temp(:)),repmat((1:numel(burst_range))',1,window_size));
    w.FaceColor = cm_rise(1+(i-1)*inc,:);
    w.FaceAlpha = .8;
    w.EdgeColor = 'black';
end
xlabel('offset (min)')
set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)    
zlabel('sna activity (au)')    
view(-15,20)
grid on
saveas(hmm_rise_dur_wt, [figPath 'burst_dur_rise_waterfall_hmm.tif'])

if ~strcmp(project,'Dl-Ven_hbP2P-mCh')
    %%% burst falls
    % initialize data arrays
    burst_fall_dur_hmm_mean = NaN(numel(burst_range),window_size);
    burst_fall_dur_spot_mean = NaN(numel(burst_range),window_size);
    burst_fall_dur_swap_mean = NaN(numel(burst_range),window_size);
    for i = 1:numel(burst_range)
        burst_vec = burst_range(i)-burst_window:burst_range(i)+burst_window;
        burst_ft = feature_sign_vec == -1 & ismember(lead_dur_vec,burst_vec) & lag_dur_vec >= min_buffer_len;
        % calculate averages
        burst_fall_dur_hmm_mean(i,:) = nanmean(hmm_array(burst_ft,:));
        burst_fall_dur_spot_mean(i,:) = nanmean(spot_array(burst_ft,:));
        burst_fall_dur_swap_mean(i,:) = nanmean(swap_array(burst_ft,:));
    end

    window_vec = ((1:window_size) - ceil(window_size/2))/3;
    cm_fall = brewermap(128,'PuBu');
    inc = floor(128/numel(burst_range));

    burst_fall_dur_hm = figure;
    colormap(hm_cm)
    % colormap(jet(128)/1.1);
    pcolor(flipud(burst_fall_dur_spot_mean))
    xlabel('offset (min)')
    set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)
    ylabel('burst duration (min)')
    set(gca,'yticklabels',fliplr(round(burst_range/3,1)))
    h = colorbar;
    ylabel(h, 'Dorsal levels (au)')
    set(gca,'FontSize', 12);
    saveas(burst_fall_dur_hm, [figPath 'burst_dur_fall_hm_target.tif'])

    swap_fall_dur_hm = figure;
    colormap(hm_cm)
    pcolor(flipud(burst_fall_dur_swap_mean))
    xlabel('offset (min)')
    set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)
    ylabel('burst duration (min)')
    set(gca,'yticklabels',fliplr(round(burst_range/3,1)))
    h = colorbar;
    ylabel(h, 'Dorsal levels (au)')
    set(gca,'FontSize', 12);
    saveas(swap_fall_dur_hm, [figPath 'burst_dur_fall_hm_swap.tif'])

    hmm_fall_dur_hm = figure;
    colormap(hm_cm)
    pcolor(flipud(burst_fall_dur_hmm_mean))
    xlabel('offset (min)')
    set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)
    ylabel('burst duration (min)')
    set(gca,'yticklabels',fliplr(round(burst_range/3,1)))
    h = colorbar;
    ylabel(h, 'sna activity (au)')
    set(gca,'FontSize', 12);
    saveas(hmm_fall_dur_hm, [figPath 'burst_dur_fall_hm_hmm.tif'])

    burst_fall_dur_wt = figure;
    index_vec = 1:numel(burst_range);
    hold on
    for i = 1:numel(burst_range)
        temp = burst_fall_dur_spot_mean;
        temp(index_vec~=i,:) = NaN;
        w = waterfall(temp-nanmin(temp(:)),repmat((1:numel(burst_range))',1,window_size));
        w.FaceColor = cm_fall(1+(i-1)*inc,:);
        w.FaceAlpha = .8;
        w.EdgeColor = 'black';
    end
    xlabel('offset (min)')
    set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)    
    zlabel('Dorsal levels (au)')    
    view(15,20)
    grid on
    saveas(burst_fall_dur_wt, [figPath 'burst_dur_fall_waterfall_target.tif'])

    swap_fall_dur_wt = figure;
    index_vec = 1:numel(burst_range);
    hold on
    for i = 1:numel(burst_range)
        temp = burst_fall_dur_swap_mean;
        temp(index_vec~=i,:) = NaN;
        w = waterfall(temp-nanmin(temp(:)),repmat((1:numel(burst_range))',1,window_size));
        w.FaceColor = cm_fall(1+(i-1)*inc,:);
        w.FaceAlpha = .8;
        w.EdgeColor = 'black';
    end
    xlabel('offset (min)')
    set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)    
    zlabel('Dorsal levels (au)')    
    view(15,20)
    grid on
    saveas(swap_fall_dur_wt, [figPath 'burst_dur_fall_waterfall_swap.tif'])

    hmm_fall_dur_wt = figure;
    index_vec = 1:numel(burst_range);
    hold on
    for i = 1:numel(burst_range)
        temp = burst_fall_dur_hmm_mean;
        temp(index_vec~=i,:) = NaN;
        w = waterfall(temp-nanmin(temp(:)),repmat((1:numel(burst_range))',1,window_size));
        w.FaceColor = cm_fall(1+(i-1)*inc,:);
        w.FaceAlpha = .8;
        w.EdgeColor = 'black';
    end
    xlabel('offset (min)')
    set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)    
    zlabel('sna activity (au)')    
    view(15,20)
    grid on
    saveas(hmm_fall_dur_wt, [figPath 'burst_dur_fall_waterfall_hmm.tif'])

    %%% Now generate controlled burst amplitude and duration series
    lag_size_vec = results_struct.lag_size_vec;
    lead_size_vec = results_struct.lead_size_vec;
    % get rise-specific vectors
    burst_rise_size_vec = lag_size_vec;
    burst_rise_size_vec(feature_sign_vec~=1|lag_dur_vec < 2 | lead_dur_vec < min_buffer_len) = NaN;
    burst_rise_dur_vec = lag_dur_vec;
    burst_rise_dur_vec(feature_sign_vec~=1|lag_dur_vec < 2 | lead_dur_vec < min_buffer_len) = NaN;
    % get empirical distributions for each variable
    n_bins = 5;
    lb_dur = prctile(burst_rise_dur_vec,10);
    ub_dur = prctile(burst_rise_dur_vec,90);
    lb_sz = prctile(burst_rise_size_vec,10);
    ub_sz = prctile(burst_rise_size_vec,90);
    % make bin vectors
    dur_bins = [min(burst_rise_dur_vec) linspace(lb_dur,ub_dur,n_bins-1) max(burst_rise_dur_vec)+1];
    size_bins = [min(burst_rise_size_vec) linspace(lb_sz,ub_sz,n_bins-1) max(burst_rise_size_vec)+1];
    % get counts
    dur_ct_vec = histc(burst_rise_dur_vec,dur_bins);
    dur_ct_vec = dur_ct_vec / sum(dur_ct_vec);
    size_ct_vec = histc(burst_rise_size_vec,size_bins);
    size_ct_vec = size_ct_vec / sum(size_ct_vec);
    %%% initialize colormap 
    dur_ct_cm = brewermap(128,'BuGn');
    % generate burst duration series with constant(ish) burst amplitude
    burst_range = 3:12;
    burst_window = 3;
    burst_rise_dur_hmm_mean = NaN(numel(burst_range),window_size);
    burst_rise_dur_spot_mean = NaN(numel(burst_range),window_size);
    burst_rise_dur_swap_mean = NaN(numel(burst_range),window_size);
    rng(123)
    % iterate
    for i = 1:numel(burst_range)
        % filter for relevant subset 
        burst_vec = burst_range(i)-burst_window:burst_range(i)+burst_window;
        burst_ft = ismember(burst_rise_dur_vec,burst_vec);
        br_size_temp = burst_rise_size_vec(burst_ft);
        br_dur_temp = burst_rise_dur_vec(burst_ft);
        hmm_temp = hmm_array(burst_ft,:);
        spot_temp = spot_array(burst_ft,:);
        swap_temp = swap_array(burst_ft,:);
        % get distribtuion of burst amplitudes and calculate sample weights
        size_ct_vec_temp = histc(br_size_temp,size_bins);
        size_ct_vec_temp = size_ct_vec_temp / sum(size_ct_vec_temp);    
        size_weight_vec = size_ct_vec ./ size_ct_vec_temp;
        % assign weights
        sample_weight_vec =  zeros(size(br_dur_temp));
        for j = 1:numel(size_bins)-1
            sample_weight_vec(br_size_temp>=size_bins(j)&br_size_temp<size_bins(j+1)) = size_weight_vec(j);
        end
        index_vec = 1:numel(br_size_temp);
        sample_ids = randsample(index_vec,numel(index_vec),true,sample_weight_vec);
        % record sample averages
        burst_rise_dur_hmm_mean(i,:) = nanmean(hmm_temp(sample_ids,:));
        burst_rise_dur_spot_mean(i,:) = nanmean(spot_temp(sample_ids,:));
        burst_rise_dur_swap_mean(i,:) = nanmean(swap_temp(sample_ids,:));    
    end

    burst_rise_dur_hm = figure;
    colormap(hm_cm)
    % colormap(jet(128)/1.1);
    pcolor(flipud(burst_rise_dur_spot_mean))
    xlabel('offset (min)')
    set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)
    ylabel('burst duration (min)')
    set(gca,'yticklabels',fliplr(round(burst_range/3,1)))
    h = colorbar;
    ylabel(h, 'Dorsal levels (au)')
    set(gca,'FontSize', 12);
    saveas(burst_rise_dur_hm, [figPath 'burst_dur_rise_constant_r_hm_target.tif'])

    swap_rise_dur_hm = figure;
    colormap(hm_cm)
    % colormap(jet(128)/1.1);
    pcolor(flipud(burst_rise_dur_swap_mean))
    xlabel('offset (min)')
    set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)
    ylabel('burst duration (min)')
    set(gca,'yticklabels',fliplr(round(burst_range/3,1)))
    h = colorbar;
    ylabel(h, 'Dorsal levels (au)')
    set(gca,'FontSize', 12);
    saveas(swap_rise_dur_hm, [figPath 'burst_dur_rise_constant_r_hm_swap.tif'])

    hmm_rise_dur_hm = figure;
    colormap(hm_cm)
    pcolor(flipud(burst_rise_dur_hmm_mean))
    xlabel('offset (min)')
    set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)
    ylabel('burst duration (min)')
    set(gca,'yticklabels',fliplr(round(burst_range/3,1)))
    h = colorbar;
    ylabel(h, 'sna activity (au)')
    set(gca,'FontSize', 12);
    saveas(hmm_rise_dur_hm, [figPath 'burst_dur_rise_constant_r_hm_hmm.tif'])

    burst_rise_dur_wt = figure;
    index_vec = 1:numel(burst_range);
    hold on
    for i = 1:numel(burst_range)
        temp = burst_rise_dur_spot_mean;
        temp(index_vec~=i,:) = NaN;
        w = waterfall(temp-nanmin(temp(:)),repmat((1:numel(burst_range))',1,window_size));
        w.FaceColor = dur_ct_cm(1+(i-1)*inc,:);
        w.FaceAlpha = .8;
        w.EdgeColor = 'black';
    end
    xlabel('offset (min)')
    set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)    
    zlabel('Dorsal levels (au)')    
    view(-15,20)
    grid on
    saveas(burst_rise_dur_wt, [figPath 'burst_dur_rise_constant_r_waterfall_target.tif'])

    swap_rise_dur_wt = figure;
    index_vec = 1:numel(burst_range);
    hold on
    for i = 1:numel(burst_range)
        temp = burst_rise_dur_swap_mean;
        temp(index_vec~=i,:) = NaN;
        w = waterfall(temp-nanmin(temp(:)),repmat((1:numel(burst_range))',1,window_size));
        w.FaceColor = dur_ct_cm(1+(i-1)*inc,:);
        w.FaceAlpha = .8;
        w.EdgeColor = 'black';
    end
    xlabel('offset (min)')
    set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)    
    zlabel('Dorsal levels (au)')    
    view(-15,20)
    grid on
    saveas(swap_rise_dur_wt, [figPath 'burst_dur_rise_constant_r_waterfall_swap.tif'])

    hmm_rise_dur_wt = figure;
    index_vec = 1:numel(burst_range);
    hold on
    for i = 1:numel(burst_range)
        temp = burst_rise_dur_hmm_mean;
        temp(index_vec~=i,:) = NaN;
        w = waterfall(temp-nanmin(temp(:)),repmat((1:numel(burst_range))',1,window_size));
        w.FaceColor = dur_ct_cm(1+(i-1)*inc,:);
        w.FaceAlpha = .8;
        w.EdgeColor = 'black';
    end
    xlabel('offset (min)')
    set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)    
    zlabel('Dorsal levels (au)')    
    view(-15,20)
    grid on
    saveas(hmm_rise_dur_wt, [figPath 'burst_dur_rise_constant_r_waterfall_hmm.tif'])
    
    sz_ct_cm = brewermap(128,'PuRD');
    % generate burst duration series with constant(ish) burst amplitude
    size_range = .6:.1:1.5;
    sz_window = .3;
    burst_rise_size_hmm_mean = NaN(numel(size_range),window_size);
    burst_rise_size_spot_mean = NaN(numel(size_range),window_size);
    burst_rise_size_swap_mean = NaN(numel(size_range),window_size);
    rng(123)
    % iterate
    for i = 1:numel(size_range)
        % filter for relevant subset 
        lb_sz = size_range(i) - sz_window;
        ub_sz = size_range(i) + sz_window;
        burst_ft = burst_rise_size_vec>=lb_sz&burst_rise_size_vec<ub_sz;
        br_size_temp = burst_rise_size_vec(burst_ft);
        br_dur_temp = burst_rise_dur_vec(burst_ft);
        hmm_temp = hmm_array(burst_ft,:);
        spot_temp = spot_array(burst_ft,:);
        swap_temp = swap_array(burst_ft,:);
        % get distribtuion of burst amplitudes and calculate sample weights
        dur_ct_vec_temp = histc(br_dur_temp,dur_bins);
        dur_ct_vec_temp = dur_ct_vec_temp / sum(dur_ct_vec_temp);    
        dur_weight_vec = dur_ct_vec ./ dur_ct_vec_temp;
        % assign weights
        sample_weight_vec =  zeros(size(br_dur_temp));
        for j = 1:numel(dur_bins)-1
            sample_weight_vec(br_dur_temp>=dur_bins(j)&br_dur_temp<dur_bins(j+1)) = dur_weight_vec(j);
        end
        index_vec = 1:numel(br_size_temp);
        sample_ids = randsample(index_vec,numel(index_vec),true,sample_weight_vec);
        % record sample averages
        burst_rise_size_hmm_mean(i,:) = nanmean(hmm_temp(sample_ids,:));
        burst_rise_size_spot_mean(i,:) = nanmean(spot_temp(sample_ids,:));
        burst_rise_size_swap_mean(i,:) = nanmean(swap_temp(sample_ids,:));    
    end

    burst_rise_size_hm = figure;
    colormap(hm_cm)
    % colormap(jet(128)/1.1);
    pcolor(flipud(burst_rise_size_spot_mean))
    xlabel('offset (min)')
    set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)
    ylabel('burst amplitude (au)')
    set(gca,'yticklabels',fliplr(round(burst_range/3,1)))
    h = colorbar;
    ylabel(h, 'Dorsal levels (au)')
    set(gca,'FontSize', 12);
    saveas(burst_rise_size_hm, [figPath 'burst_size_rise_constant_dur_hm_target.tif'])

    swap_rise_size_hm = figure;
    colormap(hm_cm)
    % colormap(jet(128)/1.1);
    pcolor(flipud(burst_rise_size_swap_mean))
    xlabel('offset (min)')
    set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)
    ylabel('burst amplitude (au)')
    set(gca,'yticklabels',fliplr(round(burst_range/3,1)))
    h = colorbar;
    ylabel(h, 'Dorsal levels (au)')
    set(gca,'FontSize', 12);
    saveas(swap_rise_size_hm, [figPath 'burst_size_rise_constant_dur_hm_swap.tif'])

    hmm_rise_size_hm = figure;
    colormap(hm_cm)
    pcolor(flipud(burst_rise_size_hmm_mean))
    xlabel('offset (min)')
    set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)
    ylabel('burst amplitdue (au)')
    set(gca,'yticklabels',fliplr(round(burst_range/3,1)))
    h = colorbar;
    ylabel(h, 'sna activity (au)')
    set(gca,'FontSize', 12);
    saveas(hmm_rise_size_hm, [figPath 'burst_size_rise_constant_dur_hm_hmm.tif'])

    burst_rise_size_wt = figure;
    index_vec = 1:numel(burst_range);
    hold on
    for i = 1:numel(burst_range)
        temp = burst_rise_size_spot_mean;
        temp(index_vec~=i,:) = NaN;
        w = waterfall(temp-nanmin(temp(:)),repmat((1:numel(burst_range))',1,window_size));
        w.FaceColor = sz_ct_cm(1+(i-1)*inc,:);
        w.FaceAlpha = .8;
        w.EdgeColor = 'black';
    end
    xlabel('offset (min)')
    set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)    
    zlabel('Dorsal levels (au)')    
    view(-15,20)
    grid on
    saveas(burst_rise_size_wt, [figPath 'burst_size_rise_constant_dur_waterfall_target.tif'])

    swap_rise_size_wt = figure;
    index_vec = 1:numel(burst_range);
    hold on
    for i = 1:numel(burst_range)
        temp = burst_rise_size_swap_mean;
        temp(index_vec~=i,:) = NaN;
        w = waterfall(temp-nanmin(temp(:)),repmat((1:numel(burst_range))',1,window_size));
        w.FaceColor = sz_ct_cm(1+(i-1)*inc,:);
        w.FaceAlpha = .8;
        w.EdgeColor = 'black';
    end
    xlabel('offset (min)')
    set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)    
    zlabel('Dorsal levels (au)')    
    view(-15,20)
    grid on
    saveas(swap_rise_size_wt, [figPath 'burst_size_rise_constant_dur_waterfall_swap.tif'])

    hmm_rise_size_wt = figure;
    index_vec = 1:numel(burst_range);
    hold on
    for i = 1:numel(burst_range)
        temp = burst_rise_size_hmm_mean;
        temp(index_vec~=i,:) = NaN;
        w = waterfall(temp-nanmin(temp(:)),repmat((1:numel(burst_range))',1,window_size));
        w.FaceColor = sz_ct_cm(1+(i-1)*inc,:);
        w.FaceAlpha = .8;
        w.EdgeColor = 'black';
    end
    xlabel('offset (min)')
    set(gca,'xtick',1:3:window_size,'xticklabels',-5:5)    
    zlabel('Dorsal levels (au)')    
    view(-15,20)
    grid on
    saveas(hmm_rise_size_wt, [figPath 'burst_size_rise_constant_dur_waterfall_hmm.tif'])
end