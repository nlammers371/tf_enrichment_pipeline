clear
close all
clc

addpath(genpath('./lib'))

%% Initialization

projectName = 'optokni_eve4+6_WT'; 

liveProject = LiveEnrichmentProject(projectName);
resultsRoot = [liveProject.dataPath filesep];

% load data
load([resultsRoot 'spot_struct.mat'])
FigurePath = [liveProject.figurePath 'reactivation_dynamics_all' filesep];
mkdir(FigurePath)


embryo(1).expID = 1;
embryo(2).expID = 2;
embryo(3).expID = 3;
embryo(4).expID = 5;


% color to be used
k_green = brighten([38 142 75]/256,.4);
color_green = [38 143 75]/256; % color from Jake
mRNA_red = brighten([212 100 39]/256,.2);

eYFP_background = 375698.13; %prctile(double(knirps_vec_long),1);

ap_lim = 0.02; % AP range for analysis, 0.02 seems to be a reasonable number


% histogram parameters
binNum = 15;% best 14; ok, 15
binMax = 8;% best 8
%binNum = 10;%14
%binMax = 8;
edges = linspace(0,binMax,binNum);
%edges1 = linspace(0,10,11); % bin for histogram
%edges1 = linspace(0,15,16);
%edges = linspace(0,6,9); % bin for histogram

% temporary correction
%correction_factor = 1.346; % calibration factor
%correction_factor = 1.243;
%correction_factor = 1;

% timerange to analyze for response time
analysis_range = 8;


%% Figure: plot mean fluorescence vs time (not aligned)

time_full_long = [];
fluo_full_long = [];
knirps_full_long = [];
frame_full_long = [];

for i = 1:length(embryo)

    expID = embryo(i).expID;  
    nBoots = 100;

    ever_on_vec = [];
    mean_ap = [];
    time_orig_long = [];
    fluo_orig_long = [];
    frame_orig_long = [];
    off_time_long = [];
    knirps_orig_long = [];
    last_on_long = [];
    first_on_long = [];

    count = 0;

    for j = 1:length(spot_struct)

        if (spot_struct(j).TraceQCFlag == 1) && (spot_struct(j).setID == expID)
            % extract core vectors 

            % extract core vectors 
            fluo_vec_orig = spot_struct(j).fluo;
            time_vec_orig = spot_struct(j).time;
            frame_vec_orig = spot_struct(j).frames;
            knirps_vec_orig = spot_struct(j).rawNCProtein;
            ap_vec_orig = spot_struct(j).APPosNucleus;

            % calculate mean
            ever_on_orig = any(~isnan(fluo_vec_orig));
            mean_ap_orig = nanmean(ap_vec_orig);
            mean_knirps_orig = nanmean(knirps_vec_orig);


            if (mean_ap_orig > -ap_lim) && (mean_ap_orig < ap_lim)
               time_orig_long = [time_orig_long time_vec_orig];
               frame_orig_long = [frame_orig_long frame_vec_orig];
               fluo_orig_long = [fluo_orig_long fluo_vec_orig];
               knirps_orig_long = [knirps_orig_long knirps_vec_orig];


               %plot((time_vec_orig)/60,fluo_vec_orig,'Color', [175 175 175]/255);
               count = count + 1;
            end

    %         if (mean_ap(end) > -0.01) && (mean_ap(end) < 0.01)
    %             plot(time_vec-time_vec(end),fluo_vec);
    %             count = count + 1
    %         end

        end

    end
    
    
    % calculate mean knirps and fraction on (before/after perturbation)

    fluo_orig_long(isnan(fluo_orig_long)) = 0;

    frame_len = max(frame_orig_long);
    
    time_vec = zeros(1,frame_len);
    fluo_vec_mean = zeros(frame_len,1);
    fluo_vec_ste = zeros(frame_len,1);
    knirps_vec_mean = zeros(frame_len,1);
    knirps_vec_ste = zeros(frame_len,1);
    frac_on = zeros(frame_len,1);
    frac_on_ste = zeros(frame_len,1);

    fluo_orig_long_zero = fluo_orig_long;
    fluo_orig_long_zero(isnan(fluo_orig_long)) = 0;

    fluo_orig_long_binary = fluo_orig_long;
    fluo_orig_long_binary(isnan(fluo_orig_long)) = 0;
    fluo_orig_long_binary(fluo_orig_long>0) = 1;


    for j = 1:frame_len

        time_filter_long = frame_orig_long==j;

        %if sum(time_filter_long) > 10        
        %    boot_samples_fluo = bootstrp(nBoots,@nanmean,fluo_orig_long(time_filter_long));
        %    fluo_vec_mean(i) = nanmean(boot_samples_fluo);
        %    fluo_vec_ste(i) = std(boot_samples_fluo);

        time_vec(j) = mean(time_orig_long(time_filter_long))/60;

        fluo_vec_mean(j) = nanmean(fluo_orig_long_zero(time_filter_long));
        fluo_vec_ste(j) = std(fluo_orig_long_zero(time_filter_long),'omitnan');

        knirps_vec_mean(j) = nanmean(knirps_orig_long(time_filter_long));
        knirps_vec_ste(j) = std(knirps_orig_long(time_filter_long),'omitnan');

        frac_on(j) = nanmean(fluo_orig_long_binary(time_filter_long));
        frac_on_ste(j) = std(fluo_orig_long_binary(time_filter_long));

    end

    time_vec_on = time_vec;

    %knirps_vec_mean(time_vec_on>=0) = knirps_vec_mean(time_vec_on>=0)/correction_factor;
    %knirps_vec_mean(time_vec_on>=0) = convert_from_458(knirps_vec_mean(time_vec_on>=0));    
    knirps_vec_mean = knirps_vec_mean-eYFP_background;
    
    
    % record results for combining the traces
    time_vec_temp = time_orig_long/60;
    time_full_long = [time_full_long time_vec_temp];
    fluo_full_long = [fluo_full_long fluo_orig_long];
    knirps_full_long = [knirps_full_long knirps_orig_long];
    frame_full_long = [frame_full_long frame_orig_long];
    
    
    temp_traj_fig  = figure('Position',[10 10 800 800]);
    tiledlayout(2,1)
    nexttile
    hold on
    %time_interp = min(time_vec_on):0.1:max(time_vec_on);
    time_interp = -10:0.1:10;
    frac_on_mean = movmean(frac_on,5);
    %frac_on_interp = interp1(time_vec_on(~isnan(frac_on_mean)),frac_on_mean(~isnan(frac_on_mean)),time_interp,'spline');
    frac_on_interp = movmean(frac_on_interp,5);
    %frac_on_interp = interp1(time_vec_on,frac_on,time_interp,'v5cubic');

    %errorbar(time_vec_on,knirps_vec_ste,'Color','k','CapSize',0);
    plot(time_vec_on,knirps_vec_mean,'-k','LineWidth',1)
    scatter(time_vec_on,knirps_vec_mean,50,'MarkerFaceColor',k_green,'MarkerEdgeColor','k')
    %xlim([-10 7])
    %ylim([3.75E5 9E5])
    xlabel(['time relative to perturbation (min)'])
    ylabel(['Knirps concentration (AU)'])
    pbaspect([3 2 1])
    
    nexttile
    hold on
    plot(time_vec_on,frac_on,'.')
    %plot(time_interp,frac_on_interp,'-','LineWidth',2);
    scatter(time_vec_on,frac_on,50,'MarkerFaceColor',mRNA_red,'MarkerEdgeColor','k')
    %xlim([-10 7])
    %ylim([0 1])
    xlabel(['time relative to perturbation (min)'])
    ylabel(['fraction of nuclei on'])
    pbaspect([3 2 1])
    
    %saveas(temp_traj_fig,[FigurePath 'figure_temporal_trajectory_' num2str(i) '.pdf'])

end

%% plot combined results

% calculate mean knirps and fraction on (before/after perturbation)

fluo_full_long(isnan(fluo_full_long)) = 0;

%time_bin_full = double((0:max(frame_full_long)))+0.5;
time_bin_full = linspace(-15,15,58);% 58,56, best
time_vec_plot = (time_bin_full(1:end-1) + time_bin_full(2:end))/2;
time_groups_full = discretize(time_aligned_full_long,time_bin_full);
%time_vec_full = zeros(1,length(time_bin_full)-1);

fluo_vec_full_mean = zeros(length(time_bin_full)-1,1);
fluo_vec_full_ste = zeros(length(time_bin_full)-1,1);

knirps_vec_full_mean = zeros(length(time_bin_full)-1,1);
knirps_vec_full_ste = zeros(length(time_bin_full)-1,1);

frac_on_full = zeros(length(time_bin_full)-1,1);
frac_on_full_ste = zeros(length(time_bin_full)-1,1);

fluo_full_long_zero = fluo_full_long;
fluo_full_long_zero(isnan(fluo_full_long)) = 0;

fluo_full_long_binary = fluo_full_long;
fluo_full_long_binary(isnan(fluo_full_long)) = 0;
fluo_full_long_binary(fluo_full_long>0) = 1;


for j = 1:length(time_bin_full)-1

    time_filter_long_full = time_groups_full==j;

    time_vec_full(j) = mean(time_aligned_full_long(time_filter_long_full));

    fluo_vec_full_mean(j) = nanmean(fluo_full_long_zero(time_filter_long_full));
    fluo_vec_full_ste(j) = std(fluo_full_long_zero(time_filter_long_full),'omitnan');

    knirps_vec_full_mean(j) = nanmean(knirps_full_long(time_filter_long_full));
    knirps_vec_full_ste(j) = std(knirps_full_long(time_filter_long_full),'omitnan');

    frac_on_full(j) = nanmean(fluo_full_long_binary(time_filter_long_full));
    frac_on_full_ste(j) = std(fluo_full_long_binary(time_filter_long_full));

end

%knirps_vec_full_mean(time_vec_plot>=0) = knirps_vec_full_mean(time_vec_plot>=0)/correction_factor;
knirps_vec_full_mean(time_vec_plot>=0) = convert_from_458(knirps_vec_full_mean(time_vec_plot>=0));
knirps_vec_full_mean = knirps_vec_full_mean - eYFP_background;

combined_traj_fig  = figure('Position',[10 10 800 800]);
tiledlayout(2,1)
nexttile

%time_full_interp = -10:0.1:10;
%knirps_vec_full_movmean = movmean(knirps_vec_full_mean,3);
%knirps_vec_full_interp = interp1(time_vec_full(~isnan(knirps_vec_full_movmean)),knirps_vec_full_movmean(~isnan(knirps_vec_full_movmean)),time_full_interp,'spline');
%knirps_vec_full_interp = movmean(knirps_vec_full_interp,1);

hold on
%plot(time_full_interp,knirps_vec_full_interp,'-k','LineWidth',2)
plot(time_vec_plot,knirps_vec_full_mean,'-k','LineWidth',1)
scatter(time_vec_plot,knirps_vec_full_mean,50,'MarkerFaceColor',k_green,'MarkerEdgeColor','k')
%plot(time_vec_plot,knirps_vec_full_mean,'o');
xlim([-10 7])
ylim([3.75E5 10E5])
xlabel(['time relative to perturbation (min)'])
ylabel(['Knirps concentration (AU)'])
pbaspect([3 2 1])

nexttile

time_full_interp = -10:0.1:10;
frac_on_full_mean = movmean(frac_on_full,3);
frac_on_full_interp = interp1(time_vec_full(~isnan(frac_on_full_mean)),frac_on_full_mean(~isnan(frac_on_full_mean)),time_full_interp,'spline');
frac_on_full_interp = movmean(frac_on_full_interp,1);

hold on
plot(time_vec_plot,frac_on_full,'o')
plot(time_full_interp,frac_on_full_interp,'-','LineWidth',2);
scatter(time_vec_plot,frac_on_full,50,'MarkerFaceColor',mRNA_red,'MarkerEdgeColor','k')
xlabel(['time relative to perturbation (min)'])
ylabel(['fraction of nuclei on'])
xlim([-10 7])
ylim([0 0.9])
pbaspect([3 2 1])