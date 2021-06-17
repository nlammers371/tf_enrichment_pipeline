clc
clear
close all

addpath(genpath('./lib'))

%% create color map
%Green color map v1
cmap_green = uint8([[247,252,245];[229,245,224];[199,233,192];[161,217,155];[116,196,118];[65,171,93];[35,139,69];[0,109,44];[0,68,27]]);

old_stepNum = size(cmap_green,1);
new_stepNum = 256;

cmap_green_1 = interp1(linspace(0,1,old_stepNum),double(cmap_green(:,1)),linspace(0,1,new_stepNum));
cmap_green_2 = interp1(linspace(0,1,old_stepNum),double(cmap_green(:,2)),linspace(0,1,new_stepNum));
cmap_green_3 = interp1(linspace(0,1,old_stepNum),double(cmap_green(:,3)),linspace(0,1,new_stepNum));

cmap_green_new = [cmap_green_1' cmap_green_2' cmap_green_3']/256;

%Green color map v2
cmap_green_v2 = [[247,252,253];[229,245,249];[204,236,230];[153,216,201];[102,194,164];[65,174,118];[35,139,69];[0,109,44];[0,68,27]];

old_stepNum = size(cmap_green_v2,1);
new_stepNum = 256;

cmap_green_1_v2 = interp1(linspace(0,1,old_stepNum),double(cmap_green_v2(:,1)),linspace(0,1,new_stepNum));
cmap_green_2_v2 = interp1(linspace(0,1,old_stepNum),double(cmap_green_v2(:,2)),linspace(0,1,new_stepNum));
cmap_green_3_v2 = interp1(linspace(0,1,old_stepNum),double(cmap_green_v2(:,3)),linspace(0,1,new_stepNum));

cmap_green_new_v2 = [cmap_green_1_v2' cmap_green_2_v2' cmap_green_3_v2']/256;

%Red color map v1
cmap_red = [[255,255,229];[255,247,188];[254,227,145];[254,196,79];[254,153,41];[236,112,20];[204,76,2];[153,52,4];[102,37,6]];

old_stepNum = size(cmap_red,1);
new_stepNum = 256;

cmap_red_1 = interp1(linspace(0,1,old_stepNum),double(cmap_red(:,1)),linspace(0,1,new_stepNum));
cmap_red_2 = interp1(linspace(0,1,old_stepNum),double(cmap_red(:,2)),linspace(0,1,new_stepNum));
cmap_red_3 = interp1(linspace(0,1,old_stepNum),double(cmap_red(:,3)),linspace(0,1,new_stepNum));

cmap_red_new = [cmap_red_1' cmap_red_2' cmap_red_3']/256;

%Color for lines
color_green = [38 143 75]/256;
color_red = [209 82 5]/256;

left_color = color_green;
right_color = color_red;


% Define some other colors  
yw = [234 194 100]/255; % yellow
bl = [115 143 193]/255; % blue
gr = [191 213 151]/255; % green

%% initialization

projectName = {'optokni_eve4+6_WT','optokni_eve4+6_ON_CONST'}; 

for i = 1:length(projectName)

liveProject = LiveEnrichmentProject(projectName{i});
resultsRoot = [liveProject.dataPath filesep];    

% load data
load([resultsRoot 'spot_struct.mat'])
FigurePath = [liveProject.figurePath 'WT_vs_CONST' filesep];
mkdir(FigurePath)

% MS2 trace quality looks atrocious. Let's look at on and off times and the
% like

ever_on_vec = zeros(size(spot_struct));
off_time_vec = NaN(size(spot_struct));
on_time_vec = NaN(size(spot_struct));
off_spot_fluo = NaN(size(spot_struct));
off_knirps_vec = NaN(size(spot_struct));
off_ap = NaN(size(spot_struct));
on_knirps_vec = NaN(size(spot_struct));
on_ap = NaN(size(spot_struct));
mean_ap = NaN(size(spot_struct));
mean_knirps = NaN(size(spot_struct));
mean_fluo = NaN(size(spot_struct));
% initialize longform vectors for regression
ap_vec_long = [];
time_vec_long = [];
knirps_vec_long_raw = [];
fluo_raw_long = [];
fluo_zeros_long = [];
mRNA_vec_long = [];

post_turn_on_flags = [];
post_turn_off_flags = [];
ever_on_flags = [];

for j = 1:length(spot_struct)
  
    % extract core vectors 
    fluo_vec = spot_struct(j).fluo;
    time_vec = spot_struct(j).time;
    knirps_vec = spot_struct(j).rawNCProtein;
    ap_vec = spot_struct(j).apPosNucleus;
    
    if time_vec(end) - time_vec(1) >=20*60
        
        % make average vectors        
      
        % get off and on indices
        ever_on_vec(j) = any(~isnan(fluo_vec));
        mean_ap(j) = nanmean(ap_vec);
        mean_knirps(j) = nanmean(knirps_vec);
        
        post_on_vec = zeros(size(ap_vec));
        post_off_vec = zeros(size(ap_vec));
        if ever_on_vec(j)
            start_i = find(~isnan(fluo_vec),1);
            stop_i = find(~isnan(fluo_vec),1,'last');
            if true%stop_i < length(fluo_vec)-10
                off_time_vec(j) = time_vec(stop_i);
                off_knirps_vec(j) = knirps_vec(stop_i);
                off_spot_fluo(j) = fluo_vec(stop_i);
                off_ap(j) = ap_vec(stop_i);
                post_off_vec(stop_i+1:end) = 1;
            end
            
            if start_i > 1
                on_time_vec(j) = time_vec(start_i);
                on_knirps_vec(j) = knirps_vec(start_i);
                on_ap(j) = ap_vec(start_i);
                post_on_vec(start_i+1:end) = 1;
            end
        end
        
        % make regression vectors
        ever_on_flags = [ever_on_flags repelem(ever_on_vec(j),length(ap_vec))];
        post_turn_on_flags = [post_turn_on_flags post_on_vec];
        post_turn_off_flags = [post_turn_off_flags post_off_vec];
        fluo_raw_long = [fluo_raw_long fluo_vec];
        
        fluo_zeros = fluo_vec;
        all_zeros = fluo_vec;
        fluo_zeros(post_on_vec&~post_off_vec&isnan(fluo_vec)) = 0;
        fluo_zeros_long = [fluo_zeros_long fluo_zeros];
        
        all_zeros(isnan(all_zeros)) = 0;
        mRNA_vec_long = [mRNA_vec_long all_zeros];
        
        mean_fluo(j) = nanmean(fluo_zeros);
        
        knirps_vec_long_raw = [knirps_vec_long_raw knirps_vec];
        ap_vec_long = [ap_vec_long ap_vec];
        time_vec_long = [time_vec_long time_vec];
        
    end
end



% Look at long vectors

%set parameters
%nBoots = 1000;
timeBins = 61;

ap_bins = linspace(-0.2,0.2,26);
ap_bins_plot = (ap_bins(1:end-1) + ap_bins(2:end))/2;
time_bins = linspace(0,35*60,timeBins);
%time_bins_plot = (time_bins(1:end-1) + time_bins(2:end))/2/60;
time_bins_plot = time_bins(1:end-1)/60;
%knirps_bins = linspace(0,15e5,nLinBins);

knirps_offset = 2.5e5;%prctile(double(knirps_vec_long),1);

% calculate mean vectors
knirps_vec_long = knirps_vec_long_raw - knirps_offset;

ap_groups = discretize(ap_vec_long,ap_bins); 
time_groups = discretize(time_vec_long,time_bins); 
%knirps_groups = discretize(knirps_vec_long,knirps_bins); 

frac_inst_on_time_array_mean = NaN(length(time_bins)-1,length(ap_bins)-1);
frac_on_time_array_mean = NaN(length(time_bins)-1,length(ap_bins)-1);
%frac_on_knirps_array_mean = NaN(length(time_bins)-1,length(ap_bins)-1);
eve_time_array_mean = NaN(length(time_bins)-1,length(ap_bins)-1);
%eve_time_array_full = NaN(length(time_bins)-1,length(ap_bins)-1,nBoots);
knirps_time_array_mean = NaN(length(time_bins)-1,length(ap_bins)-1);

frac_inst_on_time_array_std = NaN(length(time_bins)-1,length(ap_bins)-1);
frac_on_time_array_std = NaN(length(time_bins)-1,length(ap_bins)-1);
%frac_on_knirps_array_std = NaN(length(time_bins)-1,length(ap_bins)-1);
eve_time_array_std = NaN(length(time_bins)-1,length(ap_bins)-1);
knirps_time_array_std = NaN(length(time_bins)-1,length(ap_bins)-1);

frac_on_time_array_num = NaN(length(time_bins)-1,length(ap_bins)-1);
frac_inst_on_time_array_num = NaN(length(time_bins)-1,length(ap_bins)-1);
%frac_on_knirps_array_num = NaN(length(time_bins)-1,length(ap_bins)-1);
eve_time_array_num = NaN(length(time_bins)-1,length(ap_bins)-1);
knirps_time_array_num = NaN(length(time_bins)-1,length(ap_bins)-1);


frac_on_time_array_ste = NaN(length(time_bins)-1,length(ap_bins)-1);
frac_on_knirps_array_ste = NaN(length(time_bins)-1,length(ap_bins)-1);
eve_time_array_ste = NaN(length(time_bins)-1,length(ap_bins)-1);
knirps_time_array_ste = NaN(length(time_bins)-1,length(ap_bins)-1);


for t = 1:length(time_bins)
    t
    for a = 1:length(ap_bins)
        time_window_filter = time_groups==t & ap_groups==a & post_turn_on_flags;
        if sum(time_window_filter) > 10
            still_on_flags_time = ~post_turn_off_flags(time_window_filter);
            inst_on_flags_time = mRNA_vec_long(time_window_filter)>0;
            %boot_samples_time = bootstrp(nBoots,@mean,still_on_flags_time);
            %frac_on_time_array_mean(t,a) = mean(boot_samples_time);
            %frac_on_time_array_ste(t,a) = std(boot_samples_time);
            frac_on_time_array_mean(t,a) = nanmean(still_on_flags_time);
            frac_on_time_array_num(t,a) = sum(~isnan(still_on_flags_time));
            frac_inst_on_time_array_mean(t,a) = mean(inst_on_flags_time);
            frac_inst_on_time_array_num(t,a) = length(inst_on_flags_time);
            
            %frac_on_time_array_std(t,a) = std(boot_samples_time);
            %frac_on_time_array_ste(t,a) = frac_on_time_array_std(t,a)/sqrt(frac_on_time_array_num(t,a));
            
            %boot_samples_mRNA = bootstrp(nBoots,@nanmean,mRNA_vec_long(time_window_filter));
            %eve_time_array_mean(t,a) = mean(boot_samples_mRNA);
            %eve_time_array_ste(t,a) = std(boot_samples_mRNA);
            %eve_time_array_full(t,a,:) = boot_samples_mRNA; % NL: note that this is not quite right. Should really be sampling at the level of loci 
            eve_sample = mRNA_vec_long(time_window_filter);
            eve_sample(isnan(eve_sample)) = 0;
            eve_time_array_mean(t,a) = mean(eve_sample);
            eve_time_array_num(t,a) = sum(~isnan(eve_sample));
            eve_time_array_std(t,a) = std(eve_sample);
            eve_time_array_ste(t,a) = eve_time_array_std(t,a)/sqrt(eve_time_array_num(t,a));
            
            %boot_samples_knirps = bootstrp(nBoots,@nanmean,knirps_vec_long_raw(time_window_filter));
            %knirps_time_array_mean(t,a) = nanmean(boot_samples_knirps); % NL: note that this is not quite right. Should really be sampling at the level of loci 
            %knirps_time_array_ste(t,a) = nanstd(boot_samples_knirps);
            knirps_sample = knirps_vec_long(time_window_filter);
            knirps_sample(isnan(knirps_sample)) = 0;
            knirps_time_array_mean(t,a) = mean(knirps_sample);
            knirps_time_array_num(t,a) = sum(~isnan(knirps_sample));
            knirps_time_array_std(t,a) = std(knirps_sample);
            knirps_time_array_ste(t,a) = knirps_time_array_std(t,a)/sqrt(knirps_time_array_num(t,a));
            
        end
                
%         knirps_window_filter = knirps_groups==t & ap_groups==a & post_turn_on_flags;
%         if sum(knirps_window_filter)>10
%             still_on_flags_knirps = ~post_turn_off_flags(knirps_window_filter);
%             boot_samples_knirps = bootstrp(nBoots,@mean,still_on_flags_knirps);
%             frac_on_knirps_array_mean(t,a) = mean(boot_samples_knirps);
%             frac_on_knirps_array_ste(t,a) = std(boot_samples_knirps);
%         end
    end
end
knirps_time_array_mean(isnan(knirps_time_array_mean)) = 0;
eve_time_array_mean(isnan(eve_time_array_mean)) = 0;
frac_inst_on_time_array_mean(isnan(frac_inst_on_time_array_mean)) = 0;

knirps_time_array_mean = movmean(knirps_time_array_mean,5,1);
eve_time_array_mean = movmean(eve_time_array_mean,5,1);
frac_inst_on_time_array_mean = movmean(frac_inst_on_time_array_mean,5,1);

if i == 1
    WT.knirps_mean = knirps_time_array_mean;
    WT.eve_mean = eve_time_array_mean;
    WT.inst_on = frac_inst_on_time_array_mean;
else
    CONST.knirps_mean = knirps_time_array_mean;
    CONST.eve_mean = eve_time_array_mean;
    CONST.inst_on = frac_inst_on_time_array_mean;
end

end

%% plot calculated mean vector

% plot for wildtype
WT_mean_kni_fig = figure;
imagesc(ap_bins_plot,time_bins_plot,WT.knirps_mean)
colorbar
colormap(cmap_green_new_v2)
caxis([0 10E5])
pbaspect([3 2 1])
ylim([0 35])
xlabel('AP position (% embryo length)')
ylabel('time (min)')


WT_mean_eve_fig = figure;
imagesc(ap_bins_plot,time_bins_plot,WT.eve_mean)
colorbar
colormap(cmap_red_new)
caxis([0 2.1E5])
pbaspect([3 2 1])
ylim([0 35])
xlabel('AP position (% embryo length)')
ylabel('time (min)')

% % plot for const export
% CONST_mean_kni_fig = figure;
% imagesc(ap_bins_plot,time_bins_plot,CONST.knirps_mean)
% colorbar
% colormap(cmap_green_new_v2)
% caxis([0 15E5])
% pbaspect([3 2 1])
% ylim([0 35])
% xlabel('AP position (% embryo length)')
% ylabel('time (min)')


CONST_mean_eve_fig = figure;
imagesc(ap_bins_plot,time_bins_plot,CONST.eve_mean)
colorbar
colormap(cmap_red_new)
caxis([0 3.15E5])
pbaspect([3 2 1])
ylim([0 35])
xlabel('AP position (% embryo length)')
ylabel('time (min)')

WT_frac_eve_fig = figure;
imagesc(ap_bins_plot,time_bins_plot,WT.inst_on)
colorbar
colormap(cmap_red_new)
%colormap(jet)
caxis([0 1])
pbaspect([3 2 1])
ylim([0 35])
xlabel('AP position (% embryo length)')
ylabel('time (min)')

CONST_frac_eve_fig = figure;
imagesc(ap_bins_plot,time_bins_plot,CONST.inst_on)
colorbar
colormap(cmap_red_new)
%colormap(jet)
caxis([0 1])
pbaspect([3 2 1])
ylim([0 35])
xlabel('AP position (% embryo length)')
ylabel('time (min)')


%% Figure: Plot k_on, k_off vs ap

% load inference results
WT_bursting = load('P:\Jake\Pipeline\tf_enrichment_pipeline\opto_knirps_project\bursting_data\compiledResults_w7_K2_p0_ap7_t1_f2D_CONST.mat');
CONST_bursting = load('P:\Jake\Pipeline\tf_enrichment_pipeline\opto_knirps_project\bursting_data\compiledResults_w7_K2_p0_ap7_t1_f2D_WT_v1.mat');

burst_axis = (WT_bursting.compiledResults.apBins(1:end-1)+ WT_bursting.compiledResults.apBins(2:end))/2;
% data from WT
burst_freq_WT = WT_bursting.compiledResults.freq_vec_mean;
burst_freq_ste_WT = WT_bursting.compiledResults.freq_vec_ste;
burst_dur_WT = WT_bursting.compiledResults.dur_vec_mean;
burst_dur_ste_WT = WT_bursting.compiledResults.dur_vec_ste;
burst_rate_WT= WT_bursting.compiledResults.init_vec_mean*1e-5;
burst_rate_ste_WT = WT_bursting.compiledResults.init_vec_ste*1e-5;

% data from CONST
burst_freq_CONST = CONST_bursting.compiledResults.freq_vec_mean;
burst_freq_ste_CONST = CONST_bursting.compiledResults.freq_vec_ste;
burst_dur_CONST = CONST_bursting.compiledResults.dur_vec_mean;
burst_dur_ste_CONST = CONST_bursting.compiledResults.dur_vec_ste;
burst_rate_CONST= CONST_bursting.compiledResults.init_vec_mean*1e-5;
burst_rate_ste_CONST = CONST_bursting.compiledResults.init_vec_ste*1e-5;

burst_dur_center = burst_dur_WT(4);
burst_rate_center = burst_rate_WT(4);
burst_freq_center = burst_freq_WT(4);


ap_bins = linspace(-0.12,0.12,21);
ap_axis = 100*(ap_bins(1:end-1) + diff(ap_bins)/2);

burst_k_off_center = 1/burst_dur_center;
burst_k_on_center = 1/burst_freq_center;

burst_k_off_mean_WT = 1./WT_bursting.compiledResults.dur_vec_mean;
burst_k_off_ste_WT = 1./(WT_bursting.compiledResults.dur_vec_mean.^2).*WT_bursting.compiledResults.dur_vec_ste;
burst_k_on_mean_WT = 1./WT_bursting.compiledResults.freq_vec_mean;
burst_k_on_ste_WT = 1./(WT_bursting.compiledResults.freq_vec_mean.^2).*WT_bursting.compiledResults.freq_vec_ste;

burst_k_off_mean_CONST = 1./CONST_bursting.compiledResults.dur_vec_mean;
burst_k_off_ste_CONST = 1./(CONST_bursting.compiledResults.dur_vec_mean.^2).*CONST_bursting.compiledResults.dur_vec_ste;
burst_k_on_mean_CONST = 1./CONST_bursting.compiledResults.freq_vec_mean;
burst_k_on_ste_CONST = 1./(CONST_bursting.compiledResults.freq_vec_mean.^2).*CONST_bursting.compiledResults.freq_vec_ste;


% k_off vs ap
burst_k_off_fig = figure;
hold on

set(gca,'FontSize',14)

%WT
errorbar(burst_axis*100,burst_k_off_mean_WT/burst_k_off_center,burst_k_off_ste_WT/burst_k_off_center,'Color','k','CapSize',0)
plot(burst_axis*100,burst_k_off_mean_WT/burst_k_off_center,'-k')
scatter(burst_axis*100,burst_k_off_mean_WT/burst_k_off_center,50,'MarkerFaceColor',bl,'MarkerEdgeColor','k')
%CONST
errorbar(burst_axis*100,burst_k_off_mean_CONST/burst_k_off_center,burst_k_off_ste_CONST/burst_k_off_center,'Color','k','CapSize',0)
plot(burst_axis*100,burst_k_off_mean_CONST/burst_k_off_center,'-k')
scatter(burst_axis*100,burst_k_off_mean_CONST/burst_k_off_center,50,'MarkerFaceColor',bl,'MarkerEdgeColor','k')
%set(gca,'YColor',bl);
ylabel(['k_{off} (relative to center)'])
ylim([0 2])


xlabel('AP position (% embryo length)');
xlim([ap_axis(1) ap_axis(end)])
%mRNA_fig2.InvertHardcopy = 'off';
set(gcf,'color','w'); 
pbaspect([2 1 1])


% k_on vs ap
burst_k_on_fig = figure;
hold on

set(gca,'FontSize',14)
%WT
errorbar(burst_axis*100,burst_k_on_mean_WT/burst_k_on_center,burst_k_on_ste_WT/burst_k_on_center,'Color','k','CapSize',0)
plot(burst_axis*100,burst_k_on_mean_WT/burst_k_on_center,'-k')
scatter(burst_axis*100,burst_k_on_mean_WT/burst_k_on_center,50,'MarkerFaceColor',bl,'MarkerEdgeColor','k')
%CONST
errorbar(burst_axis*100,burst_k_on_mean_CONST/burst_k_on_center,burst_k_on_ste_CONST/burst_k_on_center,'Color','k','CapSize',0)
plot(burst_axis*100,burst_k_on_mean_CONST/burst_k_on_center,'-k')
scatter(burst_axis*100,burst_k_on_mean_CONST/burst_k_on_center,50,'MarkerFaceColor',bl,'MarkerEdgeColor','k')


%set(gca,'YColor',bl);
ylabel(['k_{on} (relative to center)'])
ylim([0 2])


xlabel('AP position (% embryo length)');
xlim([ap_axis(1) ap_axis(end)])
%mRNA_fig2.InvertHardcopy = 'off';
set(gcf,'color','w'); 
pbaspect([2 1 1])


% loading rate vs ap
burst_loading_rate_ap_fig = figure;
hold on

set(gca,'FontSize',14)

errorbar(burst_axis*100,burst_rate_WT/burst_rate_center,burst_rate_ste_WT/burst_rate_center,'Color','k','CapSize',0)
plot(burst_axis*100,burst_rate_WT/burst_rate_center,'-k')
scatter(burst_axis*100,burst_rate_WT/burst_rate_center,50,'MarkerFaceColor',gr,'MarkerEdgeColor','k')
%set(gca,'YColor',gr);
ylabel(['mRNA loading rate (relative to center)'])
ylim([0 2])


xlabel('AP position (% embryo length)');
xlim([ap_axis(1) ap_axis(end)])
%mRNA_fig2.InvertHardcopy = 'off';
set(gcf,'color','w'); 
pbaspect([2 1 1])


%saveas(burst_k_on_fig,[FigurePath 'figure_k_on_vs_ap.pdf'])
%saveas(burst_k_off_fig,[FigurePath 'figure_k_off_vs_ap.pdf'])
%saveas(burst_loading_rate_ap_fig,[FigurePath 'figure_loading_rate_vs_ap.pdf'])

%% Bursting parameters with mRNA pattern
% which time to plot for knirps
% time_plot_2 = 27;
% 
% % mRNA profile vs burst duration
% burst_dur_fig = figure;
% hold on
% 
% yyaxis left
% 
% eve_mRNA_sm = imgaussfilt(predicted_eve_profile_mean(time_plot_2,:)*1e-5,1);
% f = fill([ap_axis fliplr(ap_axis)], [eve_mRNA_sm zeros(size(eve_mRNA_sm))],mRNA_red);
% f.FaceAlpha = 0.3;
% set(gca,'YColor',mRNA_red)
% ylabel('accumulated {\it eve} mRNA (au)');
% ylim([0 8])
% 
% set(gca,'FontSize',14)
% 
% yyaxis right
% errorbar(burst_axis*100,burst_dur/burst_dur_center,burst_dur_ste/burst_dur_center,'Color','k','CapSize',0)
% plot(burst_axis*100,burst_dur/burst_dur_center,'-k')
% scatter(burst_axis*100,burst_dur/burst_dur_center,50,'MarkerFaceColor',bl,'MarkerEdgeColor','k')
% set(gca,'YColor',bl);
% ylabel(['burst duration (relative to center)'])
% ylim([0 2])
% 
% 
% xlabel('AP position (% embryo length)');
% xlim([ap_axis(1) ap_axis(end)])
% %mRNA_fig2.InvertHardcopy = 'off';
% set(gcf,'color','w'); 
% pbaspect([2 1 1])
% 
% 
% % mRNA profile vs loading rate
% burst_loading_rate_fig = figure;
% hold on
% 
% yyaxis left
% 
% eve_mRNA_sm = imgaussfilt(predicted_eve_profile_mean(time_plot_2,:)*1e-5,1);
% f = fill([ap_axis fliplr(ap_axis)], [eve_mRNA_sm zeros(size(eve_mRNA_sm))],mRNA_red);
% f.FaceAlpha = 0.3;
% set(gca,'YColor',mRNA_red)
% ylabel('accumulated {\it eve} mRNA (au)');
% ylim([0 8])
% 
% set(gca,'FontSize',14)
% 
% yyaxis right
% errorbar(burst_axis*100,burst_rate/burst_rate_center,burst_rate_ste/burst_rate_center,'Color','k','CapSize',0)
% plot(burst_axis*100,burst_rate/burst_rate_center,'-k')
% scatter(burst_axis*100,burst_rate/burst_rate_center,50,'MarkerFaceColor',gr,'MarkerEdgeColor','k')
% set(gca,'YColor',gr);
% ylabel(['mRNA loading rate (relative to center)'])
% ylim([0 2])
% 
% 
% xlabel('AP position (% embryo length)');
% xlim([ap_axis(1) ap_axis(end)])
% %mRNA_fig2.InvertHardcopy = 'off';
% set(gcf,'color','w'); 
% pbaspect([2 1 1])
% 
% % mRNA profile vs burst frequency
% burst_freq_fig = figure;
% hold on
% 
% yyaxis left
% 
% eve_mRNA_sm = imgaussfilt(predicted_eve_profile_mean(time_plot_2,:)*1e-5,1);
% f = fill([ap_axis fliplr(ap_axis)], [eve_mRNA_sm zeros(size(eve_mRNA_sm))],mRNA_red);
% f.FaceAlpha = 0.3;
% set(gca,'YColor',mRNA_red)
% ylabel('accumulated {\it eve} mRNA (au)');
% ylim([0 8])
% 
% set(gca,'FontSize',14)
% 
% yyaxis right
% errorbar(burst_axis*100,burst_freq,burst_freq_ste,'Color','k','CapSize',0)
% plot(burst_axis*100,burst_freq,'-k')
% scatter(burst_axis*100,burst_freq,50,'MarkerFaceColor',bl,'MarkerEdgeColor','k')
% set(gca,'YColor',bl);
% ylabel(['burst frequency (1/min)'])
% ylim([0 4])
% 
% 
% xlabel('AP position (% embryo length)');
% xlim([ap_axis(1) ap_axis(end)])
% %mRNA_fig2.InvertHardcopy = 'off';
% set(gcf,'color','w'); 
% pbaspect([2 1 1])
% 
% %saveas(burst_loading_rate_fig,[FigurePath 'figure_loading_rate_vs_ap_with_mRNA.pdf'])
% %saveas(burst_freq_fig,[FigurePath 'figure_burst_freq_vs_ap_with_mRNA.pdf'])
% %saveas(burst_dur_fig,[FigurePath 'figure_burst_dur_vs_ap_with_mRNA.pdf'])


% %% Figure 2: compare predicted mRNA
% % generate decay kernel
% close all
% 
% % which time to plot for knirps
% time_plot_1 = 5;
% time_plot_2 = 27;
% 
% eve_half_life = 7;
% eve_decay_kernel = exp(-time_bins'/eve_half_life/60);
% eve_decay_kernel = eve_decay_kernel / eve_decay_kernel(1);
% 
% % replace missing values with zeros (Need to clean this up, should really
% % be some kind of interpolation)
% eve_time_array_full(isnan(eve_time_array_full)) = 0;
% 
% predicted_eve_profile_array = convn(eve_decay_kernel,eve_time_array_full,'full');
% predicted_eve_profile_mean = nanmean(predicted_eve_profile_array,3);
% predicted_eve_profile_ste = nanstd(predicted_eve_profile_array,[],3);
% 
% predicted_eve_profile_mean = predicted_eve_profile_mean(1:length(time_bins),:);
% predicted_eve_profile_ste = predicted_eve_profile_ste(1:length(time_bins),:);
% 
% % knirps green
% k_green = brighten([38 142 75]/256,.4);
% color_green = [38 143 75]/256; % color from Jake
% mRNA_red = brighten([212 100 39]/256,.2);
% 
% % mRNA profile plot
% mRNA_fig1 = figure;
% hold on
% 
% yyaxis left
% f = fill([ap_axis fliplr(ap_axis)], [knirps_time_array_mean(time_plot_1,:)*1e-5 zeros(size(knirps_time_array_mean(time_plot_1,:)))],color_green);
% %plot(ap_axis,knirps_time_array_mean(time_plot_1,:)*1e-5,'Color',color_green,'LineWidth',3)
% f.FaceAlpha = 0.5;
% ylabel('[Knirps] (au)');
% set(gca,'YColor',color_green)
% ylim([0 15])
% 
% yyaxis right
% % errorbar(ap_axis,imgaussfilt(predicted_eve_profile_mean(end,:)*1e-5,1),predicted_eve_profile_ste(end,:)*1e-5,'Color',mRNA_red,'CapSize',0);
% %plot(ap_axis,imgaussfilt(predicted_eve_profile_mean(time_plot_1,:)*1e-5,1),'Color',brighten(mRNA_red,0),'LineWidth',1.5) % NL: applying mild smoothing since this is illustrative (not quantitative)
% % s = scatter(ap_axis,imgaussfilt(predicted_eve_profile_mean(end,:)*1e-5,1),50,'MarkerFaceColor',mRNA_red,'MarkerEdgeColor','k');
% eve_mRNA_sm = imgaussfilt(predicted_eve_profile_mean(time_plot_1,:)*1e-5,1);
% f = fill([ap_axis fliplr(ap_axis)], [eve_mRNA_sm zeros(size(eve_mRNA_sm))],mRNA_red);
% f.FaceAlpha = 0.5;
% set(gca,'YColor',mRNA_red)
% xlabel('AP position (% embryo length)');
% ylabel('accumulated {\it eve} mRNA (au)');
% ylim([0 7.5])
% 
% 
% %grid on
% set(gca,'FontSize',14)
% %set(gca,'Color',[228,221,209]/255) 
% % ylim([0.9 1.1])
% 
% xlim([ap_axis(1) ap_axis(end)])
% mRNA_fig1.InvertHardcopy = 'off';
% set(gcf,'color','w'); 
% pbaspect([3 2 1])
% 
% saveas(mRNA_fig1,[FigurePath 'figure2_mRNA_fig_time1.png'])
% saveas(mRNA_fig1,[FigurePath 'figure2_mRNA_fig_time1.pdf'])
% 
% 
% % mRNA profile plot
% mRNA_fig2 = figure;
% hold on
% 
% yyaxis left
% f = fill([ap_axis fliplr(ap_axis)], [knirps_time_array_mean(time_plot_2,:)*1e-5 zeros(size(knirps_time_array_mean(time_plot_1,:)))],color_green);
% %plot(ap_axis,knirps_time_array_mean(time_plot_2,:)*1e-5,'Color',color_green,'LineWidth',3)
% %s = scatter(ap_axis,knirps_time_array_mean(time_plot_2,:)*1e-5,50,'MarkerFaceColor',color_green,'MarkerEdgeColor','k');
% f.FaceAlpha = 0.5;
% ylabel('[Knirps] (au)');
% set(gca,'YColor',color_green)
% ylim([0 15])
% 
% yyaxis right
% % errorbar(ap_axis,imgaussfilt(predicted_eve_profile_mean(end,:)*1e-5,1),predicted_eve_profile_ste(end,:)*1e-5,'Color',mRNA_red,'CapSize',0);
% %plot(ap_axis,imgaussfilt(predicted_eve_profile_mean(time_plot_2,:)*1e-5,1),'Color',brighten(mRNA_red,0),'LineWidth',1.5) % NL: applying mild smoothing since this is illustrative (not quantitative)
% %s = scatter(ap_axis,imgaussfilt(predicted_eve_profile_mean(time_plot_2,:)*1e-5,1),50,'MarkerFaceColor',mRNA_red,'MarkerEdgeColor','k');
% eve_mRNA_sm = imgaussfilt(predicted_eve_profile_mean(time_plot_2,:)*1e-5,1);
% f = fill([ap_axis fliplr(ap_axis)], [eve_mRNA_sm zeros(size(eve_mRNA_sm))],mRNA_red);
% f.FaceAlpha = 0.5;
% set(gca,'YColor',mRNA_red)
% xlabel('AP position (% embryo length)');
% ylabel('accumulated {\it eve} mRNA (au)');
% ylim([0 7.5])
% 
% 
% %grid on
% set(gca,'FontSize',14)
% %set(gca,'Color',[228,221,209]/255) 
% % ylim([0.9 1.1])
% 
% xlim([ap_axis(1) ap_axis(end)])
% mRNA_fig2.InvertHardcopy = 'off';
% set(gcf,'color','w'); 
% pbaspect([3 2 1])
% 
% saveas(mRNA_fig2,[FigurePath 'figure2_mRNA_fig_time2.png'])
% saveas(mRNA_fig2,[FigurePath 'figure2_mRNA_fig_time2.pdf'])