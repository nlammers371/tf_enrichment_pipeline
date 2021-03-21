clear
close all
clc

addpath(genpath('./lib'))

%% Initialization

projectName = 'optokni_eveBAC_ON'; 

liveProject = LiveEnrichmentProject(projectName);
resultsRoot = [liveProject.dataPath filesep];

% load data
load([resultsRoot 'spot_struct.mat'])
FigurePath = [liveProject.figurePath 'optogenetics' filesep];
mkdir(FigurePath)

% Choose embryo10 as a test case

%expID = 1;
%time_on = 5;
%frame_on = 12;

expID = 2;
time_on = 20.48;
frame_on = 19;

%expID = 3;
%time_on = 27.17;
%frame_on = 40;

% color to be used
k_green = brighten([38 142 75]/256,.4);
color_green = [38 143 75]/256; % color from Jake
mRNA_red = brighten([212 100 39]/256,.2);

knirps_offset = 2.5e5;%prctile(double(knirps_vec_long),1);

ap_lim = 0.02; % AP range for analysis
edges = linspace(0,6,13); % bin for histogram

correction_factor = 1.4144;

%% Figure: plot mean fluorescence vs time (not aligned)

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

for i = 1:length(spot_struct)
    
    if (spot_struct(i).TraceQCFlag == 1) && (spot_struct(i).setID == expID)
        % extract core vectors 
        
        % extract core vectors 
        fluo_vec_orig = spot_struct(i).fluo;
        time_vec_orig = spot_struct(i).time;
        frame_vec_orig = spot_struct(i).frames;
        knirps_vec_orig = spot_struct(i).rawNCProtein;
        ap_vec_orig = spot_struct(i).apPosNucleus;
             
        % calculate mean
        ever_on_orig = any(~isnan(fluo_vec_orig));
        mean_ap_orig = nanmean(ap_vec_orig);
        mean_knirps_orig = nanmean(knirps_vec_orig);
        
        if ever_on_orig
            last_on_frame = frame_vec_orig(find(~isnan(fluo_vec_orig) & (frame_vec_orig <= frame_on),1,'last'));
            first_on_frame = frame_vec_orig(find(~isnan(fluo_vec_orig) & (frame_vec_orig > frame_on),1));
        end
        
        if (mean_ap_orig > -ap_lim) && (mean_ap_orig < ap_lim)
           time_orig_long = [time_orig_long time_vec_orig];
           frame_orig_long = [frame_orig_long frame_vec_orig];
           fluo_orig_long = [fluo_orig_long fluo_vec_orig];
           knirps_orig_long = [knirps_orig_long knirps_vec_orig-knirps_offset];
           
           if ~isnan(last_on_frame) & ~isnan(first_on_frame)
               last_on_long = [last_on_long last_on_frame];
               first_on_long = [first_on_long first_on_frame];
           end
           
           %plot((time_vec_orig)/60,fluo_vec_orig,'Color', [175 175 175]/255);
           count = count + 1 
        end
        
%         if (mean_ap(end) > -0.01) && (mean_ap(end) < 0.01)
%             plot(time_vec-time_vec(end),fluo_vec);
%             count = count + 1
%         end
        
    end
    
end

%% calculate mean knirps and fraction on (before/after perturbation)

fluo_orig_long(isnan(fluo_orig_long)) = 0;

time_bin = 2:max(frame_orig_long);
time_groups = discretize(frame_orig_long,time_bin);

fluo_vec_mean = zeros(length(time_bin)-1,1);
fluo_vec_ste = zeros(length(time_bin)-1,1);

knirps_vec_mean = zeros(length(time_bin)-1,1);
knirps_vec_ste = zeros(length(time_bin)-1,1);

frac_on = zeros(length(time_bin)-1,1);
frac_on_ste = zeros(length(time_bin)-1,1);

fluo_orig_long_zero = fluo_orig_long;
fluo_orig_long_zero(isnan(fluo_orig_long)) = 0;

fluo_orig_long_binary = fluo_orig_long;
fluo_orig_long_binary(isnan(fluo_orig_long)) = 0;
fluo_orig_long_binary(fluo_orig_long>0) = 1;


for i = 1:length(time_bin)-1
 
    time_filter_long = time_groups==i;
    
    %if sum(time_filter_long) > 10        
    %    boot_samples_fluo = bootstrp(nBoots,@nanmean,fluo_orig_long(time_filter_long));
    %    fluo_vec_mean(i) = nanmean(boot_samples_fluo);
    %    fluo_vec_ste(i) = std(boot_samples_fluo);
    
    time_vec(i) = mean(time_orig_long(time_filter_long))/60;
    
    fluo_vec_mean(i) = nanmean(fluo_orig_long_zero(time_filter_long));
    fluo_vec_ste(i) = std(fluo_orig_long_zero(time_filter_long),'omitnan');
    
    knirps_vec_mean(i) = nanmean(knirps_orig_long(time_filter_long));
    knirps_vec_ste(i) = std(knirps_orig_long(time_filter_long),'omitnan');
    
    frac_on(i) = nanmean(fluo_orig_long_binary(time_filter_long));
    frac_on_ste(i) = std(fluo_orig_long_binary(time_filter_long));
    
end

time_vec_on = time_vec-time_vec(frame_on);

knirps_vec_mean(time_vec_on>=0) = knirps_vec_mean(time_vec_on>=0)/correction_factor;

fig  = figure(1);
hold on

%time_interp = min(time_vec_on):0.1:max(time_vec_on);
time_interp = -10:0.1:10;
frac_on_mean = movmean(frac_on,5);
frac_on_interp = interp1(time_vec_on(~isnan(frac_on_mean)),frac_on_mean(~isnan(frac_on_mean)),time_interp,'spline');
frac_on_interp = movmean(frac_on_interp,5);
%frac_on_interp = interp1(time_vec_on,frac_on,time_interp,'v5cubic');

yyaxis left
errorbar(time_vec_on,knirps_vec_ste,'Color','k','CapSize',0);
plot(time_vec_on,knirps_vec_mean,'-k')
scatter(time_vec_on,knirps_vec_mean,50,'MarkerFaceColor',k_green,'MarkerEdgeColor','k')

ylim([3.5E5 8.5E5])

yyaxis right
plot(time_vec_on,frac_on,'.')
plot(time_interp,frac_on_interp,'-','LineWidth',2);
scatter(time_vec_on,frac_on,50,'MarkerFaceColor',mRNA_red,'MarkerEdgeColor','k')
ylim([0 1])
    
xlim([-10 10])
xlabel(['time relative to perturbation (min)'])
ylabel(['fraction of nuclei on'])

%% calculate waiting time

data_filter = (last_on_long <= (frame_on-3));

fig = figure(2);
histogram(time_vec(first_on_long(data_filter))-time_vec(frame_on),edges)
mean(time_vec(first_on_long(data_filter))-time_vec(frame_on))
