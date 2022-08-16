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
ap_lim = 0.02;
eYFP_background = 375698.13;

projectName = {'optokni_eve4+6_WT','optokni_eve4+6_ON_CONST','optokni_eve4+6_ON_LOW_FULL'}; 

for i = 1:length(projectName)

liveProject = LiveEnrichmentProject(projectName{i});
resultsRoot = [liveProject.dataPath filesep];    

% load data
load([resultsRoot 'spot_struct.mat'])
%FigurePath = [liveProject.figurePath 'WT_vs_CONST' filesep];
%mkdir(FigurePath)

% color to be used
k_green = brighten([38 142 75]/256,.4);
color_green = [38 143 75]/256; % color from Jake
mRNA_red = brighten([212 100 39]/256,.2);

%% Figure: plot mean fluorescence vs time (not aligned)

time_full_long = [];
fluo_full_long = [];
knirps_full_long = [];
frame_full_long = [];

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

    if (spot_struct(j).TraceQCFlag == 1)
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

    end

end


% calculate mean knirps and fraction on

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

if (projectName{i} == "optokni_eve4+6_ON_CONST")
    knirps_vec_mean = convert_from_458(knirps_vec_mean);
end

knirps_vec_mean = knirps_vec_mean-eYFP_background;
frac_on_movmean = movmean(frac_on,5);
fluo_vec_movmean = movmean(fluo_vec_mean,10);

if (projectName{i} == "optokni_eve4+6_WT")
    WT.time_vec = time_vec;
    WT.knirps_vec_mean = knirps_vec_mean;
    WT.frac_on = frac_on;
    WT.frac_on_movmean = frac_on_movmean;
    WT.fluo_vec_mean = fluo_vec_mean;
    WT.fluo_vec_movmean = fluo_vec_movmean;
else
    if (projectName{i} == "optokni_eve4+6_ON_CONST")
        HIGH.time_vec = time_vec;
        HIGH.knirps_vec_mean = knirps_vec_mean;
        HIGH.frac_on = frac_on;
        HIGH.frac_on_movmean = frac_on_movmean;
        HIGH.fluo_vec_mean = fluo_vec_mean;
        HIGH.fluo_vec_movmean = fluo_vec_movmean;
    else
        if (projectName{i} == "optokni_eve4+6_ON_LOW_FULL")
            LOW.time_vec = time_vec;
            LOW.knirps_vec_mean = knirps_vec_mean;
            LOW.frac_on = frac_on;
            LOW.frac_on_movmean = frac_on_movmean;
            LOW.fluo_vec_mean = fluo_vec_mean;
            LOW.fluo_vec_movmean = fluo_vec_movmean;
        end
    end
end

%{
% plot fraction on
temp_traj_fig  = figure('Position',[10 10 800 800]);
tiledlayout(2,1)
nexttile
hold on
%time_interp = min(time_vec_on):0.1:max(time_vec_on);
time_interp = -10:0.1:10;
frac_on_mean = movmean(frac_on,5);
%frac_on_interp = interp1(time_vec(~isnan(frac_on_mean)),frac_on_mean(~isnan(frac_on_mean)),time_interp,'spline');
%frac_on_interp = movmean(frac_on_interp,5);
%frac_on_interp = interp1(time_vec_on,frac_on,time_interp,'v5cubic');

%errorbar(time_vec_on,knirps_vec_ste,'Color','k','CapSize',0);
plot(time_vec,knirps_vec_mean,'-k','LineWidth',1)
scatter(time_vec,knirps_vec_mean,50,'MarkerFaceColor',k_green,'MarkerEdgeColor','k')
%xlim([-10 7])
%ylim([3.75E5 9E5])
xlabel(['time (min) into nc14'])
ylabel(['Knirps concentration (AU)'])
pbaspect([3 2 1])

nexttile
hold on
plot(time_vec,frac_on,'.')
%plot(time_interp,frac_on_interp,'-','LineWidth',2);
scatter(time_vec,frac_on,50,'MarkerFaceColor',mRNA_red,'MarkerEdgeColor','k')
%xlim([-10 7])
%ylim([0 1])
xlabel(['time (min) into nc14'])
ylabel(['fraction of nuclei on'])
pbaspect([3 2 1])
%}


end

%% plot fraction on
%frac_on_interp = interp1(time_vec(~isnan(frac_on_mean)),frac_on_mean(~isnan(frac_on_mean)),time_interp,'spline');
%frac_on_interp = movmean(frac_on_interp,5);
%frac_on_interp = interp1(time_vec_on,frac_on,time_interp,'v5cubic');


temp_traj_fig  = figure('Position',[10 10 800 800]);
tiledlayout(2,1)
nexttile
hold on
%errorbar(time_vec_on,knirps_vec_ste,'Color','k','CapSize',0);
%plot(WT.time_vec,WT.knirps_vec_mean,'-k','LineWidth',1)
plot(WT.time_vec,WT.knirps_vec_mean,'-','LineWidth',1)
plot(HIGH.time_vec,HIGH.knirps_vec_mean,'-','LineWidth',1)
plot(LOW.time_vec,LOW.knirps_vec_mean,'-','LineWidth',1)
%scatter(WT.time_vec,WT.knirps_vec_mean,50,'MarkerFaceColor',k_green,'MarkerEdgeColor','k')
xlim([14.5 30])
%ylim([0 11E5])
xlabel(['time (min) into nc14'])
ylabel(['Knirps concentration (AU)'])
pbaspect([3 2 1])


nexttile
hold on
plot(WT.time_vec,WT.frac_on_movmean,'. -')
plot(HIGH.time_vec,HIGH.frac_on_movmean,'. -')
plot(LOW.time_vec,LOW.frac_on_movmean,'. -')
%plot(WT.time_vec,WT.frac_on,'. -')
%plot(HIGH.time_vec,HIGH.frac_on,'. -')
%plot(LOW.time_vec,LOW.frac_on,'. -')
%plot(time_interp,frac_on_interp,'-','LineWidth',2);
%scatter(WT.time_vec,WT.frac_on,50,'MarkerFaceColor',mRNA_red,'MarkerEdgeColor','k')

xlim([14.5 30])
ylim([0 1])
xlabel(['time (min) into nc14'])
ylabel(['fraction of nuclei on'])
pbaspect([3 2 1])

%% plot mean rate

temp_traj_fig  = figure('Position',[10 10 800 800]);
tiledlayout(2,1)
nexttile
hold on
%errorbar(time_vec_on,knirps_vec_ste,'Color','k','CapSize',0);
%plot(WT.time_vec,WT.knirps_vec_mean,'-k','LineWidth',1)
plot(WT.time_vec,WT.knirps_vec_mean,'LineWidth',1)
plot(HIGH.time_vec,HIGH.knirps_vec_mean,'LineWidth',1)
plot(LOW.time_vec,LOW.knirps_vec_mean,'LineWidth',1)
%scatter(WT.time_vec,WT.knirps_vec_mean,50,'MarkerFaceColor',k_green,'MarkerEdgeColor','k')
xlim([14.5 30])
%ylim([3.75E5 9E5])
xlabel(['time (min) into nc14'])
ylabel(['Knirps concentration (AU)'])
pbaspect([2 1 1])


nexttile
hold on
%plot(WT.time_vec,WT.fluo_vec_mean,'.')
%plot(HIGH.time_vec,HIGH.fluo_vec_mean,'.')
%plot(LOW.time_vec,LOW.fluo_vec_mean,'.')
scatter(WT.time_vec,WT.fluo_vec_mean,25,'filled')
scatter(HIGH.time_vec,HIGH.fluo_vec_mean,25,'filled')
scatter(LOW.time_vec,LOW.fluo_vec_mean,25,'filled')

plot(WT.time_vec,WT.fluo_vec_movmean,'-')
plot(HIGH.time_vec,HIGH.fluo_vec_movmean,'-')
plot(LOW.time_vec,LOW.fluo_vec_movmean,'-')
%plot(time_interp,frac_on_interp,'-','LineWidth',2);
%scatter(WT.time_vec,WT.frac_on,50,'MarkerFaceColor',mRNA_red,'MarkerEdgeColor','k')
%scatter(HIGH.time_vec,HIGH.frac_on,50,'MarkerFaceColor',mRNA_red,'MarkerEdgeColor','k')
%scatter(LOW.time_vec,LOW.frac_on,50,'MarkerFaceColor',mRNA_red,'MarkerEdgeColor','k')
xlim([14.5 30])
%ylim([0 1])
xlabel(['time (min) into nc14'])
ylabel(['mean transcription rate (AU)'])
pbaspect([2 1 1])

%% plot "input-output" relationship for all conditions

fluo_max = 1.8E5;
xRange = linspace(0,15E5,1000);

filter_WT = (WT.time_vec>=5) & (WT.time_vec<=30);
filter_HIGH = (HIGH.time_vec>=15) & (HIGH.time_vec<=30);
filter_LOW = (LOW.time_vec>=15) & (LOW.time_vec<=30);

X = [WT.knirps_vec_mean(filter_WT)' HIGH.knirps_vec_mean(filter_HIGH)' LOW.knirps_vec_mean(filter_LOW)'];
Y = [WT.fluo_vec_mean(filter_WT)' HIGH.fluo_vec_mean(filter_HIGH)' LOW.fluo_vec_mean(filter_LOW)'];

% fit hill function
objective_fun = @(x) x(1)^x(2) ./ (x(1)^x(2) + X.^x(2)) - Y;
options = optimoptions('lsqnonlin','Display','off');
x = lsqnonlin(objective_fun,[2E5 5],[1E5 1],[Inf Inf],options); 

%objective_fun = @(x) x(3)*(x(1)^x(2) ./ (x(1)^x(2) + X.^x(2))) - Y;
%x = lsqnonlin(objective_fun,[9E5 6 1.8E5],[6E5 5 1.6E5],[10E5 8 2E5],options); 

HM = x(1);
hill_coefficient = x(2); 
%hill_function = @(t)(x(1)^x(2) ./ (x(1)^x(2) + t.^x(2)))*x(3);
hill_function = @(t)(x(1)^x(2) ./ (x(1)^x(2) + t.^x(2)))*fluo_max;
yResult = hill_function(xRange);

mean_input_output_fig = figure;
hold on
%plot(WT.knirps_vec_mean(filter_WT),WT.fluo_vec_movmean(filter_WT),'.')
%plot(HIGH.knirps_vec_mean(filter_HIGH),HIGH.fluo_vec_movmean(filter_HIGH),'.')
%plot(LOW.knirps_vec_mean(filter_LOW),LOW.fluo_vec_movmean(filter_LOW),'.')
%plot(WT.knirps_vec_mean(filter_WT),WT.fluo_vec_mean(filter_WT),'.')
%plot(HIGH.knirps_vec_mean(filter_HIGH),HIGH.fluo_vec_mean(filter_HIGH),'.')
%plot(LOW.knirps_vec_mean(filter_LOW),LOW.fluo_vec_mean(filter_LOW),'.')
scatter(WT.knirps_vec_mean(filter_WT),WT.fluo_vec_mean(filter_WT),25,'filled')
scatter(HIGH.knirps_vec_mean(filter_HIGH),HIGH.fluo_vec_mean(filter_HIGH),25,'filled')
scatter(LOW.knirps_vec_mean(filter_LOW),LOW.fluo_vec_mean(filter_LOW),25,'filled')
plot(xRange,yResult,'-')
xlabel('Knirps concentration (AU)')
ylabel('mean transcription rate (AU)')
xlim([2E5 12E5])
pbaspect([3 2 1])


frac_on_input_output_fig = figure;
hold on
%plot(WT.knirps_vec_mean(filter_WT),WT.frac_on_movmean(filter_WT),'.')
%plot(HIGH.knirps_vec_mean(filter_HIGH),HIGH.frac_on_movmean(filter_HIGH),'.')
%plot(LOW.knirps_vec_mean(filter_LOW),LOW.frac_on_movmean(filter_LOW),'.')
plot(WT.knirps_vec_mean(filter_WT),WT.frac_on(filter_WT),'.')
plot(HIGH.knirps_vec_mean(filter_HIGH),HIGH.frac_on(filter_HIGH),'.')
plot(LOW.knirps_vec_mean(filter_LOW),LOW.frac_on(filter_LOW),'.')
xlabel('Knirps concentration (AU)')
ylabel('fraction of cells on')
xlim([2E5 12E5])
ylim([0 1])
pbaspect([3 2 1])



