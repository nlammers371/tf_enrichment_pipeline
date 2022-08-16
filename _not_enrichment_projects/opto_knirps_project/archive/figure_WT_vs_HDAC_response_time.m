clear all
close all
clc

%% Initialization

load('response_time_WT_vs_HDAC.mat')

% color to be used
k_green = brighten([38 142 75]/256,.4);
color_green = [38 143 75]/256; % color from Jake
mRNA_red = brighten([212 100 39]/256,.2);

% Define some other colors  
yw = [234 194 100]/255; % yellow
bl = [115 143 193]/255; % blue
gr = [191 213 151]/255; % green

% histogram parameters
binNum = 27;% best
binMax = 16;% best

analysis_range_ON = 9;
analysis_range_OFF = 15;

%binNum = 10;%14
%binMax = 8;

edges = linspace(0,binMax,binNum);

%% plot response distribution for ON

a1 = response_time_ON_WT(response_time_ON_WT<=analysis_range_ON);
b1 = response_time_ON_HDAC(response_time_ON_HDAC<=analysis_range_ON);

response_time_WT_fig = figure;
h1 = histogram(a1,edges,'Normalization','probability');
h1.FaceColor = mRNA_red;
%plot(x,y1,'LineWidth',3,'Color',mRNA_red)
mean(a1)

xlim([0 analysis_range_ON])
%xlim([0 8])
ylim([0 0.35])
xlabel('response time (min)')
ylabel('probability')
pbaspect([3 2 1])


response_time_HDAC_fig = figure;
h2 = histogram(b1,edges,'Normalization','probability');
h2.FaceColor = mRNA_red;
%plot(x,y1,'LineWidth',3,'Color',mRNA_red)
mean(b1)

xlim([0 analysis_range_ON])
ylim([0 0.35])
xlabel('response time (min)')
ylabel('probability')
pbaspect([3 2 1])

%saveas(response_time_WT_HDAC_fig,'figure_response_time_hist.pdf')

%% plot response distribution for OFF

%a = response_time_WT(response_time_WT<=8);
a2 = response_time_OFF_WT(response_time_OFF_WT<=analysis_range_OFF);
b2 = response_time_OFF_HDAC(response_time_OFF_HDAC<=analysis_range_OFF);

response_time_WT_fig = figure;
h1 = histogram(a2,edges,'Normalization','probability');
h1.FaceColor = mRNA_red;
%plot(x,y1,'LineWidth',3,'Color',mRNA_red)
mean(a2)

xlim([0 analysis_range_OFF])
%xlim([0 8])
ylim([0 0.35])
xlabel('response time (min)')
ylabel('probability')
pbaspect([3 2 1])


response_time_HDAC_fig = figure;
h2 = histogram(b2,edges,'Normalization','probability');
h2.FaceColor = mRNA_red;
%plot(x,y1,'LineWidth',3,'Color',mRNA_red)
mean(b2)

xlim([0 analysis_range_OFF])
ylim([0 0.35])
xlabel('response time (min)')
ylabel('probability')
pbaspect([3 2 1])

%saveas(response_time_WT_HDAC_fig,'figure_response_time_hist.pdf')

%% plot both distribution together

% violin plot
data_ON = [a1 b1];
cat_ON = [ones(1,length(a1)) 2*ones(1,length(b1))];

data_OFF = [a2 b2];
cat_OFF = [ones(1,length(a2)) 2*ones(1,length(b2))];


response_time_WT_HDAC_violin_fig = figure;

tiledlayout(1,2)
nexttile
vs1 = violinplot(data_ON,cat_ON,'ShowData',false);
ylim([0 analysis_range_ON])
ylabel('response time (min)')
pbaspect([3 4 1])

nexttile
vs2 = violinplot(data_OFF,cat_OFF,'ShowData',false);
ylim([0 analysis_range_OFF])
ylabel('response time (min)')
pbaspect([3 4 1])
%pbaspect([3 2 1])

saveas(response_time_WT_HDAC_violin_fig,'figure_response_time_WT_HDAC_violin.pdf')

%%
%{
% smoothened density
[f1,xi] = ksdensity(response_time_WT); 
[f2,x2] = ksdensity(response_time_HDAC);
figure
plot(xi,f1);
hold on
plot(x2,f2);

xlim([0 analysis_range_ON])

saveas(response_time_WT_HDAC_violin_fig,'figure_response_time_WT_HDAC_violin.pdf')
%}


%% Two-sample Kolmogorov-Smirnov test
[h1 p1] = kstest2(response_time_ON_WT,response_time_ON_HDAC);
[h2 p2] = kstest2(response_time_OFF_WT,response_time_OFF_HDAC);

%% fit gamma distribution
%{
% fit gamma function
[muhat,muci] = mle(a,'distribution','gamma'); % Generic function
%[muhat,muci] = gamfit(a); % Distribution specific function

x = 0:0.1:20;
y1 = gampdf(x,muhat(1),muhat(2))/(binNum/binMax);
%y2 = gamcdf(x,muhat(1),muhat(2));

response_time_fig = figure;
hold on
h = histogram(a,edges,'Normalization','probability');
h.FaceColor = mRNA_red;
%plot(x,y1,'LineWidth',3,'Color',mRNA_red)
mean(a)

%xlim([0 analysis_range])
xlim([0 8])
ylim([0 0.275])
xlabel('response time (min)')
ylabel('probability')
pbaspect([3 2 1])
%}