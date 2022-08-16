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

analysis_range = 9;

%binNum = 10;%14
%binMax = 8;

edges = linspace(0,binMax,binNum);

%% plot distribution
%a = response_time_WT(response_time_WT<=8);
a = response_time_WT(response_time_WT<=analysis_range);
b = response_time_HDAC(response_time_HDAC<=analysis_range);

response_time_WT_fig = figure;
h1 = histogram(a,edges,'Normalization','probability');
h1.FaceColor = mRNA_red;
%plot(x,y1,'LineWidth',3,'Color',mRNA_red)
mean(a)

xlim([0 analysis_range])
%xlim([0 8])
ylim([0 0.35])
xlabel('response time (min)')
ylabel('probability')
pbaspect([3 2 1])


response_time_HDAC_fig = figure;
h2 = histogram(b,edges,'Normalization','probability');
h2.FaceColor = mRNA_red;
%plot(x,y1,'LineWidth',3,'Color',mRNA_red)
mean(b)

xlim([0 analysis_range])
ylim([0 0.35])
xlabel('response time (min)')
ylabel('probability')
pbaspect([3 2 1])

%saveas(response_time_WT_HDAC_fig,'figure_response_time_hist.pdf')

%% plot both distribution together

% violin plot
data = [a b];
cat = [ones(1,length(a)) 2*ones(1,length(b))];

response_time_WT_HDAC_violin_fig = figure;
vs = violinplot(data,cat,'ShowData',false);
ylim([0 analysis_range])
ylabel(['response time (min)'])
pbaspect([3 4 1])
%pbaspect([3 2 1])

% smoothened density
[f1,xi] = ksdensity(response_time_WT); 
[f2,x2] = ksdensity(response_time_HDAC);
figure
plot(xi,f1);
hold on
plot(x2,f2);

xlim([0 analysis_range])

saveas(response_time_WT_HDAC_violin_fig,'figure_response_time_WT_HDAC_violin.pdf')



%% Two-sample Kolmogorov-Smirnov test
[h p] = kstest2(response_time_WT,response_time_HDAC);

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