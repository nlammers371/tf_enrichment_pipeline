clear
close all
clc

%addpath(genpath('./lib'))

% color to be used
k_green = brighten([38 142 75]/256,.4);
color_green = [38 143 75]/256; % color from Jake
mRNA_red = brighten([212 100 39]/256,.2);

%% try histogram for data from two sets

load('data.mat')

a = a(a<7);

binNum = 13;
binMax = 6;

%edges1 = linspace(0,10,11); % bin for histogram
%edges1 = linspace(0,15,16);
edges1 = linspace(0,6,binNum);

% fit gamma function
[muhat,muci] = mle(a,'distribution','gamma'); % Generic function
%[muhat,muci] = gamfit(a); % Distribution specific function

x = 0:0.1:20;
y1 = gampdf(x,muhat(1),muhat(2))/(binNum/binMax);
%y2 = gamcdf(x,muhat(1),muhat(2));

fig = figure(2);
hold on
h = histogram(a,edges1,'Normalization','probability');
h.FaceColor = mRNA_red;
plot(x,y1,'LineWidth',3,'Color',mRNA_red)
mean(a)

xlim([0 7])
%xlim([-10 10])
xlabel('response time (min)')
ylabel('probability')
