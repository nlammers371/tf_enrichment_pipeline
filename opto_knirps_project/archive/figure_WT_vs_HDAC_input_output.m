clear all
close all
clc

%%
% Define some other colors  
yw = [234 194 100]/255; % yellow
bl = [115 143 193]/255; % blue
gr = [191 213 151]/255; % green

load('frac_on_WT_vs_HDAC.mat');

%%

fig = figure;
hold on
errorbar(knirps_axis,frac_on_WT,frac_on_ste_WT,'Color','k','CapSize',0)
plot(knirps_axis,frac_on_WT,'-k')
scatter(knirps_axis,frac_on_WT,50,'MarkerFaceColor',bl,'MarkerEdgeColor','k')

errorbar(knirps_axis,frac_on_HDAC_3_5,frac_on_HDAC_ste_3_5,'Color','k','CapSize',0)
plot(knirps_axis,frac_on_HDAC_3_5,'-k')
scatter(knirps_axis,frac_on_HDAC_3_5,50,'MarkerFaceColor',yw,'MarkerEdgeColor','k')
%set(gca, 'YScale', 'log')
xlim([2 16])
pbaspect([3 2 1])

%% cdf histogram

cdf_figure = figure;
cdfplot(off_knirps_WT*1e-5)
%h1 = histogram(off_knirps_WT*1e-5,nbins,'Normalization','cdf')
hold on
cdfplot(off_knirps_HDAC_3_5*1e-5)
%h2 = histogram(off_knirps_HDAC_3_5*1e-5,nbins,'Normalization','cdf')
xlim([3 15])
pbaspect([2 1 1])


nbins = 12;
edges = linspace(3,15,nbins);

hist_figure = figure;
h1 = histogram(off_knirps_WT*1e-5,edges,'Normalization','probability')
hold on
%hist_figure_HDAC = figure;
h2 = histogram(off_knirps_HDAC_3_5*1e-5,edges,'Normalization','probability')

xlim([2 15])

%% estimate density

[f1,xi1] = ksdensity(off_knirps_WT*1e-5);
[f2,xi2] = ksdensity(off_knirps_HDAC_3_5*1e-5);

density_fig = figure;
plot(xi1,f1)
hold on
plot(xi2,f2)
xlim([3 15])

%% violin plot

% violin plot
data = [off_knirps_WT*1e-5 off_knirps_HDAC_3_5*1e-5];
cat = [ones(1,length(off_knirps_WT)) 2*ones(1,length(off_knirps_HDAC_3_5))];

violin_fig = figure;

vs1 = violinplot(data,cat,'ShowData',true);
ylim([2 15])
ylabel('[Knirps] (AU)')
pbaspect([2 1 1])

violin__flipped_fig = figure;
vs1 = violinplot(data,cat,'ShowData',false);
ylabel('[Knirps] (AU)')
set(gca,'View',[90 -90])
ylim([2 16])
pbaspect([2 3 1])

%% Two-sample Kolmogorov-Smirnov test
[h1 p1] = kstest2(off_knirps_WT,off_knirps_HDAC_3_5);

