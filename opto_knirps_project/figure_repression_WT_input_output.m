clear
close all

addpath(genpath('./lib'))

eYFP_background = 375698.13;

% knirps green
k_green = brighten([38 142 75]/256,.4);
color_green = [38 143 75]/256; % color from Jake
mRNA_red = brighten([212 100 39]/256,.2);

%% Initialization

projectName = 'optokni_eve4+6_WT'; 

liveProject = LiveEnrichmentProject(projectName);
resultsRoot = [liveProject.dataPath filesep];

% load data
load([resultsRoot 'spot_struct.mat'])
FigurePath = [liveProject.figurePath 'input_output' filesep];
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

for i = 1:length(spot_struct)
  
    % extract core vectors 
    fluo_vec = spot_struct(i).fluo;
    time_vec = spot_struct(i).time;
    knirps_vec = spot_struct(i).rawNCProtein;
    ap_vec = spot_struct(i).APPosNucleus;
    
    if time_vec(end) - time_vec(1) >=30*60
        
        % make average vectors        
      
        % get off and on indices
        ever_on_vec(i) = any(~isnan(fluo_vec));
        mean_ap(i) = nanmean(ap_vec);
        mean_knirps(i) = nanmean(knirps_vec);
        
        post_on_vec = zeros(size(ap_vec));
        post_off_vec = zeros(size(ap_vec));
        if ever_on_vec(i)
            start_i = find(~isnan(fluo_vec),1);
            stop_i = find(~isnan(fluo_vec),1,'last');
            if true%stop_i < length(fluo_vec)-10
                off_time_vec(i) = time_vec(stop_i);
                off_knirps_vec(i) = knirps_vec(stop_i);
                off_spot_fluo(i) = fluo_vec(stop_i);
                off_ap(i) = ap_vec(stop_i);
                post_off_vec(stop_i+1:end) = 1;
            end
            
            if start_i > 1
                on_time_vec(i) = time_vec(start_i);
                on_knirps_vec(i) = knirps_vec(start_i);
                on_ap(i) = ap_vec(start_i);
                post_on_vec(start_i+1:end) = 1;
            end
        end
        
        % make regression vectors
        ever_on_flags = [ever_on_flags repelem(ever_on_vec(i),length(ap_vec))];
        post_turn_on_flags = [post_turn_on_flags post_on_vec];
        post_turn_off_flags = [post_turn_off_flags post_off_vec];
        fluo_raw_long = [fluo_raw_long fluo_vec];
        
        fluo_zeros = fluo_vec;
        all_zeros = fluo_vec;
        fluo_zeros(post_on_vec&~post_off_vec&isnan(fluo_vec)) = 0;
        fluo_zeros_long = [fluo_zeros_long fluo_zeros];
        
        all_zeros(isnan(all_zeros)) = 0;
        mRNA_vec_long = [mRNA_vec_long all_zeros];
        
        mean_fluo(i) = nanmean(fluo_zeros);
        
        knirps_vec_long_raw = [knirps_vec_long_raw knirps_vec];
        ap_vec_long = [ap_vec_long ap_vec];
        time_vec_long = [time_vec_long time_vec];
        
    end
end

%% Figure 2: Plot mean vectors

knirps_vec_long = knirps_vec_long_raw - eYFP_background;

nBoots = 100;

nLinBins = 41;
ap_bins = linspace(-0.105,0.105,25); % used for input-output function
%ap_bins = linspace(-0.105,0.105,15); % used in mean figure
time_bins = linspace(5*60,40*60,nLinBins);
%ap_bins = linspace(-0.125,0.125,26);
%ap_bins = linspace(-0.12,0.12,21);

time_groups = discretize(time_vec_long,time_bins); 
ap_groups = discretize(ap_vec_long,ap_bins); 
ap_groups_mean = discretize(mean_ap,ap_bins); 
ap_groups_off = discretize(off_ap,ap_bins); 

frac_on_vec_mean = NaN(1,length(ap_bins)-1);
frac_on_vec_ste = NaN(1,length(ap_bins)-1);

off_time_vec_mean = NaN(1,length(ap_bins)-1);
off_time_vec_ste = NaN(1,length(ap_bins)-1);

fluo_vec_mean = NaN(1,length(ap_bins)-1);
fluo_vec_ste = NaN(1,length(ap_bins)-1);

knirps_vec_mean_all = NaN(1,length(ap_bins)-1);


for a = 1:length(ap_bins)-1
    ap_filter_long = ap_groups==a;
    time_filter_long = time_groups>=18;
    if sum(ap_filter_long) > 10        
        boot_samples_fluo = bootstrp(nBoots,@nanmean,fluo_zeros_long(ap_filter_long));
        fluo_vec_mean(a) = mean(boot_samples_fluo);
        fluo_vec_ste(a) = std(boot_samples_fluo);
    end  
    
    ap_filter_mean = ap_groups_mean==a;
    ap_filter_off = ap_groups_off==a;
    if sum(ap_filter_mean) > 10        
        boot_samples_off = bootstrp(nBoots,@nanmean,off_time_vec(ap_filter_off));
        off_time_vec_mean(a) = mean(boot_samples_off);
        off_time_vec_ste(a) = std(boot_samples_off);
        
        boot_samples_frac = bootstrp(nBoots,@nanmean,ever_on_vec(ap_filter_mean));
        frac_on_vec_mean(a) = mean(boot_samples_frac);
        frac_on_vec_ste(a) = std(boot_samples_frac);
    end   
    
    knirps_vec_mean_all(a) = nanmean(knirps_vec_long(ap_filter_long & time_filter_long));
             
end

% Define some colors  
yw = [234 194 100]/255; % yellow
bl = [115 143 193]/255; % blue
gr = [191 213 151]/255; % green

ap_axis = 100*(ap_bins(1:end-1) + diff(ap_bins)/2);

%ap_axis = ap_axis - 61;
% fraction on
fraction_on_fig = figure;
hold on

errorbar(ap_axis,frac_on_vec_mean,frac_on_vec_ste,'Color','k','CapSize',0);
%boundedline(ap_axis',frac_on_vec_mean',frac_on_vec_ste', '- .','nan', 'gap','alpha','cmap',brighten(yw,-0.5));
plot(ap_axis,frac_on_vec_mean,'-k')
scatter(ap_axis,frac_on_vec_mean,50,'MarkerFaceColor',yw,'MarkerEdgeColor','k')

xlabel('AP position (% embryo length)');
ylabel('fraction of active nuclei');

%grid on
set(gca,'FontSize',14)
%set(gca,'Color',[228,221,209]/255) 
xlim([ap_axis(1) ap_axis(end)])
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

fraction_on_fig.InvertHardcopy = 'off';
set(gcf,'color','w'); 
ylim([0 1.1])   
pbaspect([3 2 1])
saveas(fraction_on_fig,[FigurePath 'figure_fraction_on_vs_ap.png'])
saveas(fraction_on_fig,[FigurePath 'figure_fraction_on_vs_ap.pdf'])


% off time
off_time_fig = figure;
hold on

errorbar(ap_axis,off_time_vec_mean/60,off_time_vec_ste/60,'Color','k','CapSize',0);
%boundedline(ap_axis',off_time_vec_mean'/60,off_time_vec_ste'/60, '- .','nan', 'gap','alpha','cmap',brighten(bl,-0.5));
plot(ap_axis,off_time_vec_mean/60,'-k')
scatter(ap_axis,off_time_vec_mean/60,50,'MarkerFaceColor',bl,'MarkerEdgeColor','k')

xlabel('AP position (% embryo length)');
ylabel('average off time (minutes)');

%grid on
set(gca,'FontSize',14)
%set(gca,'Color',[228,221,209]/255) 
ylim([15 40])
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
xlim([ap_axis(1) ap_axis(end)])
off_time_fig.InvertHardcopy = 'off';
set(gcf,'color','w'); 
pbaspect([3 2 1]) 
saveas(off_time_fig,[FigurePath 'figure_off_time_vs_ap.png'])
saveas(off_time_fig,[FigurePath 'figure_off_time_vs_ap.pdf'])


mean_fig = figure;
hold on

errorbar(ap_axis,fluo_vec_mean*1e-5,fluo_vec_ste*1e-5,'Color','k','CapSize',0);
%boundedline(ap_axis',fluo_vec_mean'*1e-5,fluo_vec_ste'*1e-5, '- .','nan', 'gap','alpha','cmap',brighten(gr,-0.5));
plot(ap_axis',fluo_vec_mean'*1e-5,'-k')
scatter(ap_axis,fluo_vec_mean*1e-5,50,'MarkerFaceColor',gr,'MarkerEdgeColor','k')

xlabel('AP position (% embryo length)');
ylabel('mean spot intensity (au)');

%grid on
set(gca,'FontSize',14)
%set(gca,'Color',[228,221,209]/255) 
ylim([0 2.5])
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
xlim([ap_axis(1) ap_axis(end)])
mean_fig.InvertHardcopy = 'off';
set(gcf,'color','w'); 
pbaspect([3 2 1])    
saveas(mean_fig,[FigurePath 'figure_fluo_vs_ap.png'])
saveas(mean_fig,[FigurePath 'figure_fluo_vs_ap.pdf'])

%% plot mean vectors with average knirps

% off time
off_time_fig = figure;
hold on

yyaxis right
%plot(ap_axis, knirps_vec_mean_all);
f = fill([ap_axis fliplr(ap_axis)], [knirps_vec_mean_all*1e-5 zeros(size(knirps_vec_mean_all))],color_green);
%plot(ap_axis,knirps_time_array_mean(time_plot_1,:)*1e-5,'Color',color_green,'LineWidth',3)
f.FaceAlpha = 0.3;
ylabel('Knirps (AU)');
ylim([0 12.5])

yyaxis left
errorbar(ap_axis,off_time_vec_mean/60,off_time_vec_ste/60,'Color','k','CapSize',0);
%boundedline(ap_axis',off_time_vec_mean'/60,off_time_vec_ste'/60, '- .','nan', 'gap','alpha','cmap',brighten(bl,-0.5));
plot(ap_axis,off_time_vec_mean/60,'-k')
scatter(ap_axis,off_time_vec_mean/60,50,'MarkerFaceColor',bl,'MarkerEdgeColor','k')
ylabel('average off time (minutes)');
ylim([10 40])

xlabel('AP position (% embryo length)');

%grid on
set(gca,'FontSize',14)
%set(gca,'Color',[228,221,209]/255) 
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
ax.XAxis(1).Color = 'k';
%xlim([-12 12])
xlim([ap_axis(1) ap_axis(end)])
off_time_fig.InvertHardcopy = 'off';
set(gcf,'color','w'); 
pbaspect([3 2 1])

saveas(off_time_fig,[FigurePath 'figure_off_time_vs_ap_with_kni.png'])
saveas(off_time_fig,[FigurePath 'figure_off_time_vs_ap_with_kni.pdf'])



mean_fig = figure;
hold on

yyaxis right
%plot(ap_axis, knirps_vec_mean_all);
f = fill([ap_axis fliplr(ap_axis)], [knirps_vec_mean_all*1e-5 zeros(size(knirps_vec_mean_all))],color_green);
%plot(ap_axis,knirps_time_array_mean(time_plot_1,:)*1e-5,'Color',color_green,'LineWidth',3)
f.FaceAlpha = 0.3;
ylabel('Knirps (AU)');
ylim([0 12.5])

yyaxis left
errorbar(ap_axis,fluo_vec_mean*1e-5,fluo_vec_ste*1e-5,'Color','k','CapSize',0);
%boundedline(ap_axis',fluo_vec_mean'*1e-5,fluo_vec_ste'*1e-5, '- .','nan', 'gap','alpha','cmap',brighten(gr,-0.5));
plot(ap_axis',fluo_vec_mean'*1e-5,'-k')
scatter(ap_axis,fluo_vec_mean*1e-5,50,'MarkerFaceColor',gr,'MarkerEdgeColor','k')

xlabel('AP position (% embryo length)');
ylabel('mean spot intensity (au)');

%grid on
set(gca,'FontSize',14)
%set(gca,'Color',[228,221,209]/255) 
ylim([0 3])
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
ax.XAxis(1).Color = 'k';
xlim([ap_axis(1) ap_axis(end)])
mean_fig.InvertHardcopy = 'off';
set(gcf,'color','w'); 
pbaspect([3 2 1])    
%saveas(mean_fig,[FigurePath 'figure_fluo_vs_ap_with_kni.png'])
%saveas(mean_fig,[FigurePath 'figure_fluo_vs_ap_with_kni.pdf'])


%% Look at long vectors

nLinBins = 31;
nLinBins_kni = 31;
time_bins = linspace(5*60,40*60,nLinBins);
knirps_bins = linspace(0,15e5,nLinBins_kni);

ap_groups = discretize(ap_vec_long,ap_bins); 
time_groups = discretize(time_vec_long,time_bins); 
knirps_groups = discretize(knirps_vec_long,knirps_bins);

ap_time_groups_num = NaN(length(time_bins)-1,length(ap_bins)-1);

frac_on_time_array_mean = NaN(length(time_bins)-1,length(ap_bins)-1);
eve_time_array_mean = NaN(length(time_bins)-1,length(ap_bins)-1);
eve_time_array_full = NaN(length(time_bins)-1,length(ap_bins)-1,nBoots);
knirps_time_array_mean = NaN(length(time_bins)-1,length(ap_bins)-1);

frac_on_time_array_ste = NaN(length(time_bins)-1,length(ap_bins)-1);
eve_time_array_ste = NaN(length(time_bins)-1,length(ap_bins)-1);
knirps_time_array_ste = NaN(length(time_bins)-1,length(ap_bins)-1);

frac_on_knirps_array_mean = NaN(length(knirps_bins)-1,length(ap_bins)-1);
frac_on_knirps_array_ste = NaN(length(knirps_bins)-1,length(ap_bins)-1);


for t = 1:length(time_bins)
    for a = 1:length(ap_bins)
        time_window_filter = time_groups==t & ap_groups==a & post_turn_on_flags;
        ap_time_groups_num(t,a) = sum(time_window_filter);
        if sum(time_window_filter) > 50
            still_on_flags_time = ~post_turn_off_flags(time_window_filter);
            boot_samples_time = bootstrp(nBoots,@mean,still_on_flags_time);
            frac_on_time_array_mean(t,a) = mean(boot_samples_time);
            frac_on_time_array_ste(t,a) = std(boot_samples_time);
            
            boot_samples_mRNA = bootstrp(nBoots,@nanmean,mRNA_vec_long(time_window_filter));
            eve_time_array_mean(t,a) = mean(boot_samples_mRNA);
            eve_time_array_ste(t,a) = std(boot_samples_mRNA);
            eve_time_array_full(t,a,:) = boot_samples_mRNA; % NL: note that this is not quite right. Should really be sampling at the level of loci 
            
            boot_samples_knirps = bootstrp(nBoots,@nanmean,knirps_vec_long_raw(time_window_filter));
            knirps_time_array_mean(t,a) = nanmean(boot_samples_knirps); % NL: note that this is not quite right. Should really be sampling at the level of loci 
            knirps_time_array_ste(t,a) = nanstd(boot_samples_knirps);
        end

    end
end

for t = 1:length(knirps_bins)
    for a = 1:length(ap_bins)
      
        knirps_window_filter = knirps_groups==t & ap_groups==a & post_turn_on_flags;
        if sum(knirps_window_filter)>50
            still_on_flags_knirps = ~post_turn_off_flags(knirps_window_filter);
            boot_samples_knirps = bootstrp(nBoots,@mean,still_on_flags_knirps);
            frac_on_knirps_array_mean(t,a) = mean(boot_samples_knirps);
            frac_on_knirps_array_ste(t,a) = std(boot_samples_knirps);
        end
    end
end

%% Figure: plot mean fluorescence vs time (aligned by off time)


nBoots = 100;

ever_on_vec = [];
mean_ap = [];
time_orig_long = [];
fluo_orig_long = [];
knirps_orig_long = [];
off_time_long = [];
off_knirps_long = [];


count = 0;

fluo_silencing_fig  = figure;
hold on

for i = 1:length(spot_struct)
    
    if (spot_struct(i).TraceQCFlag == 1)
        % extract core vectors 
        
        % extract core vectors 
        fluo_vec_orig = spot_struct(i).fluo;
        time_vec_orig = spot_struct(i).time;
        knirps_vec_orig = spot_struct(i).rawNCProtein;
        ap_vec_orig = spot_struct(i).APPosNucleus;
        
        % get off and on indices
        ever_on_orig = any(~isnan(fluo_vec_orig));
        mean_ap_orig = nanmean(ap_vec_orig);
        mean_knirps_orig = nanmean(knirps_vec_orig);

        %fluo_vec = spot_struct(i).fluoInterp;
        %time_vec = spot_struct(i).timeInterp;
        %ap_vec = spot_struct(i).APPosParticleInterp;
        
        mean_ap = [mean_ap nanmean(ap_vec)];
        
        if ever_on_orig
            start_i = find(~isnan(fluo_vec_orig),1);
            stop_i = find(~isnan(fluo_vec_orig),1,'last');
            off_time_orig = time_vec_orig(stop_i);
            off_knirps_vec_orig = knirps_vec_orig(stop_i);
            off_spot_fluo_orig = fluo_vec_orig(stop_i);
            off_ap_orig = ap_vec_orig(stop_i);
        end
        
        if (mean_ap_orig > -0.02) && (mean_ap_orig < 0.02)
           % 0.01  is for the figure
           time_orig_long = [time_orig_long time_vec_orig-off_time_orig];
           fluo_orig_long = [fluo_orig_long fluo_vec_orig];
           off_time_long = [off_time_long off_time_orig];
           off_knirps_long = [off_knirps_long off_knirps_vec_orig];
           knirps_orig_long = [knirps_orig_long knirps_vec_orig];
           
           %plot((time_vec_orig-off_time_orig)/60,fluo_vec_orig,'Color', [175 175 175]/255);
           plot((time_vec_orig-off_time_orig)/60,fluo_vec_orig);
           count = count + 1 
        end
        
%         if (mean_ap(end) > -0.01) && (mean_ap(end) < 0.01)
%             plot(time_vec-time_vec(end),fluo_vec);
%             count = count + 1
%         end
        
    end
    
end

time_bin = linspace(-15,0,21);
time_groups = discretize(time_orig_long/60,time_bin);

fluo_vec_mean = zeros(length(time_bin)-1,1);
fluo_vec_ste = zeros(length(time_bin)-1,1);

for i = 1:length(time_bin)-1
 
    time_filter_long = time_groups==i;
    
    %if sum(time_filter_long) > 10        
    %    boot_samples_fluo = bootstrp(nBoots,@nanmean,fluo_orig_long(time_filter_long));
    %    fluo_vec_mean(i) = nanmean(boot_samples_fluo);
    %    fluo_vec_ste(i) = std(boot_samples_fluo);
    
    fluo_vec_mean(i) = nanmean(fluo_orig_long(time_filter_long));
    fluo_vec_ste(i) = std(fluo_orig_long(time_filter_long));
    
    knirps_vec_mean(i) = nanmean(knirps_orig_long(time_filter_long));
    knirps_vec_ste(i) = std(knirps_orig_long(time_filter_long));
    
end  
    

plot(time_bin(2:end),fluo_vec_mean,'LineWidth',5,'Color',mRNA_red)
    
xlim([-10 3])
%ylim([0 3.5E5])
ylim([0 6E5])
xlabel(['time (min) relative to the silencing event'])
ylabel(['mean activity (au)'])

pbaspect([3 2 1])

%saveas(fluo_silencing_fig,[FigurePath 'figure_fluo_silencing.png'])
%saveas(fluo_silencing_fig,[FigurePath 'figure_fluo_silencing.pdf'])


%% Figure: plot mean fluorescence vs time (not aligned)

nBoots = 100;

ever_on_vec = [];
mean_ap = [];
time_orig_long = [];
fluo_orig_long = [];
off_time_long = [];

count = 0;

fig  = figure;
hold on

for i = 1:length(spot_struct)
    
    if (spot_struct(i).TraceQCFlag == 1)
        % extract core vectors 
        
        % extract core vectors 
        fluo_vec_orig = spot_struct(i).fluo;
        time_vec_orig = spot_struct(i).time;
        knirps_vec_orig = spot_struct(i).rawNCProtein;
        ap_vec_orig = spot_struct(i).APPosNucleus;
        
        % get off and on indices
        ever_on_orig = any(~isnan(fluo_vec_orig));
        mean_ap_orig = nanmean(ap_vec_orig);
        mean_knirps_orig = nanmean(knirps_vec_orig);
        

        %fluo_vec = spot_struct(i).fluoInterp;
        %time_vec = spot_struct(i).timeInterp;
        %ap_vec = spot_struct(i).APPosParticleInterp;
        
        mean_ap = [mean_ap nanmean(ap_vec)];
        
        if ever_on_orig
            start_i = find(~isnan(fluo_vec_orig),1);
            stop_i = find(~isnan(fluo_vec_orig),1,'last');
            off_time_orig = time_vec_orig(stop_i);
            off_knirps_vec_orig = knirps_vec_orig(stop_i);
            off_spot_fluo_orig = fluo_vec_orig(stop_i);
            off_ap_orig = ap_vec_orig(stop_i);
        end
        
        if (mean_ap_orig > -0.02) && (mean_ap_orig < 0.02)
           time_orig_long = [time_orig_long time_vec_orig];
           fluo_orig_long = [fluo_orig_long fluo_vec_orig];
           off_time_long = [off_time_long off_time_orig];
           
           plot((time_vec_orig)/60,fluo_vec_orig,'Color', [175 175 175]/255);
           count = count + 1 
        end
        
%         if (mean_ap(end) > -0.01) && (mean_ap(end) < 0.01)
%             plot(time_vec-time_vec(end),fluo_vec);
%             count = count + 1
%         end
        
    end
    
end

fluo_orig_long(isnan(fluo_orig_long)) = 0;

time_bin = linspace(0,50,51);
time_groups = discretize(time_orig_long/60,time_bin);

fluo_vec_mean = zeros(length(time_bin)-1,1);
fluo_vec_ste = zeros(length(time_bin)-1,1);

for i = 1:length(time_bin)-1
 
    time_filter_long = time_groups==i;
    
    %if sum(time_filter_long) > 10        
    %    boot_samples_fluo = bootstrp(nBoots,@nanmean,fluo_orig_long(time_filter_long));
    %    fluo_vec_mean(i) = nanmean(boot_samples_fluo);
    %    fluo_vec_ste(i) = std(boot_samples_fluo);
    
    fluo_vec_mean(i) = nanmean(fluo_orig_long(time_filter_long));
    fluo_vec_ste(i) = std(fluo_orig_long(time_filter_long));
    
end  
    

plot(time_bin(2:end),fluo_vec_mean,'LineWidth',5)
    
xlim([5 35])
ylim([0 2.5E5])
xlabel(['time (min) relative to the silencing event'])

silenceTimefig = figure;
histogram(off_time_long/60)
xlim([5 35])


%% Supp Figure: plot fraction_on vs time

ap_filter = ap_axis>=-4 & ap_axis <=4; % used for input-output function
ap_indices = find(ap_filter);

%ap_cell = {ap_indices(1:3) ap_indices(1:6) ap_indices(1:end)}; %ap_indices(1:9) ap_indices(1:end)};
%ap_cell = {ap_indices(1:3) ap_indices(1:6) ap_indices(1:9) ap_indices(1:end)};
ap_cell = {ap_indices(1:end)};

close all
knirps_axis = knirps_bins(1:end-1) + diff(knirps_bins)/2;
knirps_axis = knirps_axis*1e-5;
time_axis = time_bins(1:end-1) + diff(time_bins)/2;
time_axis = time_axis/60;

cmap_bin_num = 100;

for i = 1:length(ap_cell)
    ap_indices_iter = ap_cell{i};
  
    ap_time_fig = figure;
    hold on
    %cmap1 = brewermap(length(ap_indices),'blues');
    cmap1 = brewermap(cmap_bin_num*1.2,'blues');
    cmap1 = cmap1(1:cmap_bin_num,:);
    
    colormap(cmap1)
    %iter = 1;
    for a = ap_indices_iter
        errorbar(time_axis,frac_on_time_array_mean(:,a),frac_on_time_array_ste(:,a),'Color','k','CapSize',0);
        %scatter(time_axis,frac_on_time_array_mean(:,a),'MarkerFaceColor',cmap1(iter,:),'MarkerEdgeColor','k')
        cmap_bin = round((ap_axis(a)-ap_axis(ap_indices(1)))*(cmap_bin_num-1)/(ap_axis(ap_indices(end))-ap_axis(ap_indices(1))))+1;
        scatter(time_axis,frac_on_time_array_mean(:,a),'MarkerFaceColor',cmap1(cmap_bin,:),'MarkerEdgeColor','k')
        %iter = iter + 1;
    end
    caxis([ap_axis(ap_indices(1)) ap_axis(ap_indices(end))])
    h = colorbar;
    xlim([5 37.5])
    xlabel('time since start of nc14 (minutes)');
    ylabel('fraction of {\it eve} loci still on');
    ylabel(h,'AP position')

    %grid on
    set(gca,'FontSize',14)
    %set(gca,'Color',[228,221,209]/255) 

    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.XAxis(1).Color = 'k';

    ap_time_fig.InvertHardcopy = 'off';
    set(gcf,'color','w'); 

    saveas(ap_time_fig,[FigurePath 'figure_supp_fraction_on_vs_time_' num2str(i) '.png'])
    saveas(ap_time_fig,[FigurePath 'figure_supp_fraction_on_vs_time_' num2str(i) '.pdf'])
end

%% plot fraction_on vs time for one bin

%ap_cell = {ap_indices(1:3) ap_indices(1:6) ap_indices(1:end)}; %ap_indices(1:9) ap_indices(1:end)};
%ap_cell = {ap_indices(1:3) ap_indices(1:6) ap_indices(1:9) ap_indices(1:end)};
%ap_cell = {ap_indices(1:end)};

ap_cell = {round(max(ap_indices)/2)};

close all
knirps_axis = knirps_bins(1:end-1) + diff(knirps_bins)/2;
knirps_axis = knirps_axis*1e-5;
time_axis = time_bins(1:end-1) + diff(time_bins)/2;
time_axis = time_axis/60;

cmap_bin_num = 100;

for i = 1:length(ap_cell)
    ap_indices_iter = ap_cell{i};
  
    ap_time_fig = figure;
    hold on
    %cmap1 = brewermap(length(ap_indices),'blues');
    cmap1 = brewermap(cmap_bin_num*1.2,'blues');
    cmap1 = cmap1(1:cmap_bin_num,:);
    
    colormap(cmap1)
    %iter = 1;
    for a = ap_indices_iter
        errorbar(time_axis,frac_on_time_array_mean(:,a),frac_on_time_array_ste(:,a),'Color','k','CapSize',0);
        %scatter(time_axis,frac_on_time_array_mean(:,a),'MarkerFaceColor',cmap1(iter,:),'MarkerEdgeColor','k')
        %cmap_bin = round((ap_axis(a)-ap_axis(ap_indices(1)))*(cmap_bin_num-1)/(ap_axis(ap_indices(end))-ap_axis(ap_indices(1))))+1;
        scatter(time_axis,frac_on_time_array_mean(:,a),'MarkerFaceColor',cmap1(1,:),'MarkerEdgeColor','k')
        %iter = iter + 1;
    end
    caxis([ap_axis(ap_indices(1)) ap_axis(ap_indices(end))])
    h = colorbar;
    xlim([time_bins(1)/60 time_bins(end)/60])
    xlabel('time since start of nc14 (minutes)');
    ylabel('fraction of {\it eve} loci still on');
    ylabel(h,'AP position')

    %grid on
    set(gca,'FontSize',14)
    %set(gca,'Color',[228,221,209]/255) 

    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.XAxis(1).Color = 'k';

    ap_time_fig.InvertHardcopy = 'off';
    set(gcf,'color','w'); 

    saveas(ap_time_fig,[FigurePath 'figure2_supp_fraction_on_vs_time_' num2str(i) '.png'])
    saveas(ap_time_fig,[FigurePath 'figure2_supp_fraction_on_vs_time_' num2str(i) '.pdf'])
end


%% Plot input-output function

ap_cell = {ap_indices(1:end-1)};

for i = length(ap_cell)
    ap_indices_iter = ap_cell{i};
    
    ap_knirps_plot_fig = figure;
    hold on
    cmap1 = flipud(brewermap(1.2*cmap_bin_num,'Spectral'));
    %cmap1 = flipud(brewermap(1.2*cmap_bin_num,'RdYlGn'));
    cmap1 = cmap1(0.10*cmap_bin_num:1.10*cmap_bin_num-1,:);
    
    colormap(cmap1)
    iter = 1;
    for a = ap_indices_iter
        
        cmap_bin = round((ap_axis(a)-ap_axis(ap_indices(1)))*(cmap_bin_num-1)/(ap_axis(ap_indices(end))-ap_axis(ap_indices(1))))+1;
        %plot(knirps_axis,fit_profile,'Color',brighten(cmap1(cmap_bin,:),-0.5),'LineWidth',0.75)
        %plot(knirps_axis,fit_profile,'Color','k','LineWidth',1)
        errorbar(knirps_axis,frac_on_knirps_array_mean(:,a),frac_on_knirps_array_ste(:,a),'o','Color',[0 0 0 0],'CapSize',0);
    %     plot(knirps_axis,frac_on_knirps_array_mean(:,a),'color',[0 0 0 0.2]);
        scatter(knirps_axis,frac_on_knirps_array_mean(:,a),[],'MarkerFaceColor',cmap1(cmap_bin,:),'MarkerEdgeColor','k','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1)
        
        iter = iter + 1;
    end
    caxis([ap_axis(ap_indices(1)),ap_axis(ap_indices(end))])
    h = colorbar;

    xlabel('[Knirps] (au)');
    ylabel('fraction of {\it eve} loci still on');
    ylabel(h,'AP position')

    %grid on
    set(gca,'FontSize',14)
    %set(gca,'Color',[228,221,209]/255) 
    xlim([2 12.5])
    ylim([-0.05 1])
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.XAxis(1).Color = 'k';

    ap_knirps_plot_fig.InvertHardcopy = 'off';
    set(gcf,'color','w'); 
    
    pbaspect([3 2 1])

    %saveas(ap_knirps_fig,[FigurePath 'figure2_fraction_on_vs_knirps_fit_' num2str(i) '.png'])
    %saveas(ap_knirps_fig,[FigurePath 'figure2_fraction_on_vs_knirps_fit_' num2str(i) '.pdf'])
end


%% Figure: Combine simple hill function fit with input-output
% Step 1: fit simple hill function to each ap position

nBoots = 200;
% fit curves to each ap position and compare the parameters
hill_coefficients = NaN(nBoots,length(ap_indices));
HM_points = NaN(nBoots,length(ap_indices));
options = optimoptions('lsqnonlin','Display','off');

wb = waitbar(0,'Conducting input-output fits...');
iter = 1;
for a = ap_indices
    waitbar(iter/(length(ap_indices)),wb);
    ap_window_filter = ap_groups==a & post_turn_on_flags & ~isnan(knirps_vec_long);
    still_on_flags_knirps = ~post_turn_off_flags(ap_window_filter);
    knirps_vec = double(knirps_vec_long(ap_window_filter))*1e-5;
    
    for n = 1:nBoots
        boot_indices = randsample(1:length(knirps_vec),length(knirps_vec),true);
        % define the objective function
        objective_fun = @(x) x(1)^x(2) ./ (x(1)^x(2) + knirps_vec(boot_indices).^x(2)) - still_on_flags_knirps(boot_indices);
        x = lsqnonlin(objective_fun,[1 1],[0 0],[Inf Inf],options);    
        HM_points(n,iter) = x(1);
        hill_coefficients(n,iter) = x(2);       
    end
    iter = iter + 1;
end
delete(wb);

%% Step2: Add fits to plots
hill_fun = @(x) x(1)^x(2) ./ (x(1)^x(2) + knirps_axis.^x(2));

ap_cell = {ap_indices(1:end-1)};

for i = length(ap_cell)
    ap_indices_iter = ap_cell{i};
    
    ap_knirps_fig = figure;
    hold on
    cmap1 = flipud(brewermap(1.2*cmap_bin_num,'Spectral'));
    %cmap1 = flipud(brewermap(1.2*cmap_bin_num,'RdYlGn'));
    cmap1 = cmap1(0.10*cmap_bin_num:1.10*cmap_bin_num-1,:);
    
    colormap(cmap1)
    iter = 1;
    for a = ap_indices_iter
        
        % add fit profile
        fit_hill = mean(hill_coefficients(:,iter));
        fit_Kd = mean(HM_points(:,iter));
        fit_profile = hill_fun([fit_Kd,fit_hill]);
        
        cmap_bin = round((ap_axis(a)-ap_axis(ap_indices(1)))*(cmap_bin_num-1)/(ap_axis(ap_indices(end))-ap_axis(ap_indices(1))))+1;
        plot(knirps_axis,fit_profile,'Color',brighten(cmap1(cmap_bin,:),-0.5),'LineWidth',0.75)
        %plot(knirps_axis,fit_profile,'Color','k','LineWidth',1)
        errorbar(knirps_axis,frac_on_knirps_array_mean(:,a),frac_on_knirps_array_ste(:,a),'o','Color',[0 0 0 0],'CapSize',0);
    %     plot(knirps_axis,frac_on_knirps_array_mean(:,a),'color',[0 0 0 0.2]);
        scatter(knirps_axis,frac_on_knirps_array_mean(:,a),[],'MarkerFaceColor',cmap1(cmap_bin,:),'MarkerEdgeColor','k','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1)
        
        iter = iter + 1;
    end
    caxis([ap_axis(ap_indices(1)),ap_axis(ap_indices(end))])
    h = colorbar;

    xlabel('[Knirps] (au)');
    ylabel('fraction of {\it eve} loci still on');
    ylabel(h,'AP position')

    %grid on
    set(gca,'FontSize',14)
    %set(gca,'Color',[228,221,209]/255) 
    xlim([2 12.5])
    ylim([-0.05 1])
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.XAxis(1).Color = 'k';

    ap_knirps_fig.InvertHardcopy = 'off';
    set(gcf,'color','w'); 
    
    pbaspect([3 2 1])

    saveas(ap_knirps_fig,[FigurePath 'figure2_fraction_on_vs_knirps_fit_' num2str(i) '.png'])
    saveas(ap_knirps_fig,[FigurePath 'figure2_fraction_on_vs_knirps_fit_' num2str(i) '.pdf'])
end

%% Figure 2: Hill Coeff vs AP 

ap_hill_fig = figure;
hold on
%cmap1 = brewermap([],'Set2');

cmap1 = flipud(brewermap(1.2*cmap_bin_num,'Spectral'));
%cmap1 = flipud(brewermap(1.2*cmap_bin_num,'RdYlGn'));
cmap1 = cmap1(0.10*cmap_bin_num:1.10*cmap_bin_num-1,:);

errorbar(ap_axis(ap_filter),nanmean(hill_coefficients),nanstd(hill_coefficients),'Color','k','CapSize',0);
iter = 1;
for i = ap_cell{1}
    cmap_bin = round((ap_axis(i)-ap_axis(ap_indices(1)))*(cmap_bin_num-1)/(ap_axis(ap_indices(end))-ap_axis(ap_indices(1))))+1;
    scatter(ap_axis(i),nanmean(hill_coefficients(:,iter)),75,'MarkerFaceColor',cmap1(cmap_bin,:),'MarkerEdgeColor','k')
    iter = iter + 1;
end
xlabel('AP position (% embryo length)');
ylabel('response sharpness (hill coefficient)');

%grid on
set(gca,'FontSize',14)
%set(gca,'Color',[228,221,209]/255) 
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

xlim([ap_axis(ap_indices(1))-0.5 ap_axis(ap_indices(end))+0.5])
ylim([0 12])
  
ap_hill_fig.InvertHardcopy = 'off';
set(gcf,'color','w'); 

pbaspect([3 2 1])

saveas(ap_hill_fig,[FigurePath 'figure2_ap_vs_hill.png'])
saveas(ap_hill_fig,[FigurePath 'figure2_ap_vs_hill.pdf'])
%% Figure 2: Kd vs AP

ap_HM_fig = figure;
hold on
%cmap1 = brewermap([],'Set2');

cmap1 = flipud(brewermap(1.2*cmap_bin_num,'Spectral'));
%cmap1 = flipud(brewermap(1.2*cmap_bin_num,'RdYlGn'));
cmap1 = cmap1(0.10*cmap_bin_num:1.10*cmap_bin_num-1,:);
%cmap1 = brighten(cmap1,-0.5);

errorbar(ap_axis(ap_filter),nanmean(HM_points),nanstd(HM_points),'Color','k','CapSize',0);
iter = 1;
for i = ap_cell{1}
    cmap_bin = round((ap_axis(i)-ap_axis(ap_indices(1)))*(cmap_bin_num-1)/(ap_axis(ap_indices(end))-ap_axis(ap_indices(1))))+1;
    scatter(ap_axis(i),nanmean(HM_points(:,iter)),75,'MarkerFaceColor',cmap1(cmap_bin,:),'MarkerEdgeColor','k')
    iter = iter + 1;
end

xlabel('AP position (% embryo length)');
ylabel('half-max points (K_d)');

%grid on
set(gca,'FontSize',14)
%set(gca,'Color',[228,221,209]/255) 
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
xlim([ap_axis(ap_indices(1))-0.5 ap_axis(ap_indices(end))+0.5])
ap_HM_fig.InvertHardcopy = 'off';
set(gcf,'color','w'); 

pbaspect([3 2 1])

saveas(ap_HM_fig,[FigurePath 'figure2_ap_vs_HM.png'])
saveas(ap_HM_fig,[FigurePath 'figure2_ap_vs_HM.pdf'])


%% Figure 2: plot fraction_on vs knirps + hill equation fit (only for stripe center)
% Also plot Knirps concentration when silencing happened

%Look at long vectors and calculate fraction on for stripe center
ap_bins_center = [-0.02 0.02];
%ap_bins_center = [-0.03 0.03];

% Step 1: Calculate using long vectors

knirps_vec_long = knirps_vec_long_raw - eYFP_background;

nLinBins = 26;
knirps_bins = linspace(0,17.5e5,nLinBins);
knirps_axis = knirps_bins(1:end-1) + diff(knirps_bins)/2;
knirps_axis = knirps_axis*1e-5;

hill_fun = @(x) x(1)^x(2) ./ (x(1)^x(2) + knirps_axis.^x(2));

ap_groups = discretize(ap_vec_long,ap_bins_center);  
knirps_groups = discretize(knirps_vec_long,knirps_bins); 

frac_on_knirps_array_mean = NaN(length(knirps_bins)-1,length(ap_bins)-1);
frac_on_knirps_array_ste = NaN(length(knirps_bins)-1,length(ap_bins)-1);

% Step 2: fit simple hill function to each ap position
ap_indices = length(ap_bins_center)-1;

nBoots = 100;
% fit curves to each ap position and compare the parameters
hill_coefficients = NaN(nBoots,length(ap_indices));
HM_points = NaN(nBoots,length(ap_indices));
options = optimoptions('lsqnonlin','Display','off');

wb = waitbar(0,'Conducting input-output fits...');
iter = 1;
for a = ap_indices
    waitbar(iter/(length(ap_indices)),wb);
    ap_window_filter = ap_groups==a & post_turn_on_flags & ~isnan(knirps_vec_long);
    still_on_flags_knirps = ~post_turn_off_flags(ap_window_filter);
    knirps_vec = double(knirps_vec_long(ap_window_filter))*1e-5;
    
    for n = 1:nBoots
        boot_indices = randsample(1:length(knirps_vec),length(knirps_vec),true);
        % define the objective function
        objective_fun = @(x) x(1)^x(2) ./ (x(1)^x(2) + knirps_vec(boot_indices).^x(2)) - still_on_flags_knirps(boot_indices);
        x = lsqnonlin(objective_fun,[1 1],[0 0],[Inf Inf],options);    
        HM_points(n,iter) = x(1);
        hill_coefficients(n,iter) = x(2);       
    end
    iter = iter + 1;
end
delete(wb);


%Step 3: Plot the results

for t = 1:length(knirps_bins)
    for a = 1:length(ap_bins_center)
                
        knirps_window_filter = knirps_groups==t & ap_groups==a & post_turn_on_flags;
        if sum(knirps_window_filter)>10
            still_on_flags_knirps = ~post_turn_off_flags(knirps_window_filter);
            boot_samples_knirps = bootstrp(nBoots,@mean,still_on_flags_knirps);
            frac_on_knirps_array_mean(t,a) = mean(boot_samples_knirps);
            frac_on_knirps_array_ste(t,a) = std(boot_samples_knirps);
        end
    end
end

% Step 3: plot the results

tiled_fig = figure('Position', [200 200 425 600]);
t = tiledlayout(2,1);

ax1 = nexttile(t);
hold on
for a = length(ap_bins_center)-1
    
    % add fit profile
    fit_hill = mean(hill_coefficients(:,a));
    fit_Kd = mean(HM_points(:,a));
    fit_profile = hill_fun([fit_Kd,fit_hill]);
    
    plot(knirps_axis,fit_profile,'Color',brighten(color_green,-0.5),'LineWidth',1)
    errorbar(knirps_axis,frac_on_knirps_array_mean(:,a),frac_on_knirps_array_ste(:,a),'o','Color',[0 0 0 0],'CapSize',0,'LineWidth',.5);
    %plot(knirps_axis,frac_on_knirps_array_mean(:,a),'color','k');

    scatter(knirps_axis,frac_on_knirps_array_mean(:,a),50,'MarkerFaceColor',color_green,'MarkerEdgeColor','k')
end

ylabel('fraction of cells transcribing');

%grid on
set(gca,'FontSize',14)
%set(gca,'Color',[228,221,209]/255) 
%xlim([2 11])
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

set(gcf,'color','w'); 

ylim([-0.05 1])
xlabel('[Knirps] (au)');

ax2 = nexttile(t);
set(gca,'FontSize',14)
h = histogram((off_knirps_long-eYFP_background)*1e-5,16,'FaceColor',mRNA_red,'Normalization','probability');

% fit with gaussian
edges = h.BinEdges;

x = (edges(1:end-1)+edges(2:end))/2;
y = h.Values;

f = fit(x',y','gauss1');

hold on
plot(f)

ylabel('fraction of silencing events')


linkaxes([ax1,ax2],'x');
xlabel(t,'[Knirps] (au)');
t.TileSpacing = 'compact';
xlim([-1.5 16.5])

%pbaspect([3 2 1])

%saveas(tiled_fig,[FigurePath 'figure2_fraction_on_events_vs_knirps_center.png'])
%saveas(tiled_fig,[FigurePath 'figure2_fraction_on_events_vs_knirps_center.pdf'])

%% Figure: plot off time distribution

f = figure;
set(gca,'FontSize',14)
h = histogram(off_time_long/60,16,'FaceColor',mRNA_red,'Normalization','probability');

% fit with gaussian
%edges = h.BinEdges;

%x = (edges(1:end-1)+edges(2:end))/2;
%y = h.Values;

%f = fit(x',y','gauss1');

%hold on
%plot(f)

xlim([0 40])
ylabel('fraction of silencing events')


linkaxes([ax1,ax2],'x');
xlabel(t,'[Knirps] (au)');
t.TileSpacing = 'compact';

%pbaspect([3 2 1])

%saveas(tiled_fig,[FigurePath 'figure2_fraction_on_events_vs_knirps_center.png'])
%saveas(tiled_fig,[FigurePath 'figure2_fraction_on_events_vs_knirps_center.pdf'])

