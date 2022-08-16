clear
close all

addpath(genpath('./lib'))

eYFP_background = 375698.13;

% knirps green
k_green = brighten([38 142 75]/256,.4);
color_green = [38 143 75]/256; % color from Jake
mRNA_red = brighten([212 100 39]/256,.2);

%% Initialization

expID = 5;

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
    
    if (spot_struct(i).setID == expID)
  
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
end


% from original one
knirps_vec_long = knirps_vec_long_raw - eYFP_background;

nBoots = 100;

nLinBins = 41;
ap_bins = linspace(-0.105,0.105,25); % used for input-output function
%ap_bins = linspace(-0.105,0.105,15); % used in mean figure
time_bins = linspace(5*60,40*60,nLinBins);
%ap_bins = linspace(-0.125,0.125,26);
%ap_bins = linspace(-0.12,0.12,21);

% Define some colors  
yw = [234 194 100]/255; % yellow
bl = [115 143 193]/255; % blue
gr = [191 213 151]/255; % green

ap_axis = 100*(ap_bins(1:end-1) + diff(ap_bins)/2);


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
    
% add fit profile
fit_hill = mean(hill_coefficients(:,1));
fit_Kd = mean(HM_points(:,1));
fit_profile = hill_fun([fit_Kd,fit_hill]);

plot(knirps_axis,fit_profile,'Color',brighten(color_green,-0.5),'LineWidth',1)
errorbar(knirps_axis,frac_on_knirps_array_mean(:,1),frac_on_knirps_array_ste(:,1),'o','Color',[0 0 0 0],'CapSize',0,'LineWidth',.5);
%plot(knirps_axis,frac_on_knirps_array_mean(:,a),'color','k');

scatter(knirps_axis,frac_on_knirps_array_mean(:,1),50,'MarkerFaceColor',color_green,'MarkerEdgeColor','k')

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