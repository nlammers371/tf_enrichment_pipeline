clear
close all

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
    ap_vec = spot_struct(i).apPosNucleus;
    
    if time_vec(end) - time_vec(1) >=20*60
        
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

%% Look at mean vectors
nBoots = 100;

ap_bins = linspace(0.49,0.76,28);
ap_groups = discretize(ap_vec_long,ap_bins); 
ap_groups_mean = discretize(mean_ap,ap_bins); 
ap_groups_off = discretize(off_ap,ap_bins); 

frac_on_vec_mean = NaN(1,length(ap_bins)-1);
frac_on_vec_ste = NaN(1,length(ap_bins)-1);

off_time_vec_mean = NaN(1,length(ap_bins)-1);
off_time_vec_ste = NaN(1,length(ap_bins)-1);

fluo_vec_mean = NaN(1,length(ap_bins)-1);
fluo_vec_ste = NaN(1,length(ap_bins)-1);


for a = 1:length(ap_bins)
    ap_filter_long = ap_groups==a;
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
             
end

% Define some colors  
yw = [234 194 100]/255; % yellow
bl = [115 143 193]/255; % blue
gr = [191 213 151]/255; % green

ap_axis = 100*(ap_bins(1:end-1) + diff(ap_bins)/2);
ap_filter = ap_axis>=55 & ap_axis <=67;
ap_indices = find(ap_filter);
ap_axis = ap_axis - 61;
% fraction on
fraction_on_fig = figure;
hold on

errorbar(ap_axis,frac_on_vec_mean,frac_on_vec_ste,'Color','k','CapSize',0);
scatter(ap_axis,frac_on_vec_mean,50,'MarkerFaceColor',yw,'MarkerEdgeColor','k')

xlabel('AP position (% embryo length)');
ylabel('fraction of active nuclei');

grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
xlim([ap_axis(1) ap_axis(end)])
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

fraction_on_fig.InvertHardcopy = 'off';
set(gcf,'color','w'); 
ylim([0.2 1.1])    
saveas(fraction_on_fig,[FigurePath 'fraction_on_vs_ap.png'])
saveas(fraction_on_fig,[FigurePath 'fraction_on_vs_ap.pdf'])


% off time
off_time_fig = figure;
hold on

errorbar(ap_axis,off_time_vec_mean/60,off_time_vec_ste/60,'Color','k','CapSize',0);
scatter(ap_axis,off_time_vec_mean/60,50,'MarkerFaceColor',bl,'MarkerEdgeColor','k')

xlabel('AP position (% embryo length)');
ylabel('average off time (minutes)');

grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
ylim([10 40])
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
xlim([ap_axis(1) ap_axis(end)])
off_time_fig.InvertHardcopy = 'off';
set(gcf,'color','w'); 
    
saveas(off_time_fig,[FigurePath 'off_time_vs_ap.png'])
saveas(off_time_fig,[FigurePath 'off_time_vs_ap.pdf'])


mean_fig = figure;
hold on

errorbar(ap_axis,fluo_vec_mean*1e-5,fluo_vec_ste*1e-5,'Color','k','CapSize',0);
scatter(ap_axis,fluo_vec_mean*1e-5,50,'MarkerFaceColor',gr,'MarkerEdgeColor','k')

xlabel('AP position (% embryo length)');
ylabel('mean spot intensity (au)');

grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
ylim([.2 1.4])
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
xlim([ap_axis(1) ap_axis(end)])
mean_fig.InvertHardcopy = 'off';
set(gcf,'color','w'); 
    
saveas(mean_fig,[FigurePath 'fluo_vs_ap.png'])
saveas(mean_fig,[FigurePath 'fluo_vs_ap.pdf'])
    
%% Look at long vectors
knirps_offset = 2.5e5;%prctile(double(knirps_vec_long),1);
knirps_vec_long = knirps_vec_long_raw - knirps_offset;

nLinBins = 31;
time_bins = linspace(5*60,40*60,nLinBins);
knirps_bins = linspace(0,15e5,nLinBins);

ap_groups = discretize(ap_vec_long,ap_bins); 
time_groups = discretize(time_vec_long,time_bins); 
knirps_groups = discretize(knirps_vec_long,knirps_bins); 

frac_on_time_array_mean = NaN(length(time_bins)-1,length(ap_bins)-1);
frac_on_knirps_array_mean = NaN(length(time_bins)-1,length(ap_bins)-1);
eve_time_array_mean = NaN(length(time_bins)-1,length(ap_bins)-1);
eve_time_array_full = NaN(length(time_bins)-1,length(ap_bins)-1,nBoots);
knirps_time_array_mean = NaN(length(time_bins)-1,length(ap_bins)-1);

frac_on_time_array_ste = NaN(length(time_bins)-1,length(ap_bins)-1);
frac_on_knirps_array_ste = NaN(length(time_bins)-1,length(ap_bins)-1);
eve_time_array_ste = NaN(length(time_bins)-1,length(ap_bins)-1);
knirps_time_array_ste = NaN(length(time_bins)-1,length(ap_bins)-1);


for t = 1:length(time_bins)
    for a = 1:length(ap_bins)
        time_window_filter = time_groups==t & ap_groups==a & post_turn_on_flags;
        if sum(time_window_filter) > 10
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
                
        knirps_window_filter = knirps_groups==t & ap_groups==a & post_turn_on_flags;
        if sum(knirps_window_filter)>10
            still_on_flags_knirps = ~post_turn_off_flags(knirps_window_filter);
            boot_samples_knirps = bootstrp(nBoots,@mean,still_on_flags_knirps);
            frac_on_knirps_array_mean(t,a) = mean(boot_samples_knirps);
            frac_on_knirps_array_ste(t,a) = std(boot_samples_knirps);
        end
    end
end
  
%% Generate predicted mRNA
% generate decay kernel
close all

eve_half_life = 7;
eve_decay_kernel = exp(-time_bins'/eve_half_life/60);
eve_decay_kernel = eve_decay_kernel / eve_decay_kernel(1);

% replace missing values with zeros (Need to clean this up, should really
% be some kind of interpolation)
eve_time_array_full(isnan(eve_time_array_full)) = 0;

predicted_eve_profile_array = convn(eve_decay_kernel,eve_time_array_full,'full');
predicted_eve_profile_mean = nanmean(predicted_eve_profile_array,3);
predicted_eve_profile_ste = nanstd(predicted_eve_profile_array,[],3);

predicted_eve_profile_mean = predicted_eve_profile_mean(1:length(time_bins),:);
predicted_eve_profile_ste = predicted_eve_profile_ste(1:length(time_bins),:);

% knirps green
k_green = brighten([38 142 75]/256,.4);
mRNA_red = brighten([212 100 39]/256,.2);

% mRNA profile plot
mRNA_fig = figure;
hold on

yyaxis left
% f = fill([ap_axis fliplr(ap_axis)], [knirps_time_array_mean(end,:)*1e-5 zeros(size(knirps_time_array_mean(end,:)))],k_green);
plot(ap_axis,knirps_time_array_mean(end,:)*1e-5,'Color','k','LineWidth',3)
% f.FaceAlpha = 0.2;
ylabel('[Knirps] (au)');
set(gca,'YColor',k_green)

yyaxis right
% errorbar(ap_axis,imgaussfilt(predicted_eve_profile_mean(end,:)*1e-5,1),predicted_eve_profile_ste(end,:)*1e-5,'Color',mRNA_red,'CapSize',0);
% plot(ap_axis,imgaussfilt(predicted_eve_profile_mean(end,:)*1e-5,1),'Color',brighten(mRNA_red,0),'LineWidth',1.5) % NL: applying mild smoothing since this is illustrative (not quantitative)
% s = scatter(ap_axis,imgaussfilt(predicted_eve_profile_mean(end,:)*1e-5,1),50,'MarkerFaceColor',mRNA_red,'MarkerEdgeColor','k');
eve_mRNA_sm = imgaussfilt(predicted_eve_profile_mean(end,:)*1e-5,1);
f = fill([ap_axis fliplr(ap_axis)], [eve_mRNA_sm zeros(size(eve_mRNA_sm))],mRNA_red);
f.FaceAlpha = 0.5;
set(gca,'YColor',mRNA_red)
xlabel('AP position (% embryo length)');
ylabel('accumulated {\it eve} mRNA (au)');


grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
% ylim([0.9 1.1])

xlim([ap_axis(1) ap_axis(end)])
mRNA_fig.InvertHardcopy = 'off';
set(gcf,'color','w'); 

saveas(mRNA_fig,[FigurePath 'mRNA_fig.png'])
saveas(mRNA_fig,[FigurePath 'mRNA_fig.pdf'])

%% make plots
ap_cell = {ap_indices(1:3) ap_indices(1:6) ap_indices(1:9) ap_indices(1:end)};
close all
knirps_axis = knirps_bins(1:end-1) + diff(knirps_bins)/2;
knirps_axis = knirps_axis*1e-5;
time_axis = time_bins(1:end-1) + diff(time_bins)/2;
time_axis = time_axis/60;


for i = 1:length(ap_cell)
    ap_indices_iter = ap_cell{i};
  
    ap_time_fig = figure;
    hold on
    cmap1 = brewermap(length(ap_indices),'blues');
    colormap(cmap1)
    iter = 1;
    for a = ap_indices_iter
        errorbar(time_axis,frac_on_time_array_mean(:,a),frac_on_time_array_ste(:,a),'Color','k','CapSize',0);
        scatter(time_axis,frac_on_time_array_mean(:,a),'MarkerFaceColor',cmap1(iter,:),'MarkerEdgeColor','k')
        iter = iter + 1;
    end
    caxis([ap_axis(ap_indices(1)) ap_axis(ap_indices(end))])
    h = colorbar;
    xlim([time_bins(1)/60 time_bins(end)/60])
    xlabel('time since start of nc14 (minutes)');
    ylabel('fraction of {\it eve} loci still on');
    ylabel(h,'AP position')

    grid on
    set(gca,'FontSize',14)
    set(gca,'Color',[228,221,209]/255) 

    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.XAxis(1).Color = 'k';

    ap_time_fig.InvertHardcopy = 'off';
    set(gcf,'color','w'); 

    saveas(ap_time_fig,[FigurePath 'fraction_on_vs_time_' num2str(i) '.png'])
    saveas(ap_time_fig,[FigurePath 'fraction_on_vs_time_' num2str(i) '.pdf'])
end
%%

for i = 1:length(ap_cell)
    ap_indices_iter = ap_cell{i};
    
    
    ap_knirps_fig = figure;
    hold on
    cmap1 = flipud(brewermap(length(ap_indices),'PrGn'));
    colormap(cmap1)
    iter = 1;
    for a = ap_indices_iter
        errorbar(knirps_axis,frac_on_knirps_array_mean(:,a),frac_on_time_array_ste(:,a),'o','Color',[0 0 0 0],'CapSize',0);
    %     plot(knirps_axis,frac_on_knirps_array_mean(:,a),'color',[0 0 0 0.2]);
        scatter(knirps_axis,frac_on_knirps_array_mean(:,a),[],'MarkerFaceColor',cmap1(iter,:),'MarkerEdgeColor','k')
        iter = iter + 1;
    end
    caxis([ap_axis(ap_indices(1)),ap_axis(ap_indices(end))])
    h = colorbar;

    xlabel('[Knirps] (au)');
    ylabel('fraction of {\it eve} loci still on');
    ylabel(h,'AP position')

    grid on
    set(gca,'FontSize',14)
    set(gca,'Color',[228,221,209]/255) 
    xlim([1 11])
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.XAxis(1).Color = 'k';

    ap_knirps_fig.InvertHardcopy = 'off';
    set(gcf,'color','w'); 

    saveas(ap_knirps_fig,[FigurePath 'fraction_on_vs_knirps_' num2str(i) '.png'])
    saveas(ap_knirps_fig,[FigurePath 'fraction_on_vs_knirps_' num2str(i) '.pdf'])
end

% %% Make a second version of the plots
% close all
% ap_knirps_fig = figure;
% hold on
% cmap1 = flipud(brewermap(length(ap_indices),'spectral'));
% colormap(cmap1)
% iter = 1;
% for a = ap_indices
%     
%     ub = frac_on_knirps_array_mean(:,a) + 3*frac_on_knirps_array_ste(:,a);
%     lb = frac_on_knirps_array_mean(:,a) - 3*frac_on_knirps_array_ste(:,a);
%     nan_filter = isnan(ub);
%     fill([knirps_axis(~nan_filter) fliplr(knirps_axis(~nan_filter))],[ub(~nan_filter)' ...
%                     fliplr(lb(~nan_filter)')],cmap1(iter,:),'FaceAlpha',.5,'EdgeAlpha',0.2)
%     
% %     fill([time_axis fliplr(time_axis)],[master_struct(1).br_spot_ub fliplr(master_struct(1).br_spot_lb)],cmap1(2,:)
% %     scatter(knirps_axis,frac_on_knirps_array_mean(:,a),50,'MarkerFaceColor',cmap1(iter,:),'MarkerEdgeColor','k')
%     iter = iter + 1;
% end
% 
% caxis(100*[ap_bins(ap_indices(1)),ap_bins(ap_indices(end))])
% h = colorbar;
% 
% xlabel('Knirps concentration (au)');
% ylabel('fraction of {\it eve} loci still on');
% ylabel(h,'AP position')
% 
% grid on
% set(gca,'FontSize',14)
% set(gca,'Color',[228,221,209]/255) 
% xlim([2 12])
% ax = gca;
% ax.YAxis(1).Color = 'k';
% ax.XAxis(1).Color = 'k';
% 
% ap_knirps_fig.InvertHardcopy = 'off';
% set(gcf,'color','w'); 

% %%
% caxis(100*[ap_bins(ap_indices(1)),ap_bins(ap_indices(end))])
% h = colorbar;
% 
% xlabel('[Knirps] (au)');
% ylabel('fraction of {\it eve} loci still on');
% ylabel(h,'AP position')
% 
% grid on
% set(gca,'FontSize',14)
% set(gca,'Color',[228,221,209]/255) 
% xlim([0 10])
% ax = gca;
% ax.YAxis(1).Color = 'k';
% ax.XAxis(1).Color = 'k';
% 
% ap_knirps_fig.InvertHardcopy = 'off';
% set(gcf,'color','w'); 
%     
% saveas(ap_knirps_fig,[FigurePath 'fraction_on_vs_knirps.png'])
% saveas(ap_knirps_fig,[FigurePath 'fraction_on_vs_knirps.pdf'])


%% fit simple hill function to each ap position

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

%% Add fits to plots
hill_fun = @(x) x(1)^x(2) ./ (x(1)^x(2) + knirps_axis.^x(2));


for i = length(ap_cell)
    ap_indices_iter = ap_cell{i};
    
    
    ap_knirps_fig = figure;
    hold on
    cmap1 = flipud(brewermap(length(ap_indices),'PrGn'));
    colormap(cmap1)
    iter = 1;
    for a = ap_indices_iter
        
        % add fit profile
        fit_hill = mean(hill_coefficients(:,iter));
        fit_Kd = mean(HM_points(:,iter));
        fit_profile = hill_fun([fit_Kd,fit_hill]);
        
        plot(knirps_axis,fit_profile,'Color',[cmap1(iter,:) .5],'LineWidth',1.5)
        
        
        errorbar(knirps_axis,frac_on_knirps_array_mean(:,a),frac_on_time_array_ste(:,a),'o','Color',[0 0 0 0],'CapSize',0);
    %     plot(knirps_axis,frac_on_knirps_array_mean(:,a),'color',[0 0 0 0.2]);
        scatter(knirps_axis,frac_on_knirps_array_mean(:,a),[],'MarkerFaceColor',cmap1(iter,:),'MarkerEdgeColor','k','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1)
        
        
        iter = iter + 1;
    end
    caxis([ap_axis(ap_indices(1)),ap_axis(ap_indices(end))])
    h = colorbar;

    xlabel('[Knirps] (au)');
    ylabel('fraction of {\it eve} loci still on');
    ylabel(h,'AP position')

    grid on
    set(gca,'FontSize',14)
    set(gca,'Color',[228,221,209]/255) 
    xlim([1 11])
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.XAxis(1).Color = 'k';

    ap_knirps_fig.InvertHardcopy = 'off';
    set(gcf,'color','w'); 

    saveas(ap_knirps_fig,[FigurePath 'fraction_on_vs_knirps_fit_' num2str(i) '.png'])
    saveas(ap_knirps_fig,[FigurePath 'fraction_on_vs_knirps_fit_' num2str(i) '.pdf'])
end

%%
close all

ap_hill_fig = figure;
hold on
cmap1 = brewermap([],'Set2');

errorbar(ap_axis(ap_filter),nanmean(hill_coefficients),nanstd(hill_coefficients),'Color','k','CapSize',0);
scatter(ap_axis(ap_filter),nanmean(hill_coefficients),75,'MarkerFaceColor',cmap1(3,:),'MarkerEdgeColor','k')

xlabel('AP position (% embryo length)');
ylabel('response sharpness (hill coefficient)');

grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

xlim([ap_axis(ap_indices(1))-0.5 ap_axis(ap_indices(end))+0.5])
ylim([0 12])
  
ap_hill_fig.InvertHardcopy = 'off';
set(gcf,'color','w'); 
    
saveas(ap_hill_fig,[FigurePath 'ap_vs_hill.png'])
saveas(ap_hill_fig,[FigurePath 'ap_vs_hill.pdf'])
%%


ap_HM_fig = figure;
hold on
cmap1 = brewermap([],'Set2');

errorbar(ap_axis(ap_filter),nanmean(HM_points),nanstd(HM_points),'Color','k','CapSize',0);
scatter(ap_axis(ap_filter),nanmean(HM_points),75,'MarkerFaceColor',cmap1(5,:),'MarkerEdgeColor','k')

xlabel('AP position (% embryo length)');
ylabel('half-max points (K_d)');

grid on
set(gca,'FontSize',14)
set(gca,'Color',[228,221,209]/255) 
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
xlim([ap_axis(ap_indices(1))-0.5 ap_axis(ap_indices(end))+0.5])
ap_HM_fig.InvertHardcopy = 'off';
set(gcf,'color','w'); 
    
saveas(ap_HM_fig,[FigurePath 'ap_vs_HM.png'])
saveas(ap_HM_fig,[FigurePath 'ap_vs_Hm.pdf'])