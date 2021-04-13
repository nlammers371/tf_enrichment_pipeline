% script to generate stripe dynamics figures
clear
close all

baseProject = 'NSv1';
liveProject = LiveEnrichmentProject(baseProject);
% make output paths for data and figures
slashesFig = strfind(liveProject.figurePath,'\');
FigurePath = [liveProject.figurePath(1:slashesFig(end-1)) 'eve2_recon_analyses' filesep];
mkdir(FigurePath)

slashesData = strfind(liveProject.dataPath,'\');
DataPath = [liveProject.dataPath(1:slashesData(end-1)) 'eve2_recon_analyses' filesep];

%% %%%%%%%%%%%%%%%%%%%% Make offset plot and HM plots %%%%%%%%%%%%%%%%%%%%

% load hm dataset
load([DataPath 'hm_info_struct.mat'],'hm_info_struct')
% read all fields into the workspace
fnames = fieldnames(hm_info_struct);
for f = 1:length(fnames)
    eval([fnames{f} ' = hm_info_struct.(fnames{f});'])
end

offset_mean_array = permute(nanmean(offset_dynamics_array,2),[1 3 2]);
time_axis = time_index(1:end-1) + diff(time_index)/2;

offset_fig = figure;
hold on
cmap1 = brewermap([],'Set2');

pp = [];
for p = 1:length(project_id_vec)
    pp(end+1) = plot(time_axis,offset_mean_array(:,p),'Color',cmap1(project_id_vec(p),:),'LineWidth',1.5);  
end

xlabel('time into nc14 (minutes)')
ylabel('fluorescence offset (au)')
legend(pp(legend_ids),legend_str{:})
set(gca,'Fontsize',14);
grid on
xlim([0 50])
saveas(offset_fig,[FigurePath 'fluo_offset.png'])

% %%%%%%%%%%%%%%%%%%% Make mRNA plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
mRNA_path = [FigurePath 'mRNA_heatmaps' filesep];
mkdir(mRNA_path)
max_mRNA = ceil(prctile(full_dynamics_array(:),98)/1e4)*1e4;


for p = 1:length(project_id_vec)

    heat_fig = figure('Visible','off');
    cmap = brewermap([],'YlOrRd');
    colormap(cmap);
    % make plot
    pc = pcolor(flipud(full_dynamics_array(:,:,p)));
    pc.EdgeAlpha = 0.1; 

    % set axes
    set(gca,'ytick',1:10:max_time,'yticklabel',fliplr(time_index(1:10:max_time)))
    set(gca,'xtick',1:5:length(ap_index),'xticklabel',ap_index(1:5:length(ap_index)))    
    set(gca,'Fontsize',14);

    % axes
    h = colorbar;
    ylabel('time into nc14 (minutes)')
    xlabel('% embryo length')
    ylabel(h,'average mRNA production (AU)')
    xlim([6 31])
    % save
    saveas(heat_fig,[mRNA_path 'scaled_mRNA_hm_' projectNameCell{project_id_vec(p)} '_embryo' num2str(embryo_id_vec(p)) '.png'])
    
    caxis([0 max_mRNA])
    
    saveas(heat_fig,[mRNA_path 'standardized_mRNA_hm_' projectNameCell{project_id_vec(p)} '_embryo' num2str(embryo_id_vec(p)) '_standardized.png'])
end


% %%%%%%%%%%%%%%%%%%% Make frac on plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
frac_path = [FigurePath 'frac_heatmaps' filesep];
mkdir(frac_path)
max_frac = 0.7;


for p = 1:length(project_id_vec)

    heat_fig = figure('Visible','off');
    cmap = brewermap([],'YlGnBu');
    colormap(cmap);
    % make plot
    pc = pcolor(flipud(instant_frac_dynamics_array(:,:,p)));
    pc.EdgeAlpha = 0.1; 

    % set axes
    set(gca,'ytick',1:10:max_time,'yticklabel',fliplr(time_index(1:10:max_time)))
    set(gca,'xtick',1:5:length(ap_index),'xticklabel',ap_index(1:5:length(ap_index)))    
    set(gca,'Fontsize',14);

    % axes
    h = colorbar;
    ylabel('time into nc14 (minutes)')
    xlabel('% embryo length')
    ylabel(h,'instantaneous fraction active')
    xlim([6 31])
    % save
%     saveas(heat_fig,[mRNA_path 'scaled_mRNA_hm_' projectNameCell{project_id_vec(p)} '_embryo' num2str(embryo_id_vec(p)) '.png'])
    
    caxis([0 max_frac])
    
    saveas(heat_fig,[frac_path 'standardized_frac_hm_' projectNameCell{project_id_vec(p)} '_embryo' num2str(embryo_id_vec(p)) '_standardized.png'])
end


% %%%%%%%%%%%%%%%%%%% Make fluo on plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
fluo_path = [FigurePath 'fluo_heatmaps' filesep];
mkdir(fluo_path)
max_fluo = ceil(prctile(fluo_dynamics_array(:),95)/1e4)*1e4;


for p = 1:length(project_id_vec)

    heat_fig = figure('Visible','off');
    cmap = brewermap([],'PuBu');
    colormap(cmap);
    % make plot
    pc = pcolor(flipud(fluo_dynamics_array(:,:,p)));
    pc.EdgeAlpha = 0.1; 

    % set axes
    set(gca,'ytick',1:10:max_time,'yticklabel',fliplr(time_index(1:10:max_time)))
    set(gca,'xtick',1:5:length(ap_index),'xticklabel',ap_index(1:5:length(ap_index)))    
    set(gca,'Fontsize',14);

    % axes
    h = colorbar;
    ylabel('time into nc14 (minutes)')
    xlabel('% embryo length')
    ylabel(h,'average (active) MS2 spot intensity (AU)')
    xlim([6 31])
    % save
    saveas(heat_fig,[fluo_path 'scaled_fluo_hm_' projectNameCell{project_id_vec(p)} '_embryo' num2str(embryo_id_vec(p)) '.png'])
    
    caxis([0 max_fluo])
    
    saveas(heat_fig,[fluo_path 'standardized_fluo_hm_' projectNameCell{project_id_vec(p)} '_embryo' num2str(embryo_id_vec(p)) '_standardized.png'])
end

%% Make phase space scatters
% load hm dataset
load([DataPath 'stripe_param_struct.mat'],'stripe_param_struct')
% read all fields into the workspace
fnames = fieldnames(stripe_param_struct);
for f = 1:length(fnames)
    eval([fnames{f} ' = stripe_param_struct.(fnames{f});'])
end


close all
markerAlpha = 0.5;

cw_fig = figure;
hold on
cmap1 = brewermap([],'Set2');
ss_temp = [];
for p = 1:length(project_id_vec)
    plot(fit_param_raw_array(:,2,p),fit_param_raw_array(:,3,p),'Color',cmap1(project_id_vec(p),:))
    ss_temp(end+1) = scatter(fit_param_raw_array(:,2,p),fit_param_raw_array(:,3,p),20+2*(1:size(fit_param_array,1)),...
                'MarkerFaceColor',cmap1(project_id_vec(p),:),'MarkerEdgeColor','k',...
                'MarkerFaceAlpha',markerAlpha,'MarkerEdgeAlpha',0.5);
end    
xlabel('stripe center (% embryo length)')
ylabel('stripe width (% embryo length)')
set(gca,'Fontsize',14);
grid on
saveas(cw_fig,[FigurePath 'center_v_amp.png'])

ca_fig = figure;
hold on
ss_temp = [];
for p = 1:length(project_id_vec)
    plot(fit_param_raw_array(:,2,p),fit_param_raw_array(:,1,p),'Color',cmap1(project_id_vec(p),:))
    ss_temp(end+1) = scatter(fit_param_raw_array(:,2,p),fit_param_raw_array(:,1,p),20+2*(1:size(fit_param_array,1)),...
                'MarkerFaceColor',cmap1(project_id_vec(p),:),'MarkerEdgeColor','k',...
                'MarkerFaceAlpha',markerAlpha,'MarkerEdgeAlpha',0.5);
end    
grid on
xlabel('stripe center (% embryo length)')
ylabel('stripe amplitude (au)')
set(gca,'Fontsize',14);
legend(ss_temp(legend_ids),legend_str{:})
saveas(ca_fig,[FigurePath 'center_v_amp.png'])


wa_fig = figure;
hold on
ss_temp = [];
for p = 1:length(project_id_vec)
    plot(fit_param_raw_array(:,3,p),fit_param_raw_array(:,1,p),'Color',cmap1(project_id_vec(p),:))
    ss_temp(end+1) = scatter(fit_param_raw_array(:,3,p),fit_param_raw_array(:,1,p),20+2*(1:size(fit_param_array,1)),...
                'MarkerFaceColor',cmap1(project_id_vec(p),:),'MarkerEdgeColor','k',...
                'MarkerFaceAlpha',markerAlpha,'MarkerEdgeAlpha',0.5);
end

grid on
xlabel('stripe width (% embryo length)')
ylabel('stripe amplitude (au)')
set(gca,'Fontsize',14);
legend(ss_temp(legend_ids),legend_str{:})
saveas(wa_fig,[FigurePath 'width_v_amp.png'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Temporal Dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

plot_indices = 1:3;
param_name_cell = {'stripe_amplitude','stripe_center','stripe_width'};
param_label_cell = {'stripe amplitude (AU)','stripe center (% embryo length)','stripe width (% embryo length)'};
time_to_plot = time_axis(fit_indices);
loc_str = {'northwest','southwest','northwest'};

for param_index = plot_indices
    param_name = param_name_cell{param_index};
    y_label = param_label_cell{param_index}; 
    
    param_fig = figure;
    hold on
    cmap1 = brewermap([],'Set2');

    p = [];
    for m = 1:length(projectNameCell)
        ub = mean_param_array(:,param_index,m)' + std_param_array(:,param_index,m)';
        lb = mean_param_array(:,param_index,m)' - std_param_array(:,param_index,m)';
        nf = ~isnan(ub);
        % area plot for feasible range
        ap = fill([time_to_plot(nf) fliplr(time_to_plot(nf))], [ub(nf) fliplr(lb(nf))], cmap1(m,:),'FaceAlpha',0.39,'EdgeAlpha',0.3);    
        p(m) = plot(time_to_plot,mean_param_array(:,param_index,m)','Color',cmap1(m,:),'LineWidth',1.5);

        % plot individual replicates
        p_filter = find(project_id_vec == m);
        emb_data = fit_param_raw_array(:,param_index,p_filter);
        time_long = repmat(time_to_plot',1,1,sum(p_filter));
        for em = 1:length(p_filter)
            plot(time_long(:,:,em),emb_data(:,:,em),'--','Color',[cmap1(m,:) .4],'LineWidth',1.5);
        end
    end  
    if param_index == 1
      ylim([0 16e4])
    end
    xlim([10 50])
    xlabel('minutes into nc14')
    ylabel(y_label)
    set(gca,'Fontsize',14);
    grid on
    legend(p,legend_str{:},'Location',loc_str{param_index})
        
    set(gca,'Color',[228,221,209]/255) 
          
    param_fig.InvertHardcopy = 'off';
    set(gcf,'color','w');

    saveas(param_fig,[FigurePath param_name '_v_time.png'])
end


%%%%%%%%%%%%%%%%%%%%%%%
%% Cumulative fraction ON/OFF
%%%%%%%%%%%%%%%%%%%%%%%
% load frac dataset
load([DataPath 'fraction_on_info.mat'],'fraction_on_info')
% read all fields into the workspace
fnames = fieldnames(fraction_on_info);
for f = 1:length(fnames)
    eval([fnames{f} ' = fraction_on_info.(fnames{f});'])
end
% loc_str = {'southeast','southeast','northwest'};
% make plots
close all

for i = 1:3
    frac_fig = figure;
    hold on
    cmap1 = brewermap([],'Set2');

    p = [];
    for m = 1:length(projectNameCell)
        ub = mean_frac_array(:,m,i)' + std_frac_array(:,m,i)';
        lb = mean_frac_array(:,m,i)' - std_frac_array(:,m,i)';
        nf = ~isnan(ub);

        % area plot for feasible range
        ap = fill([time_axis(nf) fliplr(time_axis(nf))], [ub(nf) fliplr(lb(nf))], cmap1(m,:),'FaceAlpha',0.15,'EdgeAlpha',0);    
        p(m) = plot(time_axis,mean_frac_array(:,m,i)','Color',cmap1(m,:),'LineWidth',2);

        % plot individual replicates
        p_filter = find(hm_info_struct.project_id_vec == m);   
        for em = 1:length(p_filter)
            plot(time_axis,frac_dynamics_array(:,p_filter(em),i),'--','Color',[cmap1(m,:) .3],'LineWidth',1.5);
        end
    end  

    xlim([0 47])    
    ylim([0 1.05])    
    xlabel('minutes into nc14')
    ylabel(param_labels{i})
    set(gca,'Fontsize',14);
    grid on
    legend(p,legend_str{:},'Location',loc_str{i})

    set(gca,'Color',[228,221,209]/255) 

    frac_fig.InvertHardcopy = 'off';
    set(gcf,'color','w');
    saveas(frac_fig,[FigurePath 'cumulative_' param_cell{i} '.png'])
end


%% plot 50% ON/OFF metrics

rng(123);
on50_fig = figure;
hold on
errorbar(1:length(projectNameCell),mean_on_50,ste_on_50,'o','Color','k','CapSize',0);
s = [];
for m = 1:length(projectNameCell)
    s(end+1) = scatter(m,mean_on_50(m),75,'s',...
              'MarkerFaceColor',cmap1(m,:),'MarkerEdgeColor','k',...
              'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);
            
    scatter(repelem(m,length(on_50{m}))+.25*rand(size(on_50{m}))-0.125,on_50{m},50,'o',...
              'MarkerFaceColor',cmap1(m,:),'MarkerEdgeColor','k',...
              'MarkerFaceAlpha',.25,'MarkerEdgeAlpha',.25);
  
end
xlim([0 5])    
ylabel('time to 50% ON (minutes)')
set(gca,'Fontsize',14);
grid on
set(gca,'xtick',1:length(projectNameCell),'xticklabels',legend_str_short)
xtickangle(-30) 
set(gca,'Color',[228,221,209]/255) 

on50_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
saveas(on50_fig,[FigurePath 'time_to_50_ON.png'])

off50_fig = figure;
hold on
errorbar(1:length(projectNameCell),mean_off_50,ste_off_50,'o','Color','k','CapSize',0);
s = [];
for m = 1:length(projectNameCell)
    s(end+1) = scatter(m,mean_off_50(m),75,'s',...
              'MarkerFaceColor',cmap1(m,:),'MarkerEdgeColor','k',...
              'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);
            
    scatter(repelem(m,length(off_50{m}))+.25*rand(size(off_50{m}))-0.125,off_50{m},50,'o',...
              'MarkerFaceColor',cmap1(m,:),'MarkerEdgeColor','k',...
              'MarkerFaceAlpha',.25,'MarkerEdgeAlpha',.25);
  
end
xlim([0 5])    
ylabel('time to 50% OFF (minutes)')
set(gca,'Fontsize',14);
grid on
set(gca,'xtick',1:length(projectNameCell),'xticklabels',legend_str_short)
xtickangle(-30) 
set(gca,'Color',[228,221,209]/255) 

on50_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
saveas(on50_fig,[FigurePath 'time_to_50_OFF.png'])