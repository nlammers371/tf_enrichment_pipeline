% Script to compile inference results 
clear
close all

addpath(genpath('./lib'))
addpath(genpath('../../utilities'))

% Load data
dataRoot = 'S:\Nick\Dropbox (Personal)\ProcessedEnrichmentData\';
FigurePath = 'S:\Nick\Dropbox (Personal)\LocalEnrichmentFigures\combinedOptoSets\';
if ~exist(dataRoot)
    dataRoot = 'C:\Users\nlamm\Dropbox (Personal)\ProcessedEnrichmentData\';
    FigurePath = 'C:\Users\nlamm\Dropbox (Personal)\LocalEnrichmentFigures\combinedOptoSets\';
end
mkdir(FigurePath)
projectName = 'combinedOptoSets_v2';
dataPath = [dataRoot projectName filesep];
load([dataPath 'inference_data.mat'],'inference_data')
load([dataPath 'inference_info.mat'],'inference_struct')

% set write path
inferencePath = [dataPath 'cpHMM_results' filesep];
mkdir(inferencePath)

% load inference results
inf_files = dir([inferencePath 'hmm_results*']);

% read raw inference files into memory
inference_results = struct;
wb = waitbar(0,'loading inference results...');
for i = 1:length(inf_files)
    waitbar(i/length(inf_files),wb);
    load([inferencePath inf_files(i).name])
    f_mean = nanmean([output.fluo_data{:}]);
    t_mean = nanmean([output.time_data{:}]);
%     kni_mean = nanmean([output.kni_data{:}]);r
    fnames = fieldnames(output);
    for f = 1:length(fnames)
        inference_results(i).(fnames{f}) = output.(fnames{f});
    end
    inference_results(i).f_mean = f_mean;
    inference_results(i).t_mean = t_mean;
%     inference_results(i).kni_mean = t_mean;
end  
delete(wb)
%%
z_index = [inference_struct.z_id];
particle_id_vec_long = inference_data.particle_id_vec;


group_id_vec = [inference_results.groupID];
knirps_id_vec = [inference_results.additionalBin];
gp_kni_array = unique([group_id_vec' knirps_id_vec'],'rows');
diff_vec = [1 diff(gp_kni_array(:,2))'];
inf_ids = cumsum(1*(diff_vec<0))+1;

% generate knirps axis vector 
% knirps_axis_vec = NaN(size(group_id_vec));
% for i = 1:length(inference_results)
%     % collect list of relevant entries (this requires a little work)
%     particle_ids = inference_results(i).particle_ids;
%     time_data = inference_results(i).time_data;
%     kni_ids = NaN(size(particle_ids));
%     for p = 1:length(particle_ids)
%         p_options = find(particle_id_vec_long==particle_ids(p));
%         for j = 1:length(p_options)
%             if inference_data.time_vec_cell{p_options(j)}(1) == time_data{p}(1)
%                 kni_ids(p) = p_options(j);
%             end
%         end
%     end        
%     knirps_axis_vec(i) = nanmean([inference_data.knirps_vec_cell{kni_ids}]);
% end  
knirps_axis_vec = [inference_results.kni_mean];
% generate full inference id vec
inf_id_vec = NaN(size(group_id_vec));
for i = 1:length(inf_id_vec)
    ind = find(group_id_vec(i)==gp_kni_array(:,1) & knirps_id_vec(i)==gp_kni_array(:,2));
    inf_id_vec(i) = inf_ids(ind);
end
z_flag_vec = z_index(inf_id_vec);

% generate result vectors and arrays
A_array = cat(3,inference_results.A_mat);
R_array = NaN(size(A_array));
use_flags = false(size(knirps_axis_vec));

for i = 1:size(A_array,3)
   A_temp = A_array(:,:,i);
   warning('') % Clear last warning message
   [R_temp, exit_flag] = logm(A_temp);
   [warnMsg, warnId] = lastwarn;
   R_array(:,:,i) = R_temp/20;
   use_flags(i) = all(isreal(R_temp(:))) && sum(R_temp(:)>=0)==2 && isempty(warnMsg);
end  
r_array = cat(2,inference_results.r);

pon_vec = reshape(A_array(2,1,:),[],1);
poff_vec = reshape(A_array(1,2,:),[],1);

kon_vec = reshape(R_array(2,1,:),[],1);
koff_vec = reshape(R_array(1,2,:),[],1);

group_index = unique(group_id_vec);
knirps_mean_vec = NaN(size(group_index));
knirps_ste_vec = NaN(size(group_index));
zeros_vec  = NaN(size(group_index));
kon_mean_vec = NaN(size(group_index));
kon_ste_vec = NaN(size(group_index));
koff_mean_vec = NaN(size(group_index));
koff_ste_vec = NaN(size(group_index));

pon_mean_vec = NaN(size(group_index));
pon_ste_vec = NaN(size(group_index));
poff_mean_vec = NaN(size(group_index));
poff_ste_vec = NaN(size(group_index));

r_mean_vec = NaN(size(group_index));
r_ste_vec = NaN(size(group_index));

for g = 1:length(group_index)
    knirps_mean_vec(g) = nanmean(knirps_axis_vec(group_id_vec==group_index(g)&use_flags));
    knirps_ste_vec(g) = nanstd(knirps_axis_vec(group_id_vec==group_index(g)&use_flags));
    
    zeros_vec(g) = unique(z_flag_vec(group_id_vec==group_index(g)&use_flags));
    
    kon_mean_vec(g) = nanmean(kon_vec(group_id_vec==group_index(g)&use_flags));
    kon_ste_vec(g) = nanstd(kon_vec(group_id_vec==group_index(g)&use_flags));
    
    koff_mean_vec(g) = nanmean(koff_vec(group_id_vec==group_index(g)&use_flags));
    koff_ste_vec(g) = nanstd(koff_vec(group_id_vec==group_index(g)&use_flags));
    
    pon_mean_vec(g) = nanmean(pon_vec(group_id_vec==group_index(g)&use_flags));
    pon_ste_vec(g) = nanstd(pon_vec(group_id_vec==group_index(g)&use_flags));
    
    poff_mean_vec(g) = nanmean(poff_vec(group_id_vec==group_index(g)&use_flags));
    poff_ste_vec(g) = nanstd(poff_vec(group_id_vec==group_index(g)&use_flags));
    
    r_mean_vec(g) = nanmean(r_array(2,group_id_vec==group_index(g)&use_flags));
    r_ste_vec(g) = nanstd(r_array(2,group_id_vec==group_index(g)&use_flags));
end    


%%
close all 
for p = 0:1
    for z = 0:1

        plot_filter = zeros_vec==z;% & ismember(inf_ids,[3]);


        kon_fig = figure;
        cmap = brewermap([],'Set2');
        hold on

        inf_indices =  find(plot_filter);
        inf_ids_short = inf_ids(inf_indices);
        inf_index = unique(inf_ids_short);
        s = [];
        for i = 1:length(inf_index)
            ft = inf_indices(inf_ids_short==inf_index(i));
            % calculate upper and lower bounds
            kni_upper = knirps_ste_vec(ft);
            kni_lower = -knirps_ste_vec(ft);

            if p == 0                
                kon_upper = kon_ste_vec(ft)*60;
                kon_lower = -kon_ste_vec(ft)*60;
                kon_mean = kon_mean_vec(ft)*60;
            else
                kon_mean = pon_mean_vec(ft);
                kon_upper = pon_ste_vec(ft);
                kon_lower = -pon_ste_vec(ft);
            end
                   
            errorbar(knirps_mean_vec(ft),kon_mean,kon_lower,kon_upper,kni_lower,kni_upper,'o','Color','k','Capsize',0)
            s(end+1) = scatter(knirps_mean_vec(ft),kon_mean,50,'MarkerFaceColor',cmap(inf_index(i)-z*1,:),'MarkerEdgeColor','k');
        end

        legend(s,'ON LOW','WT','ON CONST','Location','southwest')
        xlabel('average Knirps concentration')
        if p == 0
            ylabel('burst frequency (events per minute)')
        else
            ylabel('p_{on}')
        end
        
        xlim([6 11])
    %     ylim([0 4])
        set(gca,'FontSize',14)

        ax = gca;
        ax.YAxis(1).Color = 'k';
        ax.XAxis(1).Color = 'k';

        kon_fig.InvertHardcopy = 'off';
        set(gcf,'color','w');
        grid on
        set(gca,'Color',[228,221,209]/255) 

        if p == 0
            saveas(kon_fig,[FigurePath 'burst_freq_scatter_z' num2str(z) '.png'])
            saveas(kon_fig,[FigurePath 'burst_freq_scatter_z' num2str(z) '.pdf'])
        else
            saveas(kon_fig,[FigurePath 'pon_scatter_z' num2str(z) '.png'])
            saveas(kon_fig,[FigurePath 'pon_scatter_z' num2str(z) '.pdf'])
        end

        koff_fig = figure;
        cmap = brewermap([],'Set2');
        hold on
        s = [];
        for i = 1:length(inf_index)
            ft = inf_indices(inf_ids_short==inf_index(i));
            % calculate upper and lower bounds
            kni_upper = knirps_ste_vec(ft);
            kni_lower = -knirps_ste_vec(ft);
            
            if p == 0
                koff_mean = koff_mean_vec(ft)*60;
                koff_upper = koff_ste_vec(ft)*60;
                koff_lower = -koff_ste_vec(ft)*60;
            else
                koff_mean = poff_mean_vec(ft);
                koff_upper = poff_ste_vec(ft);
                koff_lower = -poff_ste_vec(ft);
            end

            errorbar(knirps_mean_vec(ft),koff_mean,koff_lower,koff_upper,kni_lower,kni_upper,'o','Color','k','Capsize',0)
            s(end+1) = scatter(knirps_mean_vec(ft),koff_mean,50,'MarkerFaceColor',cmap(inf_index(i)-z*1,:),'MarkerEdgeColor','k');
        end

        legend(s,'ON LOW','WT','ON CONST','Location','southwest')

        xlabel('average Knirps concentration')
        if p == 0
            ylabel('k_{off} (events per minute)')
            ylim([0 8])
        else
            ylabel('p_{off}')
        end
        xlim([6 11])
        
        set(gca,'FontSize',14)

        ax = gca;
        ax.YAxis(1).Color = 'k';
        ax.XAxis(1).Color = 'k';

        koff_fig.InvertHardcopy = 'off';
        set(gcf,'color','w');
        grid on
        set(gca,'Color',[228,221,209]/255) 

        if p == 0
            saveas(koff_fig,[FigurePath 'koff_scatter_z' num2str(z) '.png'])
            saveas(koff_fig,[FigurePath 'koff_scatter_z' num2str(z) '.pdf'])
        else
            saveas(koff_fig,[FigurePath 'poff_scatter_z' num2str(z) '.png'])
            saveas(koff_fig,[FigurePath 'poff_scatter_z' num2str(z) '.pdf'])
        end

        if p == 0

            r_fig = figure;
            cmap = brewermap([],'Set2');
            hold on

            s = [];
            for i = 1:length(inf_index)
                ft = inf_indices(inf_ids_short==inf_index(i));
                % calculate upper and lower bounds
                kni_upper = knirps_ste_vec(ft);
                kni_lower = -knirps_ste_vec(ft);

                r_upper = r_ste_vec(ft)*60;
                r_lower = -r_ste_vec(ft)*60;

                errorbar(knirps_mean_vec(ft),r_mean_vec(ft)*60,r_lower,r_upper,kni_lower,kni_upper,'o','Color','k','Capsize',0)
                s(end+1) = scatter(knirps_mean_vec(ft),r_mean_vec(ft)*60,50,'MarkerFaceColor',cmap(inf_index(i)-z*1,:),'MarkerEdgeColor','k');
            end

            legend(s,'ON LOW','WT','ON CONST','Location','southwest')

            xlabel('average Knirps concentration')
            ylabel('r (au per minute)')
            xlim([6 11])
            % ylim([0 8])
            set(gca,'FontSize',14)

            ax = gca;
            ax.YAxis(1).Color = 'k';
            ax.XAxis(1).Color = 'k';

            r_fig.InvertHardcopy = 'off';
            set(gcf,'color','w');
            grid on
            set(gca,'Color',[228,221,209]/255) 

            saveas(r_fig,[FigurePath 'r_scatter_z' num2str(z) '.png'])
            saveas(r_fig,[FigurePath 'r_scatter_z' num2str(z) '.pdf'])
        end
    end    
end    