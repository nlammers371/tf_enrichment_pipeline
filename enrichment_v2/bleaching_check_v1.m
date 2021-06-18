clear
close all

% designate project
projectName = 'Bcd-GFP_hbMS2-mCh_Airy_fast';%'Bcd-GFP_hbP2P-mCh';%
liveProject = LiveEnrichmentProject(projectName);
slashes = strfind(liveProject.figurePath,'\');
FigurePath = [liveProject.figurePath(1:slashes(end-1)) 'bleaching_checks' filesep];
mkdir(FigurePath)

% designate projects to compare
prefixList = {'Bcd-GFP_hbMS2-mCh_Airy_fast','Bcd-GFP_hbMS2-mCh_AiryHP'};

% load data
master_struct = struct;
for p = 1:length(prefixList)
    projectName = prefixList{p};
    liveProject = LiveEnrichmentExperiment(projectName);
    resultsRoot = [liveProject.dataPath];
    % load    
    load([resultsRoot 'spot_struct_protein.mat'])
    load([resultsRoot 'proteinSamplingInfo.mat'])
    % store
    master_struct(p).spot_struct_protein = spot_struct_protein;
    master_struct(p).projectName = projectName;
end

%% 
close all

ap_bounds = [38 43];
time_index = 1:30;
% sets_to_use = {[2 3], 1};
peg_time = 10;

for p = 1:length(master_struct)
    spot_struct_protein = master_struct(p).spot_struct_protein;
    bcd_trend_raw = NaN(length(time_index),length(sets_to_use{p}));
    set_vec = [spot_struct_protein.setID];
    
    for s = 1:length(sets_to_use{p})
        setID = sets_to_use{p}(s);
        % get basic vectors
        time_vec = round([spot_struct_protein(set_vec==setID).time]/60);
        if p == 1 && s==2
            time_vec = time_vec+2;
        end
        ap_vec = [spot_struct_protein(set_vec==setID).APPosParticle];
        nc_vec = [spot_struct_protein(set_vec==setID).nuclear_protein_vec];
        % iterate through time
        for t = 1:length(time_index)
            t_ft = ap_vec >= ap_bounds(1) & ap_vec <= ap_bounds(2) & time_vec==time_index(t);
            bcd_trend_raw(t,s) = nanmean(nc_vec(t_ft));
        end
    end
    master_struct(p).bcd_trend_raw = bcd_trend_raw;
end
    
bcd_raw_vec = [master_struct.bcd_trend_raw];
bcd_norm_vec = bcd_raw_vec ./ bcd_raw_vec(time_index==peg_time,:);

% make basic figures
raw_fig = figure;
hold on
cm1 = brewermap([],'Set2');
colormap(cm1);
for t = 1:size(bcd_raw_vec,2)    
    plot(time_index,bcd_raw_vec(:,t),'color',cm1(t,:),'LineWidth',1.5)
end
% formating
set(gca,'Fontsize',14);
set(gca,'Color',[228,221,209]/255) 
grid on
raw_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

% axis labels
xlabel('minutes into nc14')
ylabel('[Bcd] (au)')

% save
saveas(raw_fig,[FigurePath 'raw_bcd_trends.png'])

rel_fig = figure;
hold on
cm1 = brewermap([],'Set2');
colormap(cm1);
for t = 1:size(bcd_raw_vec,2)    
    plot(time_index,bcd_norm_vec(:,t),'color',cm1(t,:),'LineWidth',1.5)
end
% formating
set(gca,'Fontsize',14);
set(gca,'Color',[228,221,209]/255) 
grid on
rel_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

% axis labels
xlabel('minutes into nc14')
ylabel('normalized [Bcd] (au)')

% save
saveas(rel_fig,[FigurePath 'rel_bcd_trends.png'])
