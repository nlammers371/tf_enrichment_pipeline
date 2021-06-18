clear
close all

% designate project
projectName = 'Bcd-GFP_hbMS2-mCh_Airy_fast';%'Bcd-GFP_hbP2P-mCh';%
liveProject = LiveEnrichmentProject(projectName);
slashes = strfind(liveProject.figurePath,'\');
FigurePath = [liveProject.figurePath(1:slashes(end-1)) 'bleaching_checks' filesep];
mkdir(FigurePath)

% designate projects to compare
prefixList = {'2021-05-14-zeiss980_airyscan_fast_02','2021-05-26-zeiss980_airyscan_fast_03',...
              '2021-05-27-zeiss980_airyscan_fast_high_power_01'};

% load data
master_struct = struct;
for p = 1:length(prefixList)
    Prefix = prefixList{p};
    liveExperiment = LiveEnrichmentExperiment(Prefix);
    
    % load    
    load([liveExperiment.resultsFolder 'CompiledParticles.mat'])
    schnitzcells = getSchnitzcells(liveExperiment);
    
    % store
    master_struct(p).schnitzcells = schnitzcells;
    master_struct(p).Time = ElapsedTime;
    master_struct(p).nc14 = nc14;
    master_struct(p).Prefix = Prefix;
end

%% 
close all

ap_bounds = [38 43];
time_index = -3:35;
% sets_to_use = {[2 3], 1};
peg_time = 10;

for p = 1:length(master_struct)
    schnitzcells = master_struct(p).schnitzcells;
    bcd_vec = max(vertcat(schnitzcells.Fluo),[],2)';
    frame_vec = vertcat(schnitzcells.frames);
    if p~=2
        ap_vec = [schnitzcells.APpos]*100;
    else
        ap_vec = [schnitzcells.APPos]*100;
    end
    Time = round(master_struct(p).Time - master_struct(p).Time(master_struct(p).nc14));
    time_vec = Time(frame_vec);
%     if p == 1
%         time_vec = time_vec + 2;
%     elseif p == 2
    
    bcd_trend_raw = NaN(length(time_index),1);          
    % iterate through time
    for t = 1:length(time_index)
        t_ft = ap_vec >= ap_bounds(1) & ap_vec <= ap_bounds(2) & time_vec==time_index(t);
        bcd_trend_raw(t) = nanmean(bcd_vec(t_ft));
    end    
    master_struct(p).bcd_trend_raw = bcd_trend_raw;
    % find peak time
    [~,peak_ind] = nanmax(bcd_trend_raw);
    master_struct(p).time_norm = time_index - time_index(peak_ind) + 10;
    master_struct(p).bcd_trend_norm = bcd_trend_raw / bcd_trend_raw(peak_ind);
end
   

% make basic figures
raw_fig = figure;
hold on
cm1 = brewermap([],'Set2');
colormap(cm1);
for p = 1:length(master_struct)  
    plot(master_struct(p).time_norm,master_struct(p).bcd_trend_raw,'color',cm1(p,:),'LineWidth',1.5)
end
% formating
set(gca,'Fontsize',14);
legend('10% 488 (1)','10% 488 (2)','17% 488')
set(gca,'Color',[228,221,209]/255) 
grid on
raw_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

% axis labels
xlabel('minutes into nc14')
ylabel('[Bcd] at AP 38-43 (au)')
xlim([0 30])
% save
saveas(raw_fig,[FigurePath 'raw_bcd_trends.png'])

rel_fig = figure;
hold on
cm1 = brewermap([],'Set2');
colormap(cm1);
for p = 1:length(master_struct)  
    plot(master_struct(p).time_norm,master_struct(p).bcd_trend_norm,'color',cm1(p,:),'LineWidth',1.5)
end
% formating
set(gca,'Fontsize',14);
legend('10% 488 (1)','10% 488 (2)','17% 488')
set(gca,'Color',[228,221,209]/255) 
grid on
rel_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

% axis labels
xlabel('minutes into nc14')
ylabel('normalized [Bcd] at AP 38-43 (au)')
xlim([0 30])
ylim([0.2 1.05])
% save
saveas(rel_fig,[FigurePath 'rel_bcd_trends.png'])
