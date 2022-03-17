clear
close all 

projectName = 'Bcd-GFP_hbMS2-mCh_AiryscanTest';

% get paths 
liveProject = LiveEnrichmentProject(projectName);
resultsRoot = [liveProject.dataPath filesep];

% make fig path 
FigurePath = [liveProject.figurePath 'bleaching_tests' filesep];
mkdir(FigurePath)
% load data
load([resultsRoot filesep 'spot_struct.mat'])

%% First let's get a sense for the spans covered by each set
set_vec = [spot_struct.setID];
set_index = unique(set_vec);
gfp_power_index = [5 10 10 20];

time_cell = cell(1,length(set_index));
bcd_cell = cell(1,length(set_index));

for  s = 1:length(set_index)
    set_filter = set_vec == set_index(s);
    time_vec_long = round([spot_struct(set_filter).time],0);
    nc_vec = [spot_struct.nc];
    % renormalize time wrpt start of nc14
    nc14_time = nanmin(round([spot_struct(nc_vec==14&set_filter).time]));
    time_vec_long = time_vec_long - nc14_time;
    
    bcd_vec_long = round([spot_struct(set_filter).rawNCProtein],0);
    bcd_vec_long(bcd_vec_long<5e4) = NaN; % deal with intermmittent zero values...need to look into source of these
    time_index = unique(time_vec_long);
    % calculate average protein over time
    bcd_vec = NaN(size(time_index));
    for t = 1:length(time_index)
      bcd_vec(t) = nanmean(bcd_vec_long(time_vec_long==time_index(t)));
    end
    % record
    bcd_cell{s} = bcd_vec;
    time_cell{s} = time_index;
end

% plot
span_fig = figure;
hold on
for s = 1:length(set_index)
    plot(time_cell{s},bcd_cell{s});
end
xlabel('time (minutes)')
ylabel('average absolute GFP signal')
legend('5% power','10% power (1)','10% power (2)','20% power')
set(gca,'Fontsize',14)
grid on
saveas(span_fig,[FigurePath 'span_fig.png'])


%% Check Frame over Frame trends, first for a single set
set_vec = [spot_struct.setID];
min_len = 10;
sets_to_use = [2 3 4];
peg_frames = [37 8 32];
protein_struct = struct;

for s = 1:length(sets_to_use)

    protein_struct(s).set_id = sets_to_use(s);    
    set_filter = set_vec == protein_struct(s).set_id;
    peg_frame = peg_frames(s);
    nc_vec = [spot_struct.nc];
    
    % renormalize time wrpt start of nc14
%     nc14_time = nanmin(round([spot_struct(nc_vec==14&set_filter).time]));
    
    % generate array of pegged and unpegged protein trends
    protein_struct(s).time_index = unique(round([spot_struct(set_filter).time],0));
    protein_struct(s).peg_time = protein_struct(s).time_index(peg_frame);
    protein_struct(s).protein_pegged = NaN(length(protein_struct(s).time_index),sum(set_filter));
    protein_struct(s).protein_raw = NaN(length(protein_struct(s).time_index),sum(set_filter));
    
    iter = 1;
    for i = find(set_filter)
        time_temp = round(spot_struct(i).time,0);
        frames_temp = spot_struct(i).frames;
        if length(time_temp) >= min_len && any(time_temp==protein_struct(s).peg_time)                        
            protein_temp = spot_struct(i).rawNCProtein;
            protein_temp(protein_temp<=1e5) = NaN;
%             if sets_to_use(s) == 4 && min(frames_temp)<20
%               error('wtf')
%             end
            % check for artifact present in 20% set
            use_flag = 1;
%             if sets_to_use(s) == 4 && any(frames_temp==46)
%                 span_filter = 30:70;
%                 use_flag = all(diff(protein_temp(ismember(frames_temp,span_filter)))<=1e4);
%             end
            if use_flag
                time_filter = ismember(protein_struct(s).time_index,time_temp);
                protein_struct(s).protein_raw(time_filter,iter) = protein_temp;
                protein_struct(s).protein_pegged(time_filter,iter) = protein_temp / protein_temp(time_temp==protein_struct(s).peg_time);
            end
        end
        iter = iter + 1;
    end    
    protein_struct(s).time_index = protein_struct(s).time_index-protein_struct(s).peg_time;    
    
    protein_struct(s).protein_pegged_alt = nanmean(protein_struct(s).protein_raw,2);
    protein_struct(s).protein_pegged_alt = protein_struct(s).protein_pegged_alt / protein_struct(s).protein_pegged_alt(peg_frame);
end
        
        

close all
gfp_fig = figure;
hold on
lw_vec = [1 1 2.5];
for i = 1:length(protein_struct)
    plot(protein_struct(i).time_index/60,nanmean(protein_struct(i).protein_pegged,2),'LineWidth',lw_vec(i))
%     plot(protein_struct(i).time_index/60,protein_struct(i).protein_pegged_alt,'LineWidth',lw_vec(i))
end

xlabel('minutes from peg time')
ylabel('normalized nuclear GFP intensity (AU)')
legend('10% power (1)','10% power (2)' ,'20% power')
set(gca,'Fontsize',14)
grid on
xlim([-5 15])
saveas(gfp_fig,[FigurePath 'bcd_bleaching_test_10vs20.png'])


%% Compare 5 to 10%

sets_to_use = [1 2];
peg_frames = [29 1];
protein_struct = struct;

for s = 1:length(sets_to_use)

    protein_struct(s).set_id = sets_to_use(s);    
    set_filter = set_vec == protein_struct(s).set_id;
    peg_frame = peg_frames(s);
%     nc_vec = [spot_struct.nc];
    
    % renormalize time wrpt start of nc14
%     nc14_time = nanmin(round([spot_struct(nc_vec==14&set_filter).time]));
    
    % generate array of pegged and unpegged protein trends
    protein_struct(s).time_index = unique(round([spot_struct(set_filter).time],0));
    protein_struct(s).peg_time = protein_struct(s).time_index(peg_frame);
    protein_struct(s).protein_pegged = NaN(length(protein_struct(s).time_index),sum(set_filter));
    protein_struct(s).protein_raw = NaN(length(protein_struct(s).time_index),sum(set_filter));
    
    iter = 1;
    for i = find(set_filter)
        time_temp = round(spot_struct(i).time,0);
        if length(time_temp) >= min_len && any(time_temp==protein_struct(s).peg_time)    
            protein_temp = spot_struct(i).rawNCProtein;
            protein_temp(protein_temp<=1e5) = NaN;
            time_filter = ismember(protein_struct(s).time_index,time_temp);
            protein_struct(s).protein_raw(time_filter,iter) = protein_temp;
            protein_struct(s).protein_pegged(time_filter,iter) = protein_temp / protein_temp(time_temp==protein_struct(s).peg_time);
        end
        iter = iter + 1;
    end
    
    protein_struct(s).time_index = protein_struct(s).time_index-protein_struct(s).peg_time;
    
end
        
        

close all
gfp_fig = figure;
hold on
lw_vec = [1 1];
for i = 1:length(protein_struct)
    plot(protein_struct(i).time_index/60,nanmean(protein_struct(i).protein_pegged,2),'LineWidth',lw_vec(i))
end

xlabel('minutes from peg time')
ylabel('normalized nuclear GFP intensity (AU)')
legend('5% power','10% power')
set(gca,'Fontsize',14)
grid on
xlim([-5 10])
saveas(gfp_fig,[FigurePath 'bcd_bleaching_test_5v10.png'])

%% Look for Signatures of Bursting
set_id = 3;
set_filter = set_vec==set_id;
tracePath = [FigurePath filesep 'single_trace_plots_set' num2str(set_id) filesep];
mkdir(tracePath);

for i = find(set_filter)
    time_temp = spot_struct(i).timeInterp;
    fluo_temp = spot_struct(i).fluoInterp;
    fluo_raw = spot_struct(i).fluo;
    if sum(~isnan(fluo_raw)) > 20
      
        temp_fig = figure('Visible','off');        
        plot(time_temp/60,fluo_temp);
        
        set(gca,'Fontsize',14)
        xlabel('minutes into nc14')
        ylabel('MS2 spot intensity')
        
        saveas(temp_fig,[tracePath 'trace' sprintf('%03d',i) '.png'])
        close all
    end
end

%% Look for bleaching in mCherry channel
min_len = 10;
sets_to_use = [2 3];
peg_frames = [48 17];
ms2_struct = struct;

for s = 1:length(sets_to_use)

    ms2_struct(s).set_id = sets_to_use(s);    
    set_filter = set_vec == ms2_struct(s).set_id;
    peg_frame = peg_frames(s);
    nc_vec = [spot_struct.nc];
    
    % renormalize time wrpt start of nc14
%     nc14_time = nanmin(round([spot_struct(nc_vec==14&set_filter).time]));
    
    % generate array of pegged and unpegged protein trends
    ms2_struct(s).time_index = unique(round([spot_struct(set_filter).time],0));
    ms2_struct(s).peg_time = ms2_struct(s).time_index(peg_frame);
    ms2_struct(s).ms2_pegged = NaN(length(ms2_struct(s).time_index),sum(set_filter));
    ms2_struct(s).ms2_raw = NaN(length(ms2_struct(s).time_index),sum(set_filter));
    
    iter = 1;
    for i = find(set_filter)
        time_temp = round(spot_struct(i).time,0);
        frames_temp = spot_struct(i).frames;
        if length(time_temp) >= min_len && any(time_temp==ms2_struct(s).peg_time)                        
            ms2_temp = spot_struct(i).fluo;
%           
            time_filter = ismember(ms2_struct(s).time_index,time_temp);
            ms2_struct(s).ms2_raw(time_filter,iter) = ms2_temp;
            ms2_struct(s).ms2_pegged(time_filter,iter) = ms2_temp/ ms2_temp(time_temp==ms2_struct(s).peg_time);
            
        end
        iter = iter + 1;
    end    
    ms2_struct(s).time_index = ms2_struct(s).time_index-ms2_struct(s).peg_time;    
    
    ms2_struct(s).ms2_pegged_alt = nanmean(ms2_struct(s).ms2_raw,2);
    ms2_struct(s).ms2_pegged_alt = ms2_struct(s).ms2_pegged_alt / ms2_struct(s).ms2_pegged_alt(peg_frame);
end
        
        

close all
gfp_fig = figure;
hold on
lw_vec = [1 1 2.5];
for i = 1:length(ms2_struct)
%     plot(ms2_struct(i).time_index/60,nanmean(ms2_struct(i).ms2_pegged,2),'LineWidth',lw_vec(i))
    plot(ms2_struct(i).time_index/60,ms2_struct(i).ms2_pegged_alt,'LineWidth',lw_vec(i))
end

xlabel('minutes from peg time')
ylabel('normalized mCherry spot intensity (AU)')
legend('6% power','5% power')
set(gca,'Fontsize',14)
grid on
xlim([-10 15])
saveas(gfp_fig,[FigurePath 'mCherry_bleaching_test_10vs20.png'])

%%
fiveP_vec = ms2_struct(2).ms2_raw(:);
fiveP_vec = fiveP_vec(~isnan(fiveP_vec));

sixP_vec = ms2_struct(1).ms2_raw(:);
sixP_vec = sixP_vec(~isnan(sixP_vec));

hist_bins = linspace(0,6e5,50);

hist_fig = figure;
hold on
histogram(fiveP_vec,hist_bins,'Normalization','probability')
histogram(sixP_vec,hist_bins,'Normalization','probability')
legend('5% power','6% power')
xlabel('MS2 spot intensity')
ylabel('probability')
set(gca,'Fontsize',14)
saveas(hist_fig,[FigurePath 'mCherry_spot_hist.png'])

