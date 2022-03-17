clear
close all

addpath(genpath('./lib'))
addpath(genpath('../../utilities'))

% Load data
dataRoot = 'S:\Nick\Dropbox\ProcessedEnrichmentData\';

% specify director for combined dataset
readPath = [dataRoot 'combinedOptoSets' filesep];

% generate directory to save stuff in
FigurePath = 'S:\Nick\Dropbox (Personal)\OptoPresentations\DecodedTraceAnalyses\';
mkdir(FigurePath)

% load dataset
load([readPath 'spot_struct.mat'],'spot_struct')

%%
% generate longform datasets with key observables

% generate bins
time_bins = 0:60:40*60;
ap_bins = -0.15:0.01:0.15;
knirps_bins = linspace(0,16,30);

% extract indexing vectors
ap_vec = [spot_struct.apPosNucleusNew];
% x_vec = [spot_struct.xPosNucleusNew];
% y_vec = [spot_struct.yPosNucleusNew];
kni_vec = [spot_struct.rawNCProteinNew];
time_vec = [spot_struct.timeNew];
fluo_vec = [spot_struct.fluoNew];

% extract viterbi state vec and build on/off switch flags to indicate when
% next observation is a state change
state_vec = [spot_struct.path_viterbi]-1;
off_switch_flags = NaN(size(state_vec));
on_switch_flags = NaN(size(state_vec));

start_i = 1;
for s = 1:length(spot_struct)
    v_path = spot_struct(s).path_viterbi-1;    
    
    d_path = [diff(v_path) NaN];
    on_flags = NaN(size(v_path));
    on_flags(v_path==0) = 0;
    on_flags(d_path==1) = 1;
    on_switch_flags(start_i:start_i+length(v_path)-1) = on_flags;
    
    off_flags = NaN(size(v_path));
    off_flags(v_path==1) = 0;
    off_flags(d_path==-1) = 1;
    off_switch_flags(start_i:start_i+length(v_path)-1) = off_flags;
    
    start_i = start_i+length(v_path);
end    


% extract id vec
master_set_id_vec =NaN(size(state_vec));
master_particle_id_vec = NaN(size(state_vec));
project_id_vec = NaN(size(state_vec));

start_i = 1;

for s = 1:length(spot_struct)
    n = length(spot_struct(s).timeNew);
    master_set_id_vec(start_i:start_i+n-1) = repelem(spot_struct(s).masterID,n);
    master_particle_id_vec(start_i:start_i+n-1) = repelem(spot_struct(s).particleID,n);
    project_id_vec(start_i:start_i+n-1) = repelem(spot_struct(s).projectID,n);
    start_i = start_i+n;
end
master_id_index = unique(master_set_id_vec);

% combine into longform table
temp_array = [double(master_particle_id_vec') double(master_set_id_vec') double(project_id_vec') double(time_vec')...
              double(ap_vec') double(kni_vec') double(fluo_vec') double(state_vec') double(on_switch_flags') double(off_switch_flags')];
var_names = {'particle_id', 'set_id', 'project_id', 'time', 'ap', 'knirps', 'fluo', 'promoter_state', 'on_switch','off_switch'};
longTable = array2table(temp_array,'VariableNames',var_names);

% save
writetable(longTable,[readPath 'longformTable.csv'])

%% make some figures
% assign observations to bins 
ap_bin_vec = discretize(ap_vec,ap_bins);
kni_bin_vec = discretize(kni_vec,knirps_bins);
time_bin_vec = discretize(time_vec,time_bins);


% initialize arrays to store MF maps
ap_time_state_array = NaN(length(time_bins)-1,length(ap_bins)-1,length(master_id_index));
ap_time_on_array = NaN(length(time_bins)-1,length(ap_bins)-1,length(master_id_index));
ap_time_off_array = NaN(length(time_bins)-1,length(ap_bins)-1,length(master_id_index));

ap_kni_state_array = NaN(length(knirps_bins)-1,length(ap_bins)-1,length(master_id_index));
ap_kni_on_array = NaN(length(time_bins)-1,length(ap_bins)-1,length(master_id_index));
ap_kni_off_array = NaN(length(time_bins)-1,length(ap_bins)-1,length(master_id_index));

minTime = 10*60;
% AP vs TIME and AP vs Knirps
for m = 1:length(master_id_index)
    m_filter = master_set_id_vec==master_id_index(m);
    for a = 1:length(ap_bins)-1
        ap_filter = ap_bin_vec>=a-1 & ap_bin_vec<=a+1;
        for t = 1:length(time_bins)-1          
            time_filter = time_bin_vec>=t-1 & time_bin_vec<=t+1;
            mat_filter = m_filter & ap_filter & time_filter;
            if sum(mat_filter)>=25
                ap_time_state_array(t,a,m) = nanmean(state_vec(mat_filter));
            end
            if sum(~isnan(on_switch_flags(mat_filter)))>=25
                ap_time_on_array(t,a,m) = nanmean(on_switch_flags(mat_filter));
            end
            if sum(~isnan(off_switch_flags(mat_filter)))>=25
                ap_time_off_array(t,a,m) = nanmean(off_switch_flags(mat_filter));
            end
        end
        
        for k = 1:length(knirps_bins)-1          
            kni_filter = kni_bin_vec>=k-1 & kni_bin_vec<=k+1 & time_vec>minTime;
            mak_filter = m_filter & ap_filter & kni_filter;
            if sum(mak_filter)>=25
                ap_kni_state_array(k,a,m) = nanmean(state_vec(mak_filter));
            end
            if sum(~isnan(on_switch_flags(mak_filter)))>=25
                ap_kni_on_array(k,a,m) = nanmean(on_switch_flags(mak_filter));
            end
            if sum(~isnan(off_switch_flags(mak_filter)))>=25
                ap_kni_off_array(k,a,m) = nanmean(off_switch_flags(mak_filter));
            end
        end
    end
end
%%

master_id_vec_short = [spot_struct.masterID];

ap_axis = ap_bins(1:end-1) + diff(ap_bins)/2;
knirps_axis = knirps_bins(1:end-1) + diff(knirps_bins)/2;
time_axis = (time_bins(1:end-1) + diff(time_bins)/2)/60;
close all
for m = 1:length(master_id_index)    
      
    kni_array = ap_kni_state_array(:,:,m);
    kni_fig = figure;
    cmap = flipud(brewermap([],'Spectral'));
    colormap(cmap)
    p = pcolor(flipud(kni_array));
    p.EdgeAlpha = 0.1;
    xt = xticks;
    yt = yticks;
    set(gca,'xticklabels',ap_axis(xt),'yticklabels',round(fliplr(knirps_axis(yt)),1))
    xlabel('relative ap position')
    ylabel('knirps concentration (AU)')
    h = colorbar;
    ylabel(h,'fraction active')
    saveas(kni_fig,[FigurePath 'frac_vs_kni' sprintf('%03d',m) '.png'])
    
end   
    
%%
close all
for m = 1:length(master_id_index)    
      
    time_array = ap_time_state_array(:,:,m);
    
    time_fig = figure;
    cmap = flipud(brewermap([],'Spectral'));
    colormap(cmap)
    p = pcolor(flipud(time_array));
    p.EdgeAlpha = 0.1;
    xt = xticks;
    yt = yticks;
    set(gca,'xticklabels',ap_axis(xt),'yticklabels',round(fliplr(time_axis(yt)),1))
    xlabel('relative ap position')
    ylabel('time (AU)')
    h = colorbar;
    ylabel(h,'fraction active')
    saveas(time_fig,[FigurePath 'frac_vs_time' sprintf('%03d',m) '.png'])
    
end   
