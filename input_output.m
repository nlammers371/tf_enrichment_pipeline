% Script to contuct exploratory analyses regarding relationship between
% local concentration and locus activity
clear 
close all
% ID variable
project = 'Dl_Venus_snaBAC_mCherry';
ControlType = 'edge';
dropboxFolder = 'E:\Nick\Dropbox (Garcia Lab)\ProcessedEnrichmentData\';
dataPath = [dropboxFolder project '\'];
writePath = ['E:\Nick\Dropbox (Garcia Lab)\LocalEnrichmentFigures\' project '\'];
mkdir(writePath);
% load data
load([dataPath 'nucleus_struct_protein.mat'])


% extract protein, gene, fluorophore info
underscores = strfind(project,'_');
protein_name = project(1:underscores(1)-1);
protein_fluor = project(underscores(1)+1:underscores(2)-1);
gene_name = project(underscores(2)+1:underscores(3)-1);
if numel(underscores) == 3
    ind = numel(project);
else
    ind = underscores(4)-1;
end
gene_fluor = project(underscores(3)+1:end);

% basic info params
pixelSize = nucleus_struct_protein(1).PixelSize;
ROIRadius = .25; % radus (um) of region used to query and compare TF concentrations
distLim = .6;
n_boots = 100;
n_xcorr_lags = 15; % time steps

%%% Look at bleaching profiles for Protein and MCP channels
set_vec = [nucleus_struct_protein.setID];
set_index = unique([nucleus_struct_protein.setID]);

% generate time-averaged vectors for each set
time_index = unique(round([nucleus_struct_protein.time]/60));
fluo_avg = NaN(numel(time_index),numel(set_index));
protein_avg = NaN(numel(time_index),numel(set_index));
time_sigma = 1;

for i = 1:numel(set_index)    
    time = [nucleus_struct_protein(set_vec==set_index(i)).time]/60;
    fluo = [nucleus_struct_protein(set_vec==set_index(i)).fluo];
    protein = [nucleus_struct_protein(set_vec==set_index(i)).protein];
    for j = 1:numel(time_index)
        t_diffs = time - time_index(j);
        t_weights = exp(-.5*(t_diffs / time_sigma).^2);
        fluo_avg(j,i) = nansum(t_weights.*fluo) / nansum(t_weights);
        protein_avg(j,i) = nansum(t_weights.*protein) / nansum(t_weights);
    end
end

% make figure
time_trends = figure;
subplot(2,1,1)
hold on
lgd_str ={};
for i = 1:numel(set_index)
    plot(time_index,fluo_avg(:,i))    
    lgd_str = [lgd_str{:} {['set ' num2str(set_index(i))]}];
end
legend(lgd_str{:})
title(['Time-averaged trends: ' gene_fluor ' channel'])
grid on
xlabel('minutes')
ylabel('au')

subplot(2,1,2)
hold on
lgd_str ={};
for i = 1:numel(set_index)
    plot(time_index,protein_avg(:,i))    
    lgd_str = [lgd_str{:} {['set ' num2str(set_index(i))]}];
end
legend(lgd_str{:})
title(['Time-averaged trends: ' protein_fluor ' channel'])
grid on
xlabel('minutes')
ylabel('au')

saveas(time_trends,[writePath 'time_trend_fig.png'])

%%% Make some useful data vectors
% Generate distance vector for filtering snip stacks
dist_vec = [nucleus_struct_protein.(['spot_' ControlType '_dist_vec'])];
time_vec = [nucleus_struct_protein.time];
frame_vec = [nucleus_struct_protein.frames];
ap_vec = [nucleus_struct_protein.ap_vector];
fluo_vec = [nucleus_struct_protein.fluo];
mf_protein_vec = [nucleus_struct_protein.protein];
% make distance_filter
dist_filter = dist_vec*pixelSize >= distLim;
set_vec = [];
particle_vec = [];
snip_filter = [];
for i = 1:numel(nucleus_struct_protein)
    snip_frame_vec = nucleus_struct_protein(i).snip_frame_vec;  
    nc_frames = nucleus_struct_protein(i).frames;
    snip_filter = [snip_filter ismember(nc_frames,snip_frame_vec)];
    set_vec = [set_vec repelem(nucleus_struct_protein(i).setID,sum(ismember(nc_frames,snip_frame_vec)))];
    particle_vec = [particle_vec repelem(nucleus_struct_protein(i).ParticleID,sum(ismember(nc_frames,snip_frame_vec)))];
end  
snip_filter = snip_filter == 1;
% filter vectors
dist_vec = dist_vec(snip_filter&dist_filter);
fluo_vec = fluo_vec(snip_filter&dist_filter);
frame_vec = frame_vec(snip_filter&dist_filter);
time_vec = time_vec(snip_filter&dist_filter);
mf_protein_vec = mf_protein_vec(snip_filter&dist_filter);

% Snip stacks
spot_protein_snips = cat(3,nucleus_struct_protein.spot_protein_snips);
null_protein_snips = cat(3,nucleus_struct_protein.([ControlType '_null_protein_snips']));
spot_mcp_snips = cat(3,nucleus_struct_protein.spot_mcp_snips);
null_mcp_snips = cat(3,nucleus_struct_protein.([ControlType '_null_mcp_snips']));

% Make r reference array
snip_size = size(spot_protein_snips,1);
[y_ref, x_ref] = meshgrid(1:snip_size,1:snip_size);
r_ref = sqrt((x_ref-ceil(snip_size/2)).^2 + (y_ref-ceil(snip_size/2)).^2)*pixelSize;


r_nan = ones(size(r_ref));
r_nan(r_ref>ROIRadius) = NaN;
r_ft = repmat(r_nan,1,1,size(spot_protein_snips(:,:,dist_filter(snip_filter)),3));
spot_protein_vec = reshape(nanmean(nanmean(r_ft.*spot_protein_snips(:,:,dist_filter(snip_filter)),1),2),1,[]);
null_protein_vec = reshape(nanmean(nanmean(r_ft.*null_protein_snips(:,:,dist_filter(snip_filter)),1),2),1,[]);
% rescale mf reference vec
mf_protein_vec = mf_protein_vec * nanmean(null_protein_vec) / nanmean(mf_protein_vec) ;
particle_vec = particle_vec(dist_filter(snip_filter));

%%% Cross-correlation
% generate de-trended fluroescence and protein vectors
minute_vec = 1:50;
fluo_vec_norm = NaN(size(fluo_vec));
fluo_mean_vec = NaN(size(minute_vec));
fluo_std_vec = NaN(size(minute_vec));

spot_pt_norm = NaN(size(fluo_vec));
null_pt_norm = NaN(size(fluo_vec));
protein_spot_mean_vec = NaN(size(minute_vec));
protein_null_mean_vec = NaN(size(minute_vec));
protein_spot_std_vec = NaN(size(minute_vec));
protein_null_std_vec = NaN(size(minute_vec));
for i = 1:numel(fluo_mean_vec)
    t_filter = round(time_vec/60)==minute_vec(i);
    fluo_mean_vec(i) = nanmean(fluo_vec(t_filter));
    
    fluo_std_vec(i) = nanstd(fluo_vec(t_filter));
    protein_spot_mean_vec(i) = nanmean(spot_protein_vec(t_filter));
    protein_null_mean_vec(i) = nanmean(null_protein_vec(t_filter));
    protein_spot_std_vec(i) = nanstd(spot_protein_vec(t_filter));
    protein_null_std_vec(i) = nanstd(null_protein_vec(t_filter));
    
    fluo_vec_norm(t_filter) = (fluo_vec(t_filter) - fluo_mean_vec(i));%/ fluo_std_vec(i);
    spot_pt_norm(t_filter) = (spot_protein_vec(t_filter) - mf_protein_vec(t_filter)) ;%/ protein_spot_std_vec(i);
    null_pt_norm(t_filter) = (null_protein_vec(t_filter) - mf_protein_vec(t_filter)) ;%/ protein_null_std_vec(i);
end    

% compile array of trace and protein fragments
spot_pt_array = [];
null_pt_array = [];
% particle ID vectors
% particle_vec_snip = particle_vec(snip_filter==1);
particle_index = unique(particle_vec);
% iterate through particles and xcorr signatures for sufficiently long
% fragments
for i = 1:numel(particle_index)
    pt_ft = particle_vec == particle_index(i);
    frames = frame_vec(pt_ft);
    
    sp_pt_vec = spot_pt_norm(pt_ft);
    nl_pt_vec = null_pt_norm(pt_ft);
    fl_pt_vec = fluo_vec_norm(pt_ft);
    
    frames_d = diff(frames);
    % identify gaps
    gap_ids = unique([1 find(frames_d>1) numel(frames)]);
    xcorr_ids = find(diff(gap_ids) + 1 >= n_xcorr_lags);
    for j = xcorr_ids
        sppt = sp_pt_vec(gap_ids(j):gap_ids(j+1));
        nlpt = nl_pt_vec(gap_ids(j):gap_ids(j+1));
        flpt = fl_pt_vec(gap_ids(j):gap_ids(j+1));
    
        spot_pt_array = [spot_pt_array xcov(sppt,flpt,n_xcorr_lags,'unbiased')'];
        null_pt_array = [null_pt_array xcov(nlpt,flpt,n_xcorr_lags,'unbiased')'];        
    end
end
% use bootstrapping to estimate standard errors

id_vec = 1:size(spot_pt_array,2);
spot_pt_boots = NaN(size(spot_pt_array,1),n_boots);
null_pt_boots = NaN(size(spot_pt_array,1),n_boots);
for i = 1:n_boots
    boot_ids = randsample(id_vec,numel(id_vec),true);
    spot_boot = spot_pt_array(:,boot_ids);
    null_boot = null_pt_array(:,boot_ids);
    
    spot_pt_boots(:,i) = nanmean(spot_boot,2);
    null_pt_boots(:,i) = nanmean(null_boot,2);
end

spot_corr_mean = nanmean(spot_pt_boots,2);
spot_corr_se = nanstd(spot_pt_boots,[],2);

null_corr_mean = nanmean(null_pt_boots,2);
null_corr_se = nanstd(null_pt_boots,[],2);

t_res = nanmedian(diff(time_vec));

xcorr_fig = figure;
hold on
e1 = errorbar(t_res*[-fliplr(1:n_xcorr_lags) 0 1:n_xcorr_lags], spot_corr_mean, spot_corr_se);
e1.CapSize = 0;
e2 = errorbar(t_res*[-fliplr(1:n_xcorr_lags) 0 1:n_xcorr_lags], null_corr_mean, null_corr_se);
e2.CapSize = 0;
% xlim([-n_lags*t_res,0])
% ylim([-.15, .1])
grid on
xlabel('lag (seconds)')
ylabel('cross-correlation')
legend('active locus','control')
saveas(xcorr_fig,[writePath 'crosscorrelation.png'])

% Interesting. Cross-correlation indicates that there is a time-lagged
% relationship wherein local protein concentration preceds corresponding
% observed fluorescence by about 120 seconds

%%% Look at local enrichment preceding sustained high and low periods
% first iterate through traces to generate convolved
fluo_kernel = ones(1,5);
fluo_interp_vec = NaN(size(fluo_vec_norm));
fluo_conv_vec = NaN(size(fluo_vec_norm));
fluo_diff_conv_vec = NaN(size(fluo_vec_norm));
for i = 1:numel(particle_index)
    pt_filter = particle_vec==particle_index(i);
    fluo = fluo_vec_norm(pt_filter);
    time = time_vec(pt_filter);
    fluo_nan = conv(isnan(fluo),fluo_kernel,'same');
    nan_ft = fluo_nan>0&fluo_nan<=3&isnan(fluo);
    fluo_interp = fluo;
    if sum(~isnan(fluo)) > 3 
        fluo_inbterp = fluo;
        if sum(nan_ft) > 1
            fluo_interp(nan_ft) = interp1(time(~nan_ft),fluo(~nan_ft),time(nan_ft));
        end        
        fluo_conv = conv(fluo_interp,fluo_kernel,'same');
        fluo_diff_conv = conv([0 diff(fluo_interp)],fluo_kernel,'same');        
    else
        fluo_interp = NaN(size(fluo));
        fluo_conv = NaN(size(fluo));
        fluo_diff_conv = NaN(size(fluo));
    end
    fluo_interp_vec(pt_filter) = fluo_interp;
    fluo_conv_vec(pt_filter) = fluo_conv;
    fluo_diff_conv_vec(pt_filter) = fluo_diff_conv;    
end
%%%
% calculate percentiles
prctile_vec = 0:10:100;
prctile_vals_abs = NaN(size(prctile_vec));
prctile_vals_diff = NaN(size(prctile_vec));
n_xcorr_lags = 20;
for p = 1:numel(prctile_vec)
    prctile_vals_abs(p) = prctile(fluo_conv_vec,prctile_vec(p));
    prctile_vals_diff(p) = prctile(fluo_diff_conv_vec,prctile_vec(p));
end
% initialize vectors
null_rise_lookback = [];
spot_rise_lookback = [];
null_fall_lookback = [];
spot_fall_lookback = [];
null_high_lookback = [];
spot_high_lookback = [];
null_low_lookback = [];
spot_low_lookback = [];

for i = 1:numel(particle_index)
    % basic info
    ParticleID = particle_index(i); 
    pt_filter = particle_vec==ParticleID;
    frames = frame_vec(pt_filter);
    snip_frames = nucleus_struct_protein(i).snip_frame_vec;
    % do absolute numbers first
    fluo_conv = fluo_conv_vec(pt_filter);
    fluo_diff_conv = fluo_diff_conv_vec(pt_filter);
    % find frames with unusually high fluorescence
    high_ids = find(fluo_conv > prctile_vals_abs(end-1));
    for j = high_ids        
        lookback_frames = frames(j)-n_xcorr_lags:frames(j);
        null = NaN(size(lookback_frames));
        spot = NaN(size(lookback_frames));
        frame_ft = ismember(frame_vec,lookback_frames) & particle_vec==ParticleID;
        
        null(ismember(lookback_frames,frame_vec(frame_ft))) = null_pt_norm(frame_ft);
        spot(ismember(lookback_frames,frame_vec(frame_ft))) = spot_pt_norm(frame_ft);
        
        null_high_lookback = [null_high_lookback ; null];
        spot_high_lookback = [spot_high_lookback ; spot];
    end
     % find frames with unusually low fluorescence
    low_ids = find(fluo_conv <= prctile_vals_abs(2));
    for j = low_ids        
        lookback_frames = frames(j)-n_xcorr_lags:frames(j);
        null = NaN(size(lookback_frames));
        spot = NaN(size(lookback_frames));
        frame_ft = ismember(frame_vec,lookback_frames) & particle_vec==ParticleID;

        null(ismember(lookback_frames,frame_vec(frame_ft))) = null_pt_norm(frame_ft);
        spot(ismember(lookback_frames,frame_vec(frame_ft))) = spot_pt_norm(frame_ft);

        null_low_lookback = [null_low_lookback ; null];
        spot_low_lookback = [spot_low_lookback ; spot];
    end
    
    rise_ids = find(fluo_diff_conv > prctile_vals_diff(end-1));  
    for j = rise_ids        
        lookback_frames = frames(j)-n_xcorr_lags:frames(j);
        null = NaN(size(lookback_frames));
        spot = NaN(size(lookback_frames));
        frame_ft = ismember(frame_vec,lookback_frames) & particle_vec==ParticleID;

        null(ismember(lookback_frames,frame_vec(frame_ft))) = null_pt_norm(frame_ft);
        spot(ismember(lookback_frames,frame_vec(frame_ft))) = spot_pt_norm(frame_ft);

        null_rise_lookback = [null_rise_lookback ; null];
        spot_rise_lookback = [spot_rise_lookback ; spot];
    end
    
    fall_ids = find(fluo_diff_conv <= prctile_vals_diff(2));
    for j = fall_ids        
        lookback_frames = frames(j)-n_xcorr_lags:frames(j);
        null = NaN(size(lookback_frames));
        spot = NaN(size(lookback_frames));
        frame_ft = ismember(frame_vec,lookback_frames) & particle_vec==ParticleID;

        null(ismember(lookback_frames,frame_vec(frame_ft))) = null_pt_norm(frame_ft);
        spot(ismember(lookback_frames,frame_vec(frame_ft))) = spot_pt_norm(frame_ft);

        null_fall_lookback = [null_fall_lookback ; null];
        spot_fall_lookback = [spot_fall_lookback ; spot];
    end
end

% perform bootstrap sampling
rise_spot_mat = NaN(n_boots,n_xcorr_lags+1);
fall_spot_mat = NaN(n_boots,n_xcorr_lags+1);
high_spot_mat = NaN(n_boots,n_xcorr_lags+1);
low_spot_mat = NaN(n_boots,n_xcorr_lags+1);

rise_null_mat = NaN(n_boots,n_xcorr_lags+1);
fall_null_mat = NaN(n_boots,n_xcorr_lags+1);
high_null_mat = NaN(n_boots,n_xcorr_lags+1);
low_null_mat = NaN(n_boots,n_xcorr_lags+1);

for n = 1:n_boots
    rise_ids = randsample(1:size(spot_rise_lookback,1),size(spot_rise_lookback,1),true);
    rise_spot_mat(n,:) = nanmean(spot_rise_lookback(rise_ids,:));
    rise_null_mat(n,:) = nanmean(null_rise_lookback(rise_ids,:));   
    
    fall_ids = randsample(1:size(spot_fall_lookback,1),size(spot_fall_lookback,1),true);
    fall_spot_mat(n,:) = nanmean(spot_fall_lookback(fall_ids,:));
    fall_null_mat(n,:) = nanmean(null_fall_lookback(fall_ids,:));
    
    high_ids = randsample(1:size(spot_high_lookback,1),size(spot_high_lookback,1),true);
    high_spot_mat(n,:) = nanmean(spot_high_lookback(high_ids,:));
    high_null_mat(n,:) = nanmean(null_high_lookback(high_ids,:));
    
    low_ids = randsample(1:size(spot_low_lookback,1),size(spot_low_lookback,1),true);
    low_spot_mat(n,:) = nanmean(spot_low_lookback(low_ids,:));
    low_null_mat(n,:) = nanmean(null_low_lookback(low_ids,:));
end
% rise
rise_spot_mean = nanmean(rise_spot_mat) - nanmean(rise_spot_mat(:));
rise_spot_ste = nanstd(rise_spot_mat);

rise_null_mean = nanmean(rise_null_mat) - nanmean(rise_null_mat(:));
rise_null_ste = nanstd(rise_null_mat);

rise_delta_mean = nanmean(rise_spot_mat-rise_null_mat);
rise_delta_ste = nanstd(rise_spot_mat-rise_null_mat);

% fall
fall_spot_mean = nanmean(fall_spot_mat)-nanmean(fall_spot_mat(:));
fall_spot_ste = nanstd(fall_spot_mat);

fall_null_mean = nanmean(fall_null_mat)-nanmean(fall_null_mat(:));
fall_null_ste = nanstd(fall_null_mat);

fall_delta_mean = nanmean(fall_spot_mat-fall_null_mat);
fall_delta_ste = nanstd(fall_spot_mat-fall_null_mat);

% high
high_spot_mean = nanmean(high_spot_mat)-nanmean(high_spot_mat(:));
high_spot_ste = nanstd(high_spot_mat);

high_null_mean = nanmean(high_null_mat)-nanmean(high_null_mat(:));
high_null_ste = nanstd(high_null_mat);

high_delta_mean = nanmean(high_spot_mat-high_null_mat);
high_delta_ste = nanstd(high_spot_mat-high_null_mat);

% low
low_spot_mean = nanmean(low_spot_mat)-nanmean(low_spot_mat(:));
low_spot_ste = nanstd(low_spot_mat);

low_null_mean = nanmean(low_null_mat)-nanmean(low_null_mat(:));
low_null_ste = nanstd(low_null_mat);

low_delta_mean = nanmean(low_spot_mat-low_null_mat);
low_delta_ste = nanstd(low_spot_mat-low_null_mat);

plot1_cell = {'low_spot','high_spot','high_delta','fall_spot','rise_spot','rise_delta'};
plot2_cell = {'low_null','high_null','low_delta','fall_null','rise_null','fall_delta'};

for i = 1:numel(plot1_cell)
    fig = figure;
    hold on
    e1 = errorbar(t_res*(-n_xcorr_lags:0),eval([plot1_cell{i} '_mean']),eval([plot1_cell{i} '_ste']));
    e1.CapSize = 0;
    e2 = errorbar(t_res*(-n_xcorr_lags:0),eval([plot2_cell{i} '_mean']),eval([plot2_cell{i} '_ste']));
    e2.CapSize = 0;
    xlabel('lag (seconds)')
    ylabel([protein_name '-' protein_fluor ' (locus minus control) (au)'])
    grid on
    legend(plot1_cell{i},plot2_cell{i},'Interpreter', 'none')
    title([protein_name '-' protein_fluor ' Enrichment Preceding Activity'])
    saveas(fig,[writePath plot1_cell{i} '_' plot2_cell{i} '_lookback.png'])
end


%%% Look at enrichment at start and end of trace lifetime
latest_start = 10*60;
earliest_end = 30*60;

spot_start_lookback = [];
null_start_lookback = [];

spot_stop_lookback = [];
null_stop_lookback = [];

nc_pt_index = [nucleus_struct_protein.ParticleID];
for i = 1:numel(particle_index)
    % basic info
    ParticleID = particle_index(i); 
    pt_filter = particle_vec==ParticleID;
    frames = frame_vec(pt_filter);
    all_time = nucleus_struct_protein(nc_pt_index==ParticleID).time;
    % do absolute numbers first
    time = time_vec(pt_filter);
    fluo = fluo_vec_norm(pt_filter);
    fluo_conv = fluo_conv_vec(pt_filter);
    fluo_diff_conv = fluo_diff_conv_vec(pt_filter);
    % get timing info
    start_time = time(find(~isnan(fluo),1));
    if isempty(start_time)
        continue
    end
    start_lag = start_time - all_time(1);
    end_lag = all_time(end) - time(find(~isnan(fluo),1,'last'));

    if start_lag > 120 && start_time <=latest_start
        lookahead_frames = frames(find(~isnan(fluo),1)):frames(find(~isnan(fluo),1))+n_xcorr_lags; 
        null = NaN(size(lookahead_frames));
        spot = NaN(size(lookahead_frames));
        frame_ft = ismember(frame_vec,lookahead_frames) & particle_vec==ParticleID;
        
        null(ismember(lookahead_frames,frame_vec(frame_ft))) = null_pt_norm(frame_ft);
        spot(ismember(lookahead_frames,frame_vec(frame_ft))) = spot_pt_norm(frame_ft);
        
        null_start_lookback = [null_start_lookback ; null];
        spot_start_lookback = [spot_start_lookback ; spot];
    end
    if end_lag >= 120
        lookahead_frames = frames(find(~isnan(fluo),1,'last'))-n_xcorr_lags:frames(find(~isnan(fluo),1,'last')); 
        null = NaN(size(lookahead_frames));
        spot = NaN(size(lookahead_frames));
        frame_ft = ismember(frame_vec,lookahead_frames) & particle_vec==ParticleID;
        
        null(ismember(lookahead_frames,frame_vec(frame_ft))) = null_pt_norm(frame_ft);
        spot(ismember(lookahead_frames,frame_vec(frame_ft))) = spot_pt_norm(frame_ft);
        
        null_stop_lookback = [null_stop_lookback ; null];
        spot_stop_lookback = [spot_stop_lookback ; spot];
    end
end


start_spot_mat = NaN(n_boots,n_xcorr_lags+1);
start_null_mat = NaN(n_boots,n_xcorr_lags+1);

stop_spot_mat = NaN(n_boots,n_xcorr_lags+1);
stop_null_mat = NaN(n_boots,n_xcorr_lags+1);

for n = 1:n_boots
    start_ids = randsample(1:size(spot_start_lookback,1),size(spot_start_lookback,1),true);
    start_spot_mat(n,:) = nanmean(spot_start_lookback(start_ids,:));
    start_null_mat(n,:) = nanmean(null_start_lookback(start_ids,:));
    
    stop_ids = randsample(1:size(spot_stop_lookback,1),size(spot_stop_lookback,1),true);
    stop_spot_mat(n,:) = nanmean(spot_stop_lookback(stop_ids,:));
    stop_null_mat(n,:) = nanmean(null_stop_lookback(stop_ids,:));
end

start_spot_mean = nanmean(start_spot_mat)-nanmean(start_spot_mat(:));
start_spot_ste = nanstd(start_spot_mat);
start_null_mean = nanmean(start_null_mat)-nanmean(start_null_mat(:));
start_null_ste = nanstd(start_null_mat);
start_delta_mean = nanmean(start_spot_mat-start_null_mat);
start_delta_ste = nanstd(start_spot_mat-start_null_mat);

stop_spot_mean = nanmean(stop_spot_mat)-nanmean(stop_spot_mat(:));
stop_spot_ste = nanstd(stop_spot_mat);
stop_null_mean = nanmean(stop_null_mat)-nanmean(stop_null_mat(:));
stop_null_ste = nanstd(stop_null_mat);
stop_delta_mean = nanmean(stop_spot_mat-stop_null_mat);
stop_delta_ste = nanstd(stop_spot_mat-stop_null_mat);

start_delta_fig = figure;
hold on
e1 = errorbar(t_res*(0:n_xcorr_lags),start_delta_mean,start_delta_ste);
e1.CapSize = 0;
xlabel('lead (seconds)')
ylabel([protein_name '-' protein_fluor ' (locus minus control) (au)'])
grid on
title([protein_name '-' protein_fluor 'Enrichment Following Onset of Activity'])
saveas(start_delta_fig,[writePath 'start_delta_lookahead.png'])

start_abs_fig = figure;
hold on
e1 = errorbar(t_res*(0:n_xcorr_lags),start_spot_mean,start_spot_ste);
e1.CapSize = 0;
e2 = errorbar(t_res*(0:n_xcorr_lags),start_null_mean,start_null_ste);
e2.CapSize = 0;
legend('start_spot','start_null','Interpreter','none')
xlabel('lead (seconds)')
ylabel([protein_name '-' protein_fluor ' concentration (au)'])
grid on
title([protein_name '-' protein_fluor 'Enrichment Following Onset of Activity'])
saveas(start_abs_fig,[writePath 'start_abs_lookahead.png'])

stop_fig = figure;
hold on
e1 = errorbar(t_res*(-n_xcorr_lags:0),stop_delta_mean,stop_delta_ste);
e1.CapSize = 0;

xlabel('lag (seconds)')
ylabel([protein_name '-' protein_fluor ' (locus minus control) (au)'])
grid on
title([protein_name '-' protein_fluor 'Enrichment Preceding Cessation of Activity'])
saveas(stop_fig,[writePath 'stop_delta_lookback.png'])

stop_abs_fig = figure;
hold on
e1 = errorbar(t_res*(0:n_xcorr_lags),stop_spot_mean,stop_spot_ste);
e1.CapSize = 0;
e2 = errorbar(t_res*(0:n_xcorr_lags),stop_null_mean,stop_null_ste);
e2.CapSize = 0;
legend('stop_spot','stop_null','Interpreter','none')
xlabel('lag (seconds)')
ylabel([protein_name '-' protein_fluor ' concentration (au)'])
grid on
title([protein_name '-' protein_fluor 'Enrichment Preceding Cessation of Activity'])
saveas(stop_abs_fig,[writePath 'stop_abs_lookahead.png'])

%%

% %% Break data up into 10 bins by spot fluorescence and look at snips 7 steps preceding
% prctile_vec = 0:10:100;
% prctile_vals = NaN(size(prctile_vec));
% prctile_grp_vec = NaN(size(fluo_vec_snip));
% for p = 1:numel(prctile_vec)
%     prctile_vals(p) = prctile(fluo_vec_snip,prctile_vec(p));
% end
% for p = 1:numel(prctile_vec)-1
%     prctile_grp_vec(fluo_vec_snip>=prctile_vals(p)&fluo_vec_snip<prctile_vals(p+1)) = p;
% end
% 
% % build map between fluo frames and lagged snips
% spot_protein_snips_norm = spot_protein_snips ./ reshape(mf_protein_vec(snip_filter==1),1,1,[]);
% spot_protein_snips_norm2 = spot_protein_snips ./ nanmean(nanmean(null_protein_snips,1),2);
% null_protein_snips_norm = null_protein_snips ./ reshape(mf_protein_vec(snip_filter==1),1,1,[]);
% frame_lag = 7;
% frame_index_vec = 1:numel(prctile_grp_vec);
% snip_map_vec = NaN(size(frame_index_vec));
% for i = 1:numel(frame_index_vec)
%     frame = frame_vec_snip(i);
%     particle = particle_vec(i);
%     ind = find(particle_vec==particle&frame_vec_snip==frame+frame_lag);
%     if ~isempty(ind)
%         snip_map_vec(i) = ind;
%     end
% end
% %%%
% % calcualte average radial profiles for each percentile
% r_index = 0:.1:2;
% r_sigma = .1;
% 
% spot_mean_radius = NaN(numel(prctile_vec)-1,numel(r_index),n_boots);
% spot_mean_radius2 = NaN(numel(prctile_vec)-1,numel(r_index),n_boots);
% null_mean_radius = NaN(numel(prctile_vec)-1,numel(r_index),n_boots);
% tic   
% for p = 1:numel(prctile_vec)-1
%     indices = find(prctile_grp_vec==p);    
%     spot_stack = spot_protein_snips_norm(:,:,ismember(snip_map_vec,indices));
%     spot_stack2 = spot_protein_snips_norm2(:,:,ismember(snip_map_vec,indices));
%     null_stack = null_protein_snips_norm(:,:,ismember(snip_map_vec,indices));
%     for r = 1:numel(r_index)
%         r_weights = repmat(exp(-.5*((r_ref-r_index(r))/r_sigma).^2),1,1,size(spot_stack,3));
%         denom = nansum(r_weights(:));
%         spot_mean_radius(p,r,n) = nansum(r_weights(:).*spot_stack(:)) / denom;
%         spot_mean_radius2(p,r,n) = nansum(r_weights(:).*spot_stack2(:)) / denom;
%         null_mean_radius(p,r,n) = nansum(r_weights(:).*null_stack(:)) / denom;
%     end
% end
% toc

% spot_mean_stack = NaN(size(spot_protein_snips,1),size(spot_protein_snips,1),numel(prctile_vec));
% null_mean_stack = NaN(size(spot_protein_snips,1),size(spot_protein_snips,1),numel(prctile_vec));
% tic   
% for p = 1:numel(prctile_vec)
%     indices = find(prctile_grp_vec==p);    
%     spot_stack = spot_protein_snips_norm(:,:,ismember(snip_map_vec,indices));
%     null_stack = null_protein_snips_norm(:,:,ismember(snip_map_vec,indices));
% 
%     spot_mean_stack(:,:,p) = nanmean(spot_stack,3);
%     null_mean_stack(:,:,p) = nanmean(null_stack,3);   
% end


% % find avverage and ste
% spot_r_mean = nanmean(spot_mean_radius,3);
% spoit_r_ste = nanstd(spot_mean_radius,[],3);
% 
% spot_r_mean2 = nanmean(spot_mean_radius2,3);
% spoit_r_ste2 = nanstd(spot_mean_radius2,[],3);
% 
% null_r_mean = nanmean(null_mean_radius,3);
% null_r_ste = nanstd(null_mean_radius,[],3);

% figure;
% hold on
% inc_vec = 1:13:128;
% cm = jet(128);
% for i = 1:numel(prctile_vec)-1
%     plot(r_index,spot_r_mean2(i,:),'Color',cm(inc_vec(i),:))
% end
% colormap(cm(inc_vec,:))
% colorbar
