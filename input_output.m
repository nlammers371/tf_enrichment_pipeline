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
ROIRadius = .3; % radus (um) of region used to query and compare TF concentrations
distLim = .6;
%% Look at bleaching profiles for Protein and MCP channels
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

%% Make some useful data vectors
% Generate distance vector for filtering snip stacks
dist_vec = [nucleus_struct_protein.(['spot_' ControlType '_dist_vec'])];
time_vec = [nucleus_struct_protein.time];
frame_vec = [nucleus_struct_protein.frames];
ap_vec = [nucleus_struct_protein.ap_vector];
fluo_vec = [nucleus_struct_protein.fluo];
mf_protein_vec = [nucleus_struct_protein.protein];

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
r_ft = repmat(r_nan,1,1,size(spot_protein_snips,3));
spot_protein_vec = reshape(nanmean(nanmean(r_ft.*spot_protein_snips,1),2),1,[]);
null_protein_vec = reshape(nanmean(nanmean(r_ft.*null_protein_snips,1),2),1,[]);
% make some regression variables


%% Run some basic regressions

[beta,~,~,CovB] = mvregress([ones(sum(snip_filter),1) dist_vec(snip_filter==1)' ...
    time_vec(snip_filter==1)' null_protein_vec' spot_protein_vec'],fluo_vec(snip_filter==1)');
errors = sqrt(diag(CovB));
coeff1 = beta ./ errors;

% so now enrichment signature and fluo are negatively correlated...

% try a normalized protein vectors
spot_pt_norm = spot_protein_vec ./ mf_protein_vec(snip_filter==1);
[beta,Sigma,E,CovB] = mvregress([ones(sum(snip_filter),1) dist_vec(snip_filter==1)' ...
    spot_pt_norm'],fluo_vec(snip_filter==1)');
errors = sqrt(diag(CovB));
coeff2 = beta ./ errors;

% Okay then. negative correlation it is

%% Cross-correlation
n_lags = 15; % time steps

spot_pt_norm = spot_protein_vec ./ mf_protein_vec(snip_filter==1);
null_pt_norm = null_protein_vec ./ mf_protein_vec(snip_filter==1);

% compile array of trace and protein fragments
spot_pt_array = [];
null_pt_array = [];
% particle ID vectors
% particle_vec_snip = particle_vec(snip_filter==1);
fluo_vec_snip = fluo_vec(snip_filter==1);
time_vec_snip = time_vec(snip_filter==1);
frame_vec_snip = frame_vec(snip_filter==1);
particle_index = unique(particle_vec);
% iterate through particles and xcorr signatures for sufficiently long
% fragments
for i = 1:numel(particle_index)
    pt_ft = particle_vec == particle_index(i);
    frames = frame_vec_snip(pt_ft);
    
    sp_pt_vec = spot_pt_norm(pt_ft);
    nl_pt_vec = null_pt_norm(pt_ft);
    fl_pt_vec = fluo_vec_snip(pt_ft);
    
    frames_d = diff(frames);
    % identify gaps
    gap_ids = unique([1 find(frames_d>1) numel(frames)]);
    xcorr_ids = find(diff(gap_ids) + 1 >= n_lags);
    for j = xcorr_ids
        sppt = sp_pt_vec(gap_ids(j):gap_ids(j+1));
        nlpt = nl_pt_vec(gap_ids(j):gap_ids(j+1));
        flpt = fl_pt_vec(gap_ids(j):gap_ids(j+1));
    
        spot_pt_array = [spot_pt_array xcov(sppt,flpt,n_lags,'unbiased')'];
        null_pt_array = [null_pt_array xcov(nlpt,flpt,n_lags,'unbiased')'];        
    end
end
% use bootstrapping to estimate standard errors
n_boots = 100;
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
e1 = errorbar(t_res*[-fliplr(1:n_lags) 0 1:n_lags], spot_corr_mean, spot_corr_se);
e1.CapSize = 0;
e2 = errorbar(t_res*[-fliplr(1:n_lags) 0 1:n_lags], null_corr_mean, null_corr_se);
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

%% Look at local enrichment preceding sustained high and low periods

% first iterate through traces to generate convolved
fluo_kernel = ones(1,5);
for i = 1:numel(nucleus_struct_protein)
    fluo = nucleus_struct_protein(i).fluo;
    time = nucleus_struct_protein(i).time;
    fluo_nan = conv(isnan(fluo),fluo_kernel,'same');
    nan_ft = fluo_nan>0&fluo_nan<=3&isnan(fluo);
    fluo_interp = fluo;
    if sum(~isnan(fluo)) > 3 && sum(nan_ft) > 1
        fluo_interp(nan_ft) = interp1(time(~nan_ft),fluo(~nan_ft),time(nan_ft));
        fluo_conv = conv(fluo_interp,fluo_kernel,'same');
        fluo_diff_conv = conv([0 diff(fluo_interp)],fluo_kernel,'same');        
    else
        fluo_interp = NaN(size(fluo));
        fluo_conv = NaN(size(fluo));
        fluo_diff_conv = NaN(size(fluo));
    end
    nucleus_struct_protein(i).fluo_interp = fluo_interp;
    nucleus_struct_protein(i).fluo_conv = fluo_conv;
    nucleus_struct_protein(i).fluo_diff_conv = fluo_diff_conv;    
end

% remove spots that are too close to edge
dist_ft = dist_vec(snip_filter==1)*pixelSize >= distLim;
null_protein_vec_dist = null_protein_vec;
null_protein_vec_dist(~dist_ft) = NaN;
spot_protein_vec_dist = spot_protein_vec;
spot_protein_vec_dist(~dist_ft) = NaN;

% calculate percentiles
snip_ref_frame_vec = [nucleus_struct_protein.snip_frame_vec];
prctile_vec = 0:10:100;
prctile_vals_abs = NaN(size(prctile_vec));
prctile_vals_diff = NaN(size(prctile_vec));
n_lags = 20;
for p = 1:numel(prctile_vec)
    prctile_vals_abs(p) = prctile([nucleus_struct_protein.fluo_conv],prctile_vec(p));
    prctile_vals_diff(p) = prctile([nucleus_struct_protein.fluo_diff_conv],prctile_vec(p));
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

for i = 1:numel(nucleus_struct_protein)
    % basic info
    ParticleID = nucleus_struct_protein(i).ParticleID; 
    frames = nucleus_struct_protein(i).frames;
    snip_frames = nucleus_struct_protein(i).snip_frame_vec;
    % do absolute numbers first
    fluo_conv = nucleus_struct_protein(i).fluo_conv;
    fluo_diff_conv = nucleus_struct_protein(i).fluo_diff_conv;
    high_ids = find(fluo_conv > prctile_vals_abs(end-1));
    for j = high_ids        
        lookback_frames = frames(j)-n_lags:frames(j);
        null = NaN(size(lookback_frames));
        spot = NaN(size(lookback_frames));
        frame_ft = ismember(snip_ref_frame_vec,lookback_frames) & particle_vec==ParticleID;
        
        null(ismember(lookback_frames,snip_ref_frame_vec(frame_ft))) = null_protein_vec_dist(frame_ft);
        spot(ismember(lookback_frames,snip_ref_frame_vec(frame_ft))) = spot_protein_vec_dist(frame_ft);
        
        null_high_lookback = [null_high_lookback ; null];
        spot_high_lookback = [spot_high_lookback ; spot];
    end
    low_ids = find(fluo_conv <= prctile_vals_abs(2));
    for j = low_ids        
        lookback_frames = frames(j)-n_lags:frames(j);
        null = NaN(size(lookback_frames));
        spot = NaN(size(lookback_frames));
        frame_ft = ismember(snip_ref_frame_vec,lookback_frames) & particle_vec==ParticleID;

        null(ismember(lookback_frames,snip_ref_frame_vec(frame_ft))) = null_protein_vec_dist(frame_ft);
        spot(ismember(lookback_frames,snip_ref_frame_vec(frame_ft))) = spot_protein_vec_dist(frame_ft);

        null_low_lookback = [null_low_lookback ; null];
        spot_low_lookback = [spot_low_lookback ; spot];
    end
    
    rise_ids = find(fluo_diff_conv > prctile_vals_diff(end-1));  
    for j = rise_ids        
        lookback_frames = frames(j)-n_lags:frames(j);
        null = NaN(size(lookback_frames));
        spot = NaN(size(lookback_frames));
        frame_ft = ismember(snip_ref_frame_vec,lookback_frames) & particle_vec==ParticleID;

        null(ismember(lookback_frames,snip_ref_frame_vec(frame_ft))) = null_protein_vec_dist(frame_ft);
        spot(ismember(lookback_frames,snip_ref_frame_vec(frame_ft))) = spot_protein_vec_dist(frame_ft);

        null_rise_lookback = [null_rise_lookback ; null];
        spot_rise_lookback = [spot_rise_lookback ; spot];
    end
    
    fall_ids = find(fluo_diff_conv <= prctile_vals_diff(2));
    for j = fall_ids        
        lookback_frames = frames(j)-n_lags:frames(j);
        null = NaN(size(lookback_frames));
        spot = NaN(size(lookback_frames));
        frame_ft = ismember(snip_ref_frame_vec,lookback_frames) & particle_vec==ParticleID;

        null(ismember(lookback_frames,snip_ref_frame_vec(frame_ft))) = null_protein_vec_dist(frame_ft);
        spot(ismember(lookback_frames,snip_ref_frame_vec(frame_ft))) = spot_protein_vec_dist(frame_ft);

        null_fall_lookback = [null_fall_lookback ; null];
        spot_fall_lookback = [spot_fall_lookback ; spot];
    end
end

% perform bootstrap sampling
rise_delta_mat = NaN(n_boots,n_lags+1);
fall_delta_mat = NaN(n_boots,n_lags+1);
high_delta_mat = NaN(n_boots,n_lags+1);
low_delta_mat = NaN(n_boots,n_lags+1);
for n = 1:n_boots
    rise_ids = randsample(1:size(spot_rise_lookback,1),size(spot_rise_lookback,1),true);
    rise_delta_mat(n,:) = nanmean(spot_rise_lookback(rise_ids,:)-null_rise_lookback(rise_ids,:));
    
    fall_ids = randsample(1:size(spot_fall_lookback,1),size(spot_fall_lookback,1),true);
    fall_delta_mat(n,:) = nanmean(spot_fall_lookback(fall_ids,:)-null_fall_lookback(fall_ids,:));
    
    high_ids = randsample(1:size(spot_high_lookback,1),size(spot_high_lookback,1),true);
    high_delta_mat(n,:) = nanmean(spot_high_lookback(fall_ids,:)-null_high_lookback(fall_ids,:));
    
    low_ids = randsample(1:size(spot_low_lookback,1),size(spot_low_lookback,1),true);
    low_delta_mat(n,:) = nanmean(spot_low_lookback(fall_ids,:)-null_low_lookback(fall_ids,:));
end

rise_mean = nanmean(rise_delta_mat);
rise_ste = nanstd(rise_delta_mat);

fall_mean = nanmean(fall_delta_mat);
fall_ste = nanstd(fall_delta_mat);

high_mean = nanmean(high_delta_mat);
high_ste = nanstd(high_delta_mat);

low_mean = nanmean(low_delta_mat);
low_ste = nanstd(low_delta_mat);

high_low = figure;
hold on
e1 = errorbar(t_res*(-n_lags:0),high_mean,high_ste);
e1.CapSize = 0;
e2 = errorbar(t_res*(-n_lags:0),low_mean,low_ste);
e2.CapSize = 0;
xlabel('lag (seconds)')
ylabel([protein_name '-' protein_fluor ' (locus minus control) (au)'])
grid on
legend('high','low')
title([protein_name '-' protein_fluor 'Enrichment Preceding Activity'])
saveas(high_low,[writePath 'high_low_lookback.png'])

rise_fall = figure;
hold on
e1 = errorbar(t_res*(-n_lags:0),rise_mean,rise_ste);
e1.CapSize = 0;
e2 = errorbar(t_res*(-n_lags:0),fall_mean,fall_ste);
e2.CapSize = 0;
xlabel('lag (seconds)')
ylabel([protein_name '-' protein_fluor ' (locus minus control) (au)'])
grid on
legend('rise','fall')
title([protein_name '-' protein_fluor 'Enrichment Preceding Change in Activity'])
saveas(rise_fall,[writePath 'rise_fall_lookback.png'])

%% Look at enrichment at start and end of trace lifetime
latest_start = 8*60;
earliest_end = 25*60;

spot_start_lookback = [];
null_start_lookback = [];

spot_stop_lookback = [];
null_stop_lookback = [];

for i = 1:numel(nucleus_struct_protein)
    % basic info
    ParticleID = nucleus_struct_protein(i).ParticleID; 
    frames = nucleus_struct_protein(i).frames;
    time = nucleus_struct_protein(i).time;
    fluo = nucleus_struct_protein(i).fluo;
    % get timing info
    start_time = time(find(~isnan(fluo),1));
    if isempty(start_time)
        continue
    end
    start_lag = start_time - time(1);
    end_lag = time(end) - time(find(~isnan(fluo),1,'last'));
    snip_frames = nucleus_struct_protein(i).snip_frame_vec;

    if start_lag > 120 && start_time <=latest_start
        lookahead_frames = frames(find(~isnan(fluo),1)):frames(find(~isnan(fluo),1))+n_lags; 
        null = NaN(size(lookahead_frames));
        spot = NaN(size(lookahead_frames));
        frame_ft = ismember(snip_ref_frame_vec,lookahead_frames) & particle_vec==ParticleID;
        
        null(ismember(lookahead_frames,snip_ref_frame_vec(frame_ft))) = null_protein_vec_dist(frame_ft);
        spot(ismember(lookahead_frames,snip_ref_frame_vec(frame_ft))) = spot_protein_vec_dist(frame_ft);
        
        null_start_lookback = [null_start_lookback ; null];
        spot_start_lookback = [spot_start_lookback ; spot];
    end
    if end_lag >= 120
        lookahead_frames = frames(find(~isnan(fluo),1,'last'))-n_lags:frames(find(~isnan(fluo),1,'last')); 
        null = NaN(size(lookahead_frames));
        spot = NaN(size(lookahead_frames));
        frame_ft = ismember(snip_ref_frame_vec,lookahead_frames) & particle_vec==ParticleID;
        
        null(ismember(lookahead_frames,snip_ref_frame_vec(frame_ft))) = null_protein_vec_dist(frame_ft);
        spot(ismember(lookahead_frames,snip_ref_frame_vec(frame_ft))) = spot_protein_vec_dist(frame_ft);
        
        null_stop_lookback = [null_stop_lookback ; null];
        spot_stop_lookback = [spot_stop_lookback ; spot];
    end
end


start_delta_mat = NaN(n_boots,n_lags+1);
stop_delta_mat = NaN(n_boots,n_lags+1);
for n = 1:n_boots
    start_ids = randsample(1:size(spot_start_lookback,1),size(spot_start_lookback,1),true);
    start_delta_mat(n,:) = nanmean(spot_start_lookback(start_ids,:)-null_start_lookback(start_ids,:));
    
    stop_ids = randsample(1:size(spot_stop_lookback,1),size(spot_stop_lookback,1),true);
    stop_delta_mat(n,:) = nanmean(spot_stop_lookback(stop_ids,:)-null_stop_lookback(stop_ids,:));
end

start_mean = nanmean(start_delta_mat);
start_ste = nanstd(start_delta_mat);

stop_mean = nanmean(stop_delta_mat);
stop_ste = nanstd(stop_delta_mat);


start_fig = figure;
hold on
e1 = errorbar(t_res*(0:n_lags),start_mean,start_ste);
e1.CapSize = 0;

xlabel('lead (seconds)')
ylabel([protein_name '-' protein_fluor ' (locus minus control) (au)'])
grid on
title([protein_name '-' protein_fluor 'Enrichment Following Onset of Activity'])
saveas(start_fig,[writePath 'start_lookahead.png'])

stop_fig = figure;
hold on
e1 = errorbar(t_res*(-n_lags:0),stop_mean,stop_ste);
e1.CapSize = 0;

xlabel('lag (seconds)')
ylabel([protein_name '-' protein_fluor ' (locus minus control) (au)'])
grid on
title([protein_name '-' protein_fluor 'Enrichment Preceding Cessation of Activity'])
saveas(stop_fig,[writePath 'stop_lookback.png'])

%%
for i = 1:size(spot_protein_snips,3)
    imagesc(spot_protein_snips(:,:,i))
    colorbar
    pause(1)
end
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
