% Script to generate figures establishing presence of enrichment bursts at
% start of transcription bursts
clear 
% close all
addpath('utilities')
% set ID variables
targetProject = 'Dl-Ven_snaBAC-mCh_v3';
controlProject = 'Dl-Ven_hbP2P-mCh_v2';
DropboxFolder = 'E:\Nick\LivemRNA\Dropbox (Personal)\';
[~, DataPathTarget, FigureRoot] =   header_function(DropboxFolder, targetProject); 
[~, DataPathControl, ~] =   header_function(DropboxFolder, controlProject); 

FigPath = [FigureRoot '\' targetProject '\lag_analyses\'];
mkdir(FigPath)

% load data
load([DataPathTarget 'hmm_input_output_results.mat'])
target_results_struct = results_struct;
load([DataPathControl 'hmm_input_output_results.mat'])
control_results_struct = results_struct;
clear results_struct;

Tres = 20; % seconds
% extract relevant arrays from target project 
lag_dur_vec_target = target_results_struct.lag_dur_vec;
lead_dur_vec_target = target_results_struct.lead_dur_vec;
hmm_array = target_results_struct.hmm_array;
swap_hmm_array = target_results_struct.swap_hmm_array;
spot_array_dt = target_results_struct.spot_array_dt;
% spot_array_dm = target_results_struct.spot_array_dm;
swap_array_dt = target_results_struct.swap_array_dt;
virtual_array_dt = target_results_struct.virtual_array_dt;
feature_sign_vec_target = target_results_struct.feature_sign_vec;
% extract qc vectors 
target_swap_qc = target_results_struct.swap_qc_vec;
target_virtual_qc = target_results_struct.virtual_qc_vec;
target_set_vec = floor(target_results_struct.particle_id_vec);
% extract arrays from control project
lag_dur_vec_control = control_results_struct.lag_dur_vec;
lead_dur_vec_control = control_results_struct.lead_dur_vec;
biocontrol_array_dt = control_results_struct.spot_array_dt;
biohmm_array = control_results_struct.hmm_array;
feature_sign_vec_control = control_results_struct.feature_sign_vec;

min_pause_len = 6;
min_burst_len = 2;
window_size = floor(size(hmm_array,2)/2);
time_axis = (-window_size:window_size)*Tres/60;
%%

burst_ft_primary = feature_sign_vec_target == 1&lead_dur_vec_target>=min_pause_len&lag_dur_vec_target>min_burst_len&target_virtual_qc;%...
    %& sum(~isnan(spot_array_dt),2)'==31; % filter for rise events
burst_ft_control = feature_sign_vec_control == 1&lead_dur_vec_control>=min_pause_len&lag_dur_vec_control>min_burst_len;% ...
    %& sum(~isnan(biocontrol_array_dt),2)'==31; % filter for rise events


n_lags = 2*numel(time_axis)-1;
xcov_array_target = NaN(sum(burst_ft_primary),n_lags);
xcov_array_virtual = NaN(sum(burst_ft_primary),n_lags);
xcov_array_control = NaN(sum(burst_ft_primary),n_lags);
target_indices = find(burst_ft_primary);
control_indices = find(burst_ft_control);
for t = 1:numel(target_indices)
    % target
    hmm_target = hmm_array(target_indices(t),:);
    hmm_target_v = hmm_target;
    pt_target = spot_array_dt(target_indices(t),:);
    v_target = virtual_array_dt(target_indices(t),:);
    
    nan_vec = isnan(hmm_target)|isnan(pt_target);
    nan_vec_v = isnan(hmm_target)|isnan(v_target);
    
    hmm_target(nan_vec) = mean(hmm_target(~nan_vec));
    pt_target(nan_vec) = mean(pt_target(~nan_vec));
    
    v_target(nan_vec_v) = mean(v_target(~nan_vec_v));
    hmm_target_v(nan_vec_v) = mean(hmm_target_v(~nan_vec_v));
    
    xvec = xcov(hmm_target,pt_target)./xcorr(~nan_vec);
    xvec(isinf(xvec)) = NaN;
    xcov_array_target(t,:) = xvec;        
    
    xvec_v = xcov(hmm_target_v,v_target)./xcorr(~nan_vec_v);
    xvec_v(isinf(xvec_v)) = NaN;
    xcov_array_virtual(t,:) = xvec_v;        
end

% control
for c = 1:numel(control_indices)
    hmm_control = biohmm_array(control_indices(c),:);
    pt_control = biocontrol_array_dt(control_indices(c),:);
    nan_vec = isnan(hmm_control)|isnan(pt_control);
    hmm_control(nan_vec) = mean(hmm_control(~nan_vec));
    pt_control(nan_vec) = mean(pt_control(~nan_vec));
    xvec = xcov(hmm_control,pt_control)./xcorr(~nan_vec);
    xvec(isinf(xvec)) = NaN;
    xcov_array_control(c,:) = xvec;
end

cmap1 = brewermap([],'Set2');
x_axis = [-fliplr(1:30) 0 1:30];
xcov_fig = figure;
hold on
plot(x_axis,nanmean(xcov_array_target),'Color',cmap1(2,:),'LineWidth',1.5);
plot(x_axis,nanmean(xcov_array_control),'Color',cmap1(6,:),'LineWidth',1.5);
plot(x_axis,nanmean(xcov_array_virtual),'Color',cmap1(3,:),'LineWidth',1.5);
xlabel('lags (time steps)')
ylabel('cross-covariance')
legend('target locus','control locus','virtual spot')
set(gca,'Fontsize',14)
grid on
xlim([-20 20])
saveas(xcov_fig,[FigPath 'xcov_plot.png'])
saveas(xcov_fig,[FigPath 'xcov_plot.pdf'])

%%
close all
test = rand(1,100);
c1 = xcov(test(11:100),test(1:90));
plot(c1)
grid on