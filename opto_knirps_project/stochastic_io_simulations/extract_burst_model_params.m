clear
close all

% load cpHMM results
load('S:\Jake\Dropbox\ProcessedEnrichmentData\optokni_eve4+6_MCP-GFP_Homo\cpHMM_results\#good\compiledResults_w7_K3_p0_ap1_t1_f2D.mat')
kon = compiledResults.freq_vec_mean
koff = 1./compiledResults.dur_vec_mean
r_gfp = compiledResults.init_vec_mean
noise_gfp = compiledResults.sigma_vec_mean

% get rough MCP-mCherry cross-calibration 
load('S:\Jake\Dropbox\ProcessedEnrichmentData\optokni_eve4+6_MCP-GFP_Homo\spot_struct.mat')
spot_struct_gfp = spot_struct;
time_vec_gfp = [spot_struct_gfp.time];
fluo_vec_gfp = [spot_struct_gfp.fluo];
ap_vec_gfp = [spot_struct_gfp.APPosParticle];

load('S:\Jake\Dropbox\ProcessedEnrichmentData\optokni_eve4+6_WT\spot_struct.mat')
spot_struct_mch = spot_struct;
time_vec_mch = [spot_struct_mch.time];
fluo_vec_mch = [spot_struct_mch.fluo];
ap_vec_mch = [spot_struct_mch.APPosParticle];

% apply bounds and compare mean spot fluorescence
time_bounds = [15 25]*60; % nuclei must have been around for full extent of this interval 
ap_bounds = [50 60];

fluo_vec_gfp_ft = fluo_vec_gfp(time_vec_gfp>=time_bounds(1)&time_vec_gfp<=time_bounds(2) & ...
                               ap_vec_gfp>=ap_bounds(1)&ap_vec_gfp<=ap_bounds(2));
                             
fluo_vec_mch_ft = fluo_vec_mch(time_vec_mch>=time_bounds(1)&time_vec_mch<=time_bounds(2) & ...
                               ap_vec_mch>=ap_bounds(1)&ap_vec_mch<=ap_bounds(2));
                             
gfp_to_mch_factor = nanmean(fluo_vec_mch_ft)/nanmean(fluo_vec_gfp_ft)
 

r_mch = gfp_to_mch_factor * r_gfp
r_mch_per_step = r_mch/3
noise_mch = gfp_to_mch_factor*noise_gfp