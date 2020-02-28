clear
close all

load("E:\Nick\LivemRNA\Dropbox (Personal)\processedenrichmentdata\Zld-GFP_hbP2P-ParB28mCh\nucleus_struct_protein.mat")

% time indexing variables
time_index = 0:1:30;
set_index = unique([nucleus_struct_protein.setID]);
time_vec = [nucleus_struct_protein.time]/60;
time_sigma = 1;
nBoots = 100;
PixelSize = nucleus_struct_protein(1).PixelSize;
% protein vectors
mf_protein_vec = [nucleus_struct_protein.mf_null_protein_vec]*100;
spot_protein_vec = [nucleus_struct_protein.spot_protein_vec]*100;
null_protein_vec = [nucleus_struct_protein.edge_null_protein_vec]*100;
delta_protein_vec = spot_protein_vec - null_protein_vec;
qc_vec = ~isnan(spot_protein_vec)&~isnan(null_protein_vec);

set_vec = [];
for i = 1:numel(nucleus_struct_protein)
    setID = nucleus_struct_protein(i).setID;
    set_vec = [set_vec repelem(setID, numel(nucleus_struct_protein(i).spot_protein_vec))];
end
% arrays to store time trend results
mf_protein_array = NaN(numel(time_index),numel(set_index),nBoots);
spot_protein_array = NaN(numel(time_index),numel(set_index),nBoots);
spot_protein_array_noqc = NaN(numel(time_index),numel(set_index),nBoots);
null_protein_array = NaN(numel(time_index),numel(set_index),nBoots);
delta_protein_array = NaN(numel(time_index),numel(set_index),nBoots);
time_array = NaN(numel(time_index),numel(set_index),nBoots);

for s = 1:numel(set_index)
    set_ids = find(set_vec==set_index(s)&qc_vec);
    set_ids_noqc = find(set_vec==set_index(s));
    for n = 1:nBoots
        boot_indices = randsample(set_ids,numel(set_ids),true);
        boot_times = time_vec(boot_indices);
        boot_indices_noqc = randsample(set_ids_noqc,numel(set_ids),true);
        boot_times_noqc = time_vec(boot_indices_noqc);
        for t = 1:numel(time_index)
            t_weights = exp(-.5*((boot_times-time_index(t))/time_sigma).^2);
            denom = nansum(t_weights);
            t_weights_noqc = exp(-.5*((boot_times_noqc-time_index(t))/time_sigma).^2);
            denom_noqc = nansum(t_weights_noqc);
            % calculate bootstrap means
            mf_protein_array(t,s,n) = nansum(t_weights.*mf_protein_vec(boot_indices)) / denom;
            spot_protein_array(t,s,n) = nansum(t_weights.*spot_protein_vec(boot_indices)) / denom;
            spot_protein_array_noqc(t,s,n) = nansum(t_weights_noqc.*spot_protein_vec(boot_indices_noqc)) / denom_noqc;
            null_protein_array(t,s,n) = nansum(t_weights.*null_protein_vec(boot_indices)) / denom;
            delta_protein_array(t,s,n) = nansum(t_weights.*delta_protein_vec(boot_indices)) / denom;
            time_array(t,s,n) = nansum(t_weights.*time_vec(boot_indices)) / denom;
        end
    end
end    


mf_protein_mean = nanmean(mf_protein_array,3);
mf_protein_se = nanstd(mf_protein_array,[],3);
spot_protein_mean = nanmean(spot_protein_array,3);
spot_protein_se = nanstd(spot_protein_array,[],3);
spot_protein_mean_noqc = nanmean(spot_protein_array_noqc,3);
spot_protein_se_noqc = nanstd(spot_protein_array_noqc,[],3);
null_protein_mean = nanmean(null_protein_array,3);
null_protein_se = nanstd(null_protein_array,[],3);
delta_protein_mean = nanmean(delta_protein_array,3);
delta_protein_se = nanstd(delta_protein_array,[],3);
time_eff_mean = nanmean(time_array,3);

%%
mf_fig_1 = figure;
hold on
errorbar(time_index,mf_protein_mean,mf_protein_se,'CapSize',0)
errorbar(time_index,spot_protein_mean,mf_protein_se,'CapSize',0)
