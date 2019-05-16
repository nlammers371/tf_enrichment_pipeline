% Script to probe relationship between input protein concentration and
% output transcriptional response
clear 
close all
% define ID variables
K = 3;
w = 7;
project = 'Dl-Ven x snaBAC';
dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
dataPath = [dropboxFolder '\ProcessedEnrichmentData\' project '\'];
figPath = [dropboxFolder '\LocalEnrichmentFigures\' project '\hmm_input_output_K' num2str(K) '_w' num2str(w) '\'];
mkdir(figPath)
% sampling parameters
window_size = 15; % number of lags and leads over which to track protein/fluo dynamics
nBoots = 100;
protein_name = 'Dorsal';
gene_name = 'snaBAC';
% load data set
load([dataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'])
Tres = hmm_input_output(1).Tres;
x_axis = Tres*(-window_size:window_size)/60;
pt_sm_kernel = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Examine fluorescence levels in vicinity of Dl peaks, troughs, rises, and
%%% falls 
% specify features and data types to use
feature_cell = {'high','low','rise','fall'};
feature_titles = {'peaks', 'troughs', 'rises', 'falls'};
% data_type_cell = {'pt_sm'};
data_titles = {'protein'};

% set scales for feature identification
pt_vec = [hmm_input_output.spot_protein];
pt_scale = prctile(pt_vec,10) - prctile(pt_vec,5);
% pull time series samples from vicinity of protein events
feature_struct = struct;
ref_vec = -window_size:window_size;%1:(2*window_size+1);
iter = 1;
for i = 1:numel(hmm_input_output)
    %%% core data vectors
    time_vec = hmm_input_output(i).time/60;   
    frame_vec = 1:numel(time_vec);
    fluo_vec = hmm_input_output(i).fluo;
    fluo_swap_vec = hmm_input_output(i).swap_fluo;
    pt_spot = hmm_input_output(i).spot_protein;
    pt_mf = hmm_input_output(i).mf_protein;        
    pt_sm_vec = imgaussfilt(pt_spot-pt_mf,pt_sm_kernel);
    qc_filter = hmm_input_output(i).dt_filter_gap;
    pt_nan_ft = ~isnan(pt_sm_vec);
    %%% find features
    % find fluo peaks and troughs    
    [~,pt_high_ids] = findpeaks(pt_sm_vec,frame_vec,'MinPeakProminence',pt_scale); % NL: eye-balled atm
    [~,pt_low_ids] = findpeaks(-pt_sm_vec,frame_vec,'MinPeakProminence',pt_scale); % NL: eye-balled atm

    % identify changepoints     
    pt_d_vec = sign([0 diff(pt_sm_vec)]);    
    ft_frames = frame_vec(pt_nan_ft);
    ipt_pt = findchangepts(pt_sm_vec(pt_nan_ft),'MinThreshold',3*pt_scale);
    ipt_pt = ft_frames(ipt_pt);
    pt_rise_ids = ipt_pt(pt_d_vec(ipt_pt)==1);
    pt_fall_ids = ipt_pt(pt_d_vec(ipt_pt)==-1);        
        
    % apply filters 
    pt_spot(qc_filter) = NaN;
    fluo_swap_vec(qc_filter) = NaN;
    fluo_vec(qc_filter) = NaN;
    
    % sample proiten and activity traces
    for j = 1:numel(feature_cell)
        feature_str = feature_cell{j};       
        % get variables
        eval(['ids = pt_' feature_str '_ids;'])         
        % initialize arrays
        response_temp = NaN(numel(ids),2*window_size+1);            
        swap_response_temp = NaN(numel(ids),2*window_size+1);            
        spot_protein_temp = NaN(numel(ids),2*window_size+1);            
        time_temp = NaN(1,numel(ids));
        for m = 1:numel(ids)                
            raw_ind = ref_vec+ids(m); 
            ft_vec1 = raw_ind > 0 & raw_ind <= numel(pt_sm_vec);
            ft_vec2 = raw_ind(ft_vec1);
            % record
            response_temp(m,ft_vec1) = fluo_vec(ft_vec2);   
            swap_response_temp(m,ft_vec1) = fluo_swap_vec(ft_vec2);                    
            spot_protein_temp(m,ft_vec1) = pt_sm_vec(ft_vec2);                
            time_temp(m) = time_vec(ids(m));               
        end
        % add to main cell structures
        feature_struct(iter).(['pt_' feature_str '_response']) = response_temp;
        feature_struct(iter).(['pt_' feature_str '_swap_response']) = swap_response_temp;
        feature_struct(iter).(['pt_' feature_str '_spot_protein']) = spot_protein_temp;            
        feature_struct(iter).(['pt_' feature_str '_time']) = time_temp;                    
    end
    iter = iter + 1;
end
tic

%%% Calculate bootstrap estimates of Mean and SE
results_struct = struct;
iter = 1;
for j = 1:numel(feature_cell)
    feature_str = feature_cell{j};    
    % get variables
    response_mat = vertcat(feature_struct.(['pt_' feature_str '_response']));
    swap_response_mat = vertcat(feature_struct.(['pt_' feature_str '_swap_response']));
    spot_protein_mat = vertcat(feature_struct.(['pt_' feature_str '_spot_protein']));       
    % bootstrap variables
    boot_index_vec = 1:size(response_mat,1);
    response_boot_mat = NaN(nBoots,size(response_mat,2));
    swap_response_boot_mat = NaN(nBoots,size(response_mat,2));
    spot_protein_boot_mat = NaN(nBoots,size(response_mat,2));    
    for n = 1:nBoots
        boot_ids = randsample(boot_index_vec,numel(boot_index_vec),true);
        response_boot_mat(n,:) = nanmean(response_mat(boot_ids,:));
        swap_response_boot_mat(n,:) = nanmean(swap_response_mat(boot_ids,:));
        spot_protein_boot_mat(n,:) = nanmean(spot_protein_mat(boot_ids,:));        
    end
    results_struct(iter).ID = ['protein ' feature_titles{j}];
    results_struct(iter).fn = ['protein_' feature_titles{j}];
    results_struct(iter).data_name = 'protein';
    % response
    results_struct(iter).response_mean = nanmean(response_boot_mat);
    results_struct(iter).response_ste = nanstd(response_boot_mat);
    % response
    results_struct(iter).swap_response_mean = nanmean(swap_response_boot_mat);
    results_struct(iter).swap_response_ste = nanstd(swap_response_boot_mat);
    % difference
    results_struct(iter).diff_response_mean = nanmean(response_boot_mat - swap_response_boot_mat);
    results_struct(iter).diff_response_ste = nanstd(response_boot_mat - swap_response_boot_mat);
    % spot protein
    results_struct(iter).spot_protein_mean = nanmean(spot_protein_boot_mat);
    results_struct(iter).spot_protein_ste = nanstd(spot_protein_boot_mat);  

    iter = iter + 1;    
end

%%% Make figures                     
time_axis = Tres*ref_vec / 60;
% Define some colors
yw = [234 194 100]/256; % yellow
bl = [115 143 193]/256; % blue
rd = [213 108 85]/256; % red
gr = [191 213 151]/256; % green
br = [207 178 147]/256; % brown

for i = 1:numel(results_struct)  
    input_output_fig = figure;
    hold on
    yyaxis right
    p1 = plot(time_axis,results_struct(i).spot_protein_mean,'Color',[.6 .6 .6],'LineWidth',1.5);
    ylabel([protein_name ' concentration (au)'])
    ax = gca;
    ax.YColor = 'black';
    set(gca,'ytick',[])
    % generate error range vrctors
    spot_err_top = results_struct(i).diff_response_mean + results_struct(i).diff_response_ste;
    spot_err_bottom = results_struct(i).diff_response_mean - results_struct(i).diff_response_ste;    
    
    yyaxis left
    f1 = fill([time_axis fliplr(time_axis)],[spot_err_top fliplr(spot_err_bottom)],rd);
    f1.FaceAlpha = .2;
    f1.EdgeAlpha = 0;    
    
    p2 = plot(time_axis,results_struct(i).diff_response_mean,'-','Color',rd,'LineWidth',1.5);    
    
    ylabel([gene_name ' activity (au)'])
    ax = gca;
    ax.YColor = 'black';
    
    xlabel('offset (minutes)')       
    title(results_struct(i).ID )
    results_struct(i).project = project;
    legend([p1 p2], 'protein feature','transcriptional response','Location','southwest');%,'fluorescence (control)','fluorescence (trend)')
    grid on
    saveas(input_output_fig,[figPath results_struct(i).fn '_in_out.png'])
end
% save results

save([dataPath 'input_output_results_protein.mat'],'results_struct')