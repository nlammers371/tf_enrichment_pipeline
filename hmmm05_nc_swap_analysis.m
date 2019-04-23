% Script to probe relationship between input protein concentration and
% output transcriptional response
clear 
close all
% define ID variables
K = 3;
w = 6;
project = 'Dl_Venus_snaBAC_MCPmCherry_Leica_Zoom2_7uW14uW';
dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
% dropboxFolder = 'C:\Users\nlamm\Dropbox (Garcia Lab)\';
dataPath = [dropboxFolder '\ProcessedEnrichmentData\' project '\'];
figPath = [dropboxFolder '\LocalEnrichmentFigures\' project '\hmm_input_output_K' num2str(K) '_w' num2str(w) '\'];
mkdir(figPath)

nTraces = 50; % number of individual traces to select for plotting
window_size = 15; % number of lags and leads over which to track protein/fluo dynamics
nBoots = 100;

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

% load data set
load([dataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'])

close all
% ignore possibility of second project for noq 
Tres = hmm_input_output(1).Tres;

% Find nearest mf protein match for each nucleus

% generate array of average protein levels for each nucleus
time_vec = unique([hmm_input_output.time]);
mf_array = NaN(numel(time_vec),numel(hmm_input_output));
start_time_vec = NaN(size(hmm_input_output));
stop_time_vec = NaN(size(hmm_input_output));
for i = 1:numel(hmm_input_output)
    t_vec = hmm_input_output(i).time;
    mf_vec = hmm_input_output(i).mf_protein_all;
    mf_array(ismember(time_vec,t_vec),i) = mf_vec;
    start_time_vec(i) = min(t_vec);
    stop_time_vec(i) = max(t_vec);
end
% now find closest match for each nucleus
for i = 1:numel(hmm_input_output)
    mf_i = mf_array(:,i);       
    dt_mf_vec = nanmean(abs(mf_array-mf_i));    
    dt_mf_vec(start_time_vec>start_time_vec(i)|stop_time_vec<stop_time_vec(i)) = NaN;
    dt_mf_vec(i) = NaN;
    [~, best_ind] = nanmin(dt_mf_vec);
    % record vales 
    time_swap = hmm_input_output(best_ind).time;  
    time_i = hmm_input_output(i).time; 
    % fill
    s_pt = NaN(size(time_i));
    s_pt(ismember(time_i,time_swap)) = hmm_input_output(best_ind).spot_protein(ismember(time_swap,time_i));
    s_pt_all = NaN(size(time_i));
    s_pt_all(ismember(time_i,time_swap)) = hmm_input_output(best_ind).spot_protein_all(ismember(time_swap,time_i));
    mf_pt = NaN(size(time_i));
    mf_pt(ismember(time_i,time_swap)) = hmm_input_output(best_ind).mf_protein(ismember(time_swap,time_i));
    mf_pt_all = NaN(size(time_i));
    mf_pt_all(ismember(time_i,time_swap)) = hmm_input_output(best_ind).mf_protein_all(ismember(time_swap,time_i));
    % record
    hmm_input_output(i).swap_ind = best_ind;
    hmm_input_output(i).swap_spot_protein = s_pt;
    hmm_input_output(i).swap_spot_protein_all = s_pt_all;
    hmm_input_output(i).swap_mf_protein = mf_pt;
    hmm_input_output(i).swap_mf_protein_all = mf_pt_all;
end


% First perform basic cross-correlation analysis
r_xcov_spot_mat = NaN(2*window_size+1,numel(hmm_input_output));
r_xcov_swap_mat = NaN(2*window_size+1,numel(hmm_input_output));

f_xcov_spot_mat = NaN(2*window_size+1,numel(hmm_input_output));
f_xcov_swap_mat = NaN(2*window_size+1,numel(hmm_input_output));
iter = 0;
wt_mat = NaN(2*window_size+1,numel(hmm_input_output));
for i = 1:numel(hmm_input_output)
    r_vec = hmm_input_output(i).r_vec;
    f_vec = hmm_input_output(i).fluo;
    pt_spot = hmm_input_output(i).spot_protein_all-hmm_input_output(i).mf_protein_all;
    pt_swap = hmm_input_output(i).swap_spot_protein_all-hmm_input_output(i).swap_mf_protein_all;    
    % get raw-ish vectors
%     qc_flag = hmm_input_output(i).mcp_qc_flag;
    ft = ~isnan(pt_spot)&~isnan(pt_swap);     
    
    if sum(~isnan(f_vec)) > window_size + 1 % && qc_flag
        iter = iter + 1;
        r_xcov_spot_mat(:,i) = xcov(r_vec(ft),pt_spot(ft),window_size,'coeff');        
        r_xcov_swap_mat(:,i) = xcov(r_vec(ft),pt_swap(ft),window_size,'coeff');
        
        f_xcov_spot_mat(:,i) = xcov(f_vec(ft),pt_spot(ft),window_size,'coeff');        
        f_xcov_swap_mat(:,i) = xcov(f_vec(ft),pt_swap(ft),window_size,'coeff');
        
        wt_vec = (sum(ft)-window_size):sum(ft)-1;
        wt_mat(:,i) = [fliplr(wt_vec) sum(ft) wt_vec]; 
    end
end

r_xcov_spot_mean = nansum(r_xcov_spot_mat.*wt_mat,2) ./ nansum(wt_mat,2);
r_xcov_swap_mean = nansum(r_xcov_swap_mat.*wt_mat,2) ./ nansum(wt_mat,2);

f_xcov_spot_mean = nansum(f_xcov_spot_mat.*wt_mat,2) ./ nansum(wt_mat,2);
f_xcov_swap_mean = nansum(f_xcov_swap_mat.*wt_mat,2) ./ nansum(wt_mat,2);
%%% make xcov figures

% hmm-decoded rate
hmm_xcov = figure;
hold on
x_axis = Tres*(-window_size:window_size)/60;
plot(x_axis,r_xcov_spot_mean)
plot(x_axis,r_xcov_swap_mean)
legend('locus','control (random target)')
xlabel('offset (minutes)')
ylabel('cross-covariance')
title('hmm-decoded activity vs. protein')
grid on
saveas(hmm_xcov,[figPath 'xcov_hmm_swap.png'])

% raw fluorescence
fluo_xcov = figure;
hold on
x_axis = Tres*(-window_size:window_size)/60;
plot(x_axis,f_xcov_spot_mean)
plot(x_axis,f_xcov_swap_mean)
legend('locus','control (random target)')
xlabel('offset (minutes)')
ylabel('cross-covariance')
title('spot fluorescence vs. protein')
grid on
saveas(fluo_xcov,[figPath 'xcov_fluo_swap.png'])

% comparison
comp_xcov = figure;
hold on
yyaxis left 
plot(x_axis,r_xcov_spot_mean)
ylabel('cross-covariance (hmm)')

yyaxis right
plot(x_axis,f_xcov_spot_mean)
ylabel('cross-covariance (fluorescence)')

xlabel('offset (minutes)')
legend('hmm activity','fluorescence')
grid on
saveas(comp_xcov,[figPath 'xcov_comparison_swap.png'])
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Examine protein levels in vicinity of sna peaks, troughs, rises, and
%%% falls as detected in raw fluorescence and hmm channels
% specify features and data types to use
feature_cell = {'high','low','rise','fall'};
feature_titles = {'peaks', 'troughs', 'rises', 'falls'};
data_type_cell = {'fluo','hmm'};
data_titles = {'fluorescence', 'HMM'};
% set scales for feature identification
fluo_scale = prctile([hmm_input_output.fluo],40);
r_scale = prctile(vertcat(hmm_input_output.r_vec),40);
% pull time series samples from vicinity of protein events
feature_struct = struct;
ref_vec = -window_size:window_size;%1:(2*window_size+1);
iter = 1;
for i = 1:numel(hmm_input_output)
    %%% core data vectors
    time_vec = hmm_input_output(i).time/60;   
    frame_vec = 1:numel(time_vec);
    fluo_vec = hmm_input_output(i).fluo;
%     qc_flag = hmm_input_output(i).mcp_qc_flag;
    hmm_vec = hmm_input_output(i).r_vec;
    pt_spot = hmm_input_output(i).spot_protein-hmm_input_output(i).mf_protein;      
    pt_swap = hmm_input_output(i).swap_spot_protein-hmm_input_output(i).swap_mf_protein; 
    nan_ft = ~isnan(pt_spot) & ~isnan(pt_swap);
    pt_spot(~nan_ft) = NaN;
    pt_swap(~nan_ft) = NaN;
%     if ~qc_flag
%         continue
%     end
    %%% find features
    % find fluo peaks and troughs
    [~,fluo_high_ids] = findpeaks(fluo_vec,frame_vec,'MinPeakProminence',fluo_scale); % NL: eye-balled atm
    [~,fluo_low_ids] = findpeaks(-fluo_vec,frame_vec,'MinPeakProminence',fluo_scale); % NL: eye-balled atm
    % find hmm peaks
    [~,hmm_high_ids] = findpeaks(hmm_vec,frame_vec,'MinPeakWidth',3,'MinPeakProminence',r_scale); % NL: eye-balled atm
    [~,hmm_low_ids] = findpeaks(-hmm_vec,frame_vec,'MinPeakWidth',3,'MinPeakProminence',r_scale); % NL: eye-balled atm
    % identify changepoints 
    hmm_d_vec = sign([0 diff(hmm_vec')]);
    fluo_d_vec = sign([0 diff(fluo_vec)]);
    ipt_hmm = findchangepts(hmm_vec,'MinThreshold',r_scale);
    ipt_fluo = findchangepts(fluo_vec,'MinThreshold',fluo_scale);
    
    hmm_rise_ids = ipt_hmm(hmm_d_vec(ipt_hmm)==1);
    hmm_fall_ids = ipt_hmm(hmm_d_vec(ipt_hmm)==-1);
    fluo_rise_ids = ipt_fluo(fluo_d_vec(ipt_fluo)==1);
    fluo_fall_ids = ipt_fluo(fluo_d_vec(ipt_fluo)==-1);
    
    % sample proiten and activity traces
    for j = 1:numel(feature_cell)
        feature_str = feature_cell{j};
        for k = 1:numel(data_type_cell)
            data_str = data_type_cell{k};
            % get variables
            eval(['ids =' data_str '_' feature_str '_ids;'])
            eval(['data_vec = ' data_str '_vec;'])
            
            response_temp = NaN(numel(ids),2*window_size+1);            
            spot_protein_temp = NaN(numel(ids),2*window_size+1);
            swap_protein_temp = NaN(numel(ids),2*window_size+1);       
            time_temp = NaN(1,numel(ids));
            for m = 1:numel(ids)                
                raw_ind = ref_vec+ids(m);
                ft_vec1 = raw_ind > 0 & raw_ind <= numel(data_vec);
                ft_vec2 = raw_ind(ft_vec1);
                % record
                response_temp(m,ft_vec1) = data_vec(ft_vec2);          
                spot_protein_temp(m,ft_vec1) = pt_spot(ft_vec2);
                swap_protein_temp(m,ft_vec1) = pt_swap(ft_vec2);  
                time_temp(m) = time_vec(ids(m));               
            end
            % add to main cell structures
            feature_struct(iter).([data_str '_' feature_str '_response']) = response_temp;
            feature_struct(iter).([data_str '_' feature_str '_spot_protein']) = spot_protein_temp;
            feature_struct(iter).([data_str '_' feature_str '_swap_protein']) = swap_protein_temp;
            feature_struct(iter).([data_str '_' feature_str '_time']) = time_temp;            
        end
    end
    iter = iter + 1;
end
tic

%%% Calculate bootstrap estimates of Mean and SE
results_struct = struct;

iter = 1;
for j = 1:numel(feature_cell)
    feature_str = feature_cell{j};
    for k = 1:numel(data_type_cell)
        data_str = data_type_cell{k};
        % get variables
        response_mat = vertcat(feature_struct.([data_str '_' feature_str '_response']));
        spot_protein_mat = vertcat(feature_struct.([data_str '_' feature_str '_spot_protein']));
        swap_protein_mat = vertcat(feature_struct.([data_str '_' feature_str '_swap_protein']));
        % bootstrap variables
        boot_index_vec = 1:size(response_mat,1);
        response_boot_mat = NaN(nBoots,size(response_mat,2));
        spot_protein_boot_mat = NaN(nBoots,size(response_mat,2));
        swap_protein_boot_mat = NaN(nBoots,size(response_mat,2));
        for n = 1:nBoots
            boot_ids = randsample(boot_index_vec,numel(boot_index_vec),true);
            response_boot_mat(n,:) = nanmean(response_mat(boot_ids,:));
            spot_protein_boot_mat(n,:) = nanmean(spot_protein_mat(boot_ids,:));
            swap_protein_boot_mat(n,:) = nanmean(swap_protein_mat(boot_ids,:));
        end
        results_struct(iter).ID = [data_titles{k} ' ' feature_titles{j}];
        results_struct(iter).fn = [data_titles{k} '_' feature_titles{j}];
        results_struct(iter).data_name = data_titles{k};
        % response
        results_struct(iter).response_mean = nanmean(response_boot_mat)-nanmean(response_boot_mat(:));
        results_struct(iter).response_ste = nanstd(response_boot_mat);
        % spot protein
        results_struct(iter).spot_protein_mean = nanmean(spot_protein_boot_mat)-nanmean(spot_protein_boot_mat(:));
        results_struct(iter).spot_protein_ste = nanstd(spot_protein_boot_mat);
        % swap protein
        results_struct(iter).swap_protein_mean = nanmean(swap_protein_boot_mat)-nanmean(swap_protein_boot_mat(:));
        results_struct(iter).swap_protein_ste = nanstd(swap_protein_boot_mat);
        iter = iter + 1;
    end
end

%%% Make figures                     
time_axis = Tres*ref_vec / 60;
cm = jet(128);
red = cm(110,:);
blue = cm(35,:);

for i = 1:numel(results_struct)  
    input_output_fig = figure;
    hold on
    yyaxis left
    plot(time_axis,results_struct(i).response_mean,'Color',[.6 .6 .6])
    ylabel([gene_name ' activity (' results_struct(i).data_name ')'])
    ax = gca;
    ax.YColor = 'black';
    
    yyaxis right
    errorbar(time_axis,results_struct(i).spot_protein_mean,results_struct(i).spot_protein_ste,...
        'Color',red,'CapSize',0,'LineWidth',1.2)
    plot(time_axis,results_struct(i).swap_protein_mean,'-',...results_struct(i).serial_protein_ste,...
        'Color',blue);%,'CapSize',0)    
    ylabel([protein_name ' concentration (au)'])
    ax = gca;
    ax.YColor = 'black';
    
    xlabel('offset (minutes)')       
    title(results_struct(i).ID )
    results_struct(i).project = project;
    legend('transcriptional response','protein (locus)','protein (control)','Location','southwest');%,'fluorescence (control)','fluorescence (trend)')
    grid on
    saveas(input_output_fig,[figPath results_struct(i).fn '_in_out.png'])
end
% save results

save([dataPath 'input_output_results.mat'],'results_struct')