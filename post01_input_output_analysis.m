% Script to probe relationship between input protein concentration and
% output transcriptional response
clear 
close all
% define ID variables
K = 3;
w = 7;
project = 'Dl-Ven x hbP2P';
% project = 'Dl_Venus_hbP2P_MCPmCherry_Zoom2_7uW14uW';
dropboxFolder =  'E:\Nick\Dropbox (Garcia Lab)\';
% dropboxFolder = 'C:\Users\nlamm\Dropbox (Garcia Lab)\';
dataPath = [dropboxFolder '\ProcessedEnrichmentData\' project '\'];
figPath = [dropboxFolder '\LocalEnrichmentFigures\' project '\hmm_input_output_K' num2str(K) '_w' num2str(w) '\'];
mkdir(figPath)

nTraces = 50; % number of individual traces to select for plotting
window_size = 15; % number of lags and leads over which to track protein/fluo dynamics
nBoots = 100;
make_trace_plots = 0;

% % extract protein, gene, fluorophore info
% underscores = strfind(project,'_');
% protein_name = project(1:underscores(1)-1);
% protein_fluor = project(underscores(1)+1:underscores(2)-1);
% gene_name = project(underscores(2)+1:underscores(3)-1);
% if numel(underscores) == 3
%     ind = numel(project);
% else
%     ind = underscores(4)-1;
% end
% gene_fluor = project(underscores(3)+1:end);
protein_name = 'Dorsal';
gene_name = 'hbP2P';

% load data set
load([dataPath 'hmm_input_output_w' num2str(w) '_K' num2str(K) '.mat'])
% Define some colors
yw = [234 194 100]/256; % yellow
bl = [115 143 193]/256; % blue
rd = [213 108 85]/256; % red
% gr = [191 213 151]/256; % green
% br = [207 178 147]/256; % brown
% first make figures to ensure that hmmm results have been properly
% concatenated with protein data
if make_trace_plots
    disp('making single-trace plots...')
    qcPath = [figPath '\qc_' project '\'];
    mkdir(qcPath);
    s_index = 1:numel(hmm_input_output);
    rng(123);
    plot_indices = randsample(s_index,min([nTraces,numel(s_index)]),false);
    for j = 1:-1%numel(plot_indices)
        % MCP channel checks
        mcp_check = hmm_input_output(plot_indices(j)).mcp_check;
        fluo_check = hmm_input_output(plot_indices(j)).fluo_check;
        fluo = hmm_input_output(plot_indices(j)).fluo_hmm;
        time = hmm_input_output(plot_indices(j)).time;
        r_vec = sum(hmm_input_output(plot_indices(j)).r_mat,2);
        % make figure
        qc_fig = figure('Visible','off');
        hold on
        plot(time,fluo / nanmean(fluo))
        plot(time,fluo_check / nanmean(fluo_check));
        plot(time,mcp_check / nanmean(mcp_check))
        plot(time,r_vec / nanmean(r_vec))
        legend('fluo (HMM)', 'fluo (data)','raw mcp','activity state (HMM)')
        xlabel('time')
        ylabel([gene_name ' activity (au)'])
        saveas(qc_fig,[qcPath 'mcp_check_nc_' sprintf('%03d',plot_indices(j)) '.png'])

        % Protein Channel checks
        spot_protein = hmm_input_output(plot_indices(j)).spot_protein;
        spot_protein_3D = hmm_input_output(plot_indices(j)).spot_protein_3D;
        null_protein = hmm_input_output(plot_indices(j)).serial_protein;
        null_protein_3D = hmm_input_output(plot_indices(j)).serial_protein_3D;
        mf_protein = hmm_input_output(plot_indices(j)).mf_protein;        
        % make figure
        qc_fig = figure('Visible','off');
        hold on
        plot(time,spot_protein)
        plot(time,null_protein)
        plot(time,spot_protein_3D)
        plot(time,null_protein_3D)
        plot(time,mf_protein)
        legend('protein (spot)', 'protein (control spot)','protein (mf control)','protein 3D (control spot)','protein 3D (mf control)')
        xlabel('time')
        ylabel([protein_name ' - ' protein_fluor ' (au)'])
        saveas(qc_fig,[qcPath 'protein_check_nc_' sprintf('%03d',plot_indices(j)) '.png'])
    end

    % Make single trace input-output plots
    tracePath = [figPath '\single_trace_' project '\'];
    mkdir(tracePath);
    s_index = 1:numel(hmm_input_output);
    rng(321);
    plot_indices = randsample(s_index,min([nTraces,numel(s_index)]),false);
    for j = 1:numel(plot_indices)
        % MCP channel checks
        time = hmm_input_output(plot_indices(j)).time;
        fluo = hmm_input_output(plot_indices(j)).fluo;
        spot_protein = hmm_input_output(plot_indices(j)).spot_protein_3D;
        serial_protein = hmm_input_output(plot_indices(j)).serial_protein_3D;
        mf_protein = hmm_input_output(plot_indices(j)).mf_protein;
        % make figure
        trace_fig = figure('Visible','off');
        hold on

        yyaxis left
        area(time,fluo,'FaceColor',[.7 .7 .7])
        ylabel([gene_name ' activity (au)'])

        yyaxis right
        plot(time,spot_protein,'-','Color',rd,'LineWidth',1.5);
        plot(time,serial_protein,'-','Color',bl,'LineWidth',1.5);
        plot(time,mf_protein,'-','Color',yw,'LineWidth',1.5);
        ylabel(['absolute ' protein_name ' enrichment(au)'])

        ax = gca;
        ax.YAxis(1).Color = 'black';
        ax.YAxis(2).Color = 'black';

%         legend('transcriptional activity', 'protein (active locus)', 'protein (virtual spot)','protein (nucleus average)')
        xlabel('time')
%         error('afs')
        ylabel([gene_name ' activity (au)'])
        saveas(trace_fig,[tracePath 'input_output_nc' sprintf('%03d',plot_indices(j)) '.png'])               
    end
end

close all
Tres = hmm_input_output(1).Tres;

% First perform basic cross-correlation analysis
r_xcov_spot_mat = NaN(2*window_size+1,numel(hmm_input_output));
r_xcov_null_mat = NaN(2*window_size+1,numel(hmm_input_output));

f_xcov_spot_mat = NaN(2*window_size+1,numel(hmm_input_output));
f_xcov_null_mat = NaN(2*window_size+1,numel(hmm_input_output));
iter = 0;
wt_mat = NaN(2*window_size+1,numel(hmm_input_output));
for i = 1:numel(hmm_input_output)
    r_vec = hmm_input_output(i).r_vec;
    f_vec = hmm_input_output(i).fluo;
    pt_spot = hmm_input_output(i).spot_protein-hmm_input_output(i).mf_protein;
    pt_serial = hmm_input_output(i).serial_protein-hmm_input_output(i).mf_protein;
    % get raw-ish vectors    
    ft = ~isnan(pt_spot)&~isnan(pt_serial);     
    
    if numel(f_vec) >= window_size + 1 && sum(ft) / numel(f_vec) > .8
        iter = iter + 1;        
        r_xcov_spot_mat(:,i) = xcov(r_vec(ft),pt_spot(ft),window_size,'coeff');
        r_xcov_null_mat(:,i) = xcov(r_vec(ft),pt_serial(ft),window_size,'coeff');
        
        f_xcov_spot_mat(:,i) = xcov(f_vec(ft),pt_spot(ft),window_size,'coeff');
        f_xcov_null_mat(:,i) = xcov(f_vec(ft),pt_serial(ft),window_size,'coeff');
        
        wt_vec = (sum(ft)-window_size):sum(ft)-1;
        wt_mat(:,i) = [fliplr(wt_vec) sum(ft) wt_vec]; 
    end
end

r_xcov_spot_mean = nansum(r_xcov_spot_mat.*wt_mat,2) ./ nansum(wt_mat,2);
r_xcov_null_mean = nansum(r_xcov_null_mat.*wt_mat,2) ./ nansum(wt_mat,2);

f_xcov_spot_mean = nansum(f_xcov_spot_mat.*wt_mat,2) ./ nansum(wt_mat,2);
f_xcov_null_mean = nansum(f_xcov_null_mat.*wt_mat,2) ./ nansum(wt_mat,2);
%%% make xcov figures

% hmm-decoded rate
hmm_xcov = figure;
hold on
x_axis = Tres*(-window_size:window_size)/60;
plot(x_axis,r_xcov_spot_mean)
plot(x_axis,r_xcov_null_mean)
legend('locus','control (virtual spot)')
xlabel('offset (minutes)')
ylabel('cross-covariance')
title('hmm-decoded activity vs. protein')
grid on
saveas(hmm_xcov,[figPath 'xcov_hmm.png'])

% raw fluorescence
fluo_xcov = figure;
hold on
x_axis = Tres*(-window_size:window_size)/60;
plot(x_axis,f_xcov_spot_mean)
plot(x_axis,f_xcov_null_mean)
legend('locus','control (virtual spot)')
xlabel('offset (minutes)')
ylabel('cross-covariance')
title('spot fluorescence vs. protein')
grid on
saveas(fluo_xcov,[figPath 'xcov_fluo.png'])

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
saveas(comp_xcov,[figPath 'xcov_comparison.png'])
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
r_scale = prctile([hmm_input_output.r_vec],40);
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
    pt_serial = hmm_input_output(i).serial_protein-hmm_input_output(i).mf_protein; 
    
    %%% find features
    % find fluo peaks and troughs
    [~,fluo_high_ids] = findpeaks(fluo_vec,frame_vec,'MinPeakProminence',fluo_scale); % NL: eye-balled atm
    [~,fluo_low_ids] = findpeaks(-fluo_vec,frame_vec,'MinPeakProminence',fluo_scale); % NL: eye-balled atm
    % find hmm peaks
    [~,hmm_high_ids] = findpeaks(hmm_vec,frame_vec,'MinPeakWidth',3,'MinPeakProminence',r_scale); % NL: eye-balled atm
    [~,hmm_low_ids] = findpeaks(-hmm_vec,frame_vec,'MinPeakWidth',3,'MinPeakProminence',r_scale); % NL: eye-balled atm
    % identify changepoints 
    hmm_d_vec = sign([0 diff(hmm_vec)]);
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
            serial_protein_temp = NaN(numel(ids),2*window_size+1);       
            time_temp = NaN(1,numel(ids));
            for m = 1:numel(ids)                
                raw_ind = ref_vec+ids(m);
                ft_vec1 = raw_ind > 0 & raw_ind <= numel(data_vec);
                ft_vec2 = raw_ind(ft_vec1);
                % record
                response_temp(m,ft_vec1) = data_vec(ft_vec2);          
                spot_protein_temp(m,ft_vec1) = pt_spot(ft_vec2);
                serial_protein_temp(m,ft_vec1) = pt_serial(ft_vec2);  
                time_temp(m) = time_vec(ids(m));               
            end
            % add to main cell structures
            feature_struct(iter).([data_str '_' feature_str '_response']) = response_temp;
            feature_struct(iter).([data_str '_' feature_str '_spot_protein']) = spot_protein_temp;
            feature_struct(iter).([data_str '_' feature_str '_serial_protein']) = serial_protein_temp;
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
        serial_protein_mat = vertcat(feature_struct.([data_str '_' feature_str '_serial_protein']));
        % bootstrap variables
        boot_index_vec = 1:size(response_mat,1);
        response_boot_mat = NaN(nBoots,size(response_mat,2));
        spot_protein_boot_mat = NaN(nBoots,size(response_mat,2));
        serial_protein_boot_mat = NaN(nBoots,size(response_mat,2));
        for n = 1:nBoots
            boot_ids = randsample(boot_index_vec,numel(boot_index_vec),true);
            response_boot_mat(n,:) = nanmean(response_mat(boot_ids,:));
            spot_protein_boot_mat(n,:) = nanmean(spot_protein_mat(boot_ids,:));
            serial_protein_boot_mat(n,:) = nanmean(serial_protein_mat(boot_ids,:));
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
        % serial protein
        results_struct(iter).serial_protein_mean = nanmean(serial_protein_boot_mat)-nanmean(serial_protein_boot_mat(:));
        results_struct(iter).serial_protein_ste = nanstd(serial_protein_boot_mat);
        iter = iter + 1;
    end
end
%%
%%% Make figures                     
time_axis = Tres*ref_vec / 60;
% Define some colors
yw = [234 194 100]/256; % yellow
bl = [115 143 193]/256; % blue
rd = [213 108 85]/256; % red
gr = [191 213 151]/256; % green
br = [207 178 147]/256; % brown

if strcmp(project,'Dl_Venus_hbP2P_MCPmCherry_Zoom2_7uW14uW')
    target_color = yw;
else
    target_color = rd;
end

for i = 1:numel(results_struct)  
    input_output_fig = figure;
    hold on
    yyaxis left
    p1 = plot(time_axis,results_struct(i).response_mean,'Color',[.6 .6 .6],'LineWidth',1.5);
    ylabel([gene_name ' activity (' results_struct(i).data_name ')'])
    ax = gca;
    ax.YColor = 'black';
    % generate error range vrctors
    spot_err_top = results_struct(i).spot_protein_mean + results_struct(i).spot_protein_ste;
    spot_err_bottom = results_struct(i).spot_protein_mean - results_struct(i).spot_protein_ste;
    null_err_top = results_struct(i).serial_protein_mean + results_struct(i).serial_protein_ste;
    null_err_bottom = results_struct(i).serial_protein_mean - results_struct(i).serial_protein_ste;
    
    yyaxis right
    f1 = fill([time_axis fliplr(time_axis)],[spot_err_top fliplr(spot_err_bottom)],target_color);
    f1.FaceAlpha = .2;
    f1.EdgeAlpha = 0;
    f2 = fill([time_axis fliplr(time_axis)],[null_err_top fliplr(null_err_bottom)],bl);
    f2.FaceAlpha = .2;
    f2.EdgeAlpha = 0;
    
    p2 = plot(time_axis,results_struct(i).spot_protein_mean,'-','Color',target_color,'LineWidth',1.5);
    p3 = plot(time_axis,results_struct(i).serial_protein_mean,'-','LineWidth',1.5,...
        'Color',bl);%,'CapSize',0)    
    ylabel([protein_name ' enrichment (au)'])
    ax = gca;
    ax.YColor = 'black';
    
    xlabel('offset (minutes)')       
    title(results_struct(i).ID )
    results_struct(i).project = project;
    legend([p1 p2 p3], 'transcriptional response','protein (locus)','protein (control)','Location','southwest');%,'fluorescence (control)','fluorescence (trend)')
    grid on
    saveas(input_output_fig,[figPath results_struct(i).fn '_in_out.png'])
end
% save results

save([dataPath 'input_output_results.mat'],'results_struct')