clear
close all
addpath('../utilities');

% type_name = 'hbv1_sets';
type_name = 'hb_enrichment';
% type_name = 'sna_enrichment';
% specify quantiles to use for analysis
quantiles_in = [15:15:90 98]/100;
tlims_in = [5 25];
nc14_flag_in = true;
% define basic ID variables
DropboxFolder =  'E:\Nick\LivemRNA\Dropbox (Personal)\';
% prefix_cell = {...
%     '2017-09-14-P2P-MCP-NoNLS-mCherry-doubledosage2',...
%     '2017-09-14-P2P-MCP-NoNLS-mCherry-doubledosage3',...
%     '2017-12-13-P2P-MCP-NoNLS-mCherry-doubledosage',...
%     '2017-09-10-P2P-MCP-NoNLS-mCherry-threeovertwodosage',...
%     '2017-09-15-P2P-MCP-NoNLS-mCherry-threeovertwodosage',...
%     '2017-09-17-P2P-MCP-NoNLS-mCherry-threeovertwodosage',...
%     '2017-12-03-P2P-MCP-NoNLS-mCherry',...
%     '2017-12-03-P2P-MCP-NoNLS-mCherry_15',...
%     '2019-10-31-F_F_F_dorsal_synthetics_saturation_check',...
%     '2019-11-07-F_F_F_MCP-NoNLS_x_HbP2Pv1MS2_01',...
%     '2019-11-04-F_F_3_Dlxsna_mcp_check_10umStack',...
%     '2017-09-14-P2P-MCP-NoNLS-mCherry-doubledosage2_rerun',...
%     '2019-11-11-original_Dl_Venus_MCP-NoNLS_x_HbP2Pv1MS2_01'
%     };
prefix_cell = {};
if isempty(prefix_cell)
%     project_cell = {'Dl-Ven_snaBAC-mCh','Bcd-GFP_snaBAC-mCh'};   
    project_cell = {'Dl-Ven_hbP2P-mCh','Bcd-GFP_hbP2P-mCh'};   
    % initialize cell array to store Prefix names
    prefix_cell = {};
    % iterate through list of projects
    for p = 1:numel(project_cell)
        project = project_cell{p};
        % get list of prefixes for each project 
        [RawResultsRoot, ~, FigRoot] = header_function(DropboxFolder, project);
        if p == 1
            FigPath = [FigRoot 'mcp_saturation_analyses/'];
            mkdir(FigPath);
        end
        % find sheet
        sheet_path = [RawResultsRoot 'DataStatus.xlsx'];
        [~,sheet_names]=xlsfinfo(sheet_path);
        sheet_index = find(ismember(sheet_names,project));
        if isempty(sheet_index)
            error('no tab matching "DropboxTab" string found in DataStatus')
        end
        [~,~,sheet_cell] = xlsread(sheet_path,sheet_index);
        name_col = sheet_cell(1:33,1); % hard coded for now
        ready_ft = contains(name_col,'ReadyForEnrichment');
        ready_cols = 1 + find([sheet_cell{ready_ft,2:end}]==1);
        sheet_cell = sheet_cell(:,[1 ready_cols]);
        % get list of project names
        prefix_ft = contains(name_col,'Prefix');
        prefix_cell_raw = sheet_cell(prefix_ft,2:end);    
        for i = 1:numel(prefix_cell_raw)
            if ~isempty(prefix_cell_raw{i})
                eval([prefix_cell_raw{i} ';'])
                prefix_cell = [prefix_cell{:} {Prefix}];
            end
        end
    end
else
    [RawResultsRoot, ~, FigRoot] = header_function(DropboxFolder, 'a');
    FigPath = [FigRoot 'mcp_saturation_analyses/'];
    mkdir(FigPath);
end
% intialize data arrays
spot_mean_array = NaN(numel(quantiles_in),numel(prefix_cell));
spot_se_array = NaN(numel(quantiles_in),numel(prefix_cell));

offset_mean_array = NaN(numel(quantiles_in),numel(prefix_cell));
offset_se_array = NaN(numel(quantiles_in),numel(prefix_cell));

% call function to calculate offset and fluo quantiles for each prefix

disp('calculating quantiles...')
for p = 4%1:numel(prefix_cell)    
    Prefix = prefix_cell{p};
    disp(['analyzing ' Prefix '...'])
    [spot_mean_array(:,p), spot_se_array(:,p), offset_mean_array(:,p), offset_se_array(:,p)...
        , nc14_flag, quantiles, tlims] = ...
                    mcp_saturation_analysis(Prefix,'quantiles',quantiles_in,'nc14_flag',nc14_flag_in,'tlims',tlims_in);
end

% Make Figure
nrows = size(offset_mean_array,1);
inc = 1/nrows;
close all
pct_fig = figure;
hold on
hm_cm = flipud(brewermap(nrows,'Spectral'));
colormap(hm_cm);
plot(offset_mean_array,spot_mean_array,'-o','Color',[0 0 0 .3])
for i = 1:nrows
    scatter(offset_mean_array(i,:),spot_mean_array(i,:),'MarkerFaceColor',hm_cm(i,:),'MarkerEdgeColor','black')
end
xlabel('MCP offset (au)')
ylabel('spot fluorescence (au)')
h = colorbar('YTick',inc/2:inc:1,'YTickLabel',quantiles*100);
ylabel(h,'quantile (%)')
grid on
box on
set(gca,'Fontsize',14)
saveas(pct_fig,[FigPath type_name '_prctile_scatters.png'])

% save data structure with analysis details
detail_struct.nc14_flag = nc14_flag;
detail_struct.tlims = tlims;
detail_struct.prefix_cell = prefix_cell;
detail_struct.quantiles = quantiles;

detail_struct.offset_mean_array = offset_mean_array;
detail_struct.offset_se_array = offset_se_array;
detail_struct.spot_mean_array = spot_mean_array;
detail_struct.spot_se_array = spot_se_array;

save([FigPath type_name '_data.mat'],'detail_struct')