clear
close all

% set path to analysis folder
analysisDir = 'E:\LocalEnrichment\Data\RawDynamicsData\2019-11-26\';
% define project roots
MCPonly = 'enosx2_MCPmCherry';
MCPVenus = '2xDl_Venus_MCPmCherry';
% get list of files
MCPonlyFiles = dir([analysisDir '*' MCPonly '*510*']);
MCPVenusFiles = dir([analysisDir '*' MCPVenus '*510*']);

% make figure write path
DropboxFolder =  'E:\Nick\LivemRNA\Dropbox (Personal)\';
FigPath = [DropboxFolder 'LocalEnrichmentFigures\PipelineOutput\mcp_saturation\'];
mkdir(FigPath)

% make data write path
OutPath = [DropboxFolder 'ProcessedEnrichmentData\bleedthrough_controls\'];
mkdir(OutPath)

% get dimensions
nChannels = 3;
fname = MCPVenusFiles(1).name;
tiffList = dir([analysisDir  fname '\*tif']);
% data = bfopen([analysisDir  MCPonlyFiles(1).name '\' MCPonlyFiles(1).name '.lif']);
nSlices = numel(tiffList);
imDims = size(imread([analysisDir  fname '\' tiffList(1).name]));
nSlicesPerStack = nSlices/nChannels;

% load files
MCPVenusStruct = struct;
iter = 1;

for f = 1:numel(MCPVenusFiles)
    fname = MCPVenusFiles(f).name;    
    folder = MCPVenusFiles(f).folder;    
    % extract power level
    pw_index = strfind(fname,'Zoom2') + 6;
    u_indices = strfind(fname,'uW');
    pw = str2double(fname(pw_index:u_indices(1)-1));
    % extract condition info (which lasers on?)
    v_index = strfind(fname,'510') + 4;
    m_index = strfind(fname,'587') + 4;
    c_vec = [strcmp(fname(v_index),'n'),strcmp(fname(m_index),'n')];
    % generate stack to store images
    im_stack = zeros(imDims(1),imDims(2),nSlicesPerStack);
    % read each channel
    for ch = 0:nChannels-1
        % iterate through slices
        for z = 0:nSlicesPerStack-1   
            tiffName = [fname '*z'  sprintf('%02d',z) '_RAW_ch' sprintf('%02d',ch) '.tif'];
            % srach for name 
            im_file = dir([folder '/' fname '/' tiffName]);            
            im = imread([folder '/' fname '/' im_file(1).name]);                
            im_stack(:,:,z+1) = im;        
        end
        im_med = median(im_stack,3);
        im_mean = mean(im_stack,3);
        % record
        MCPVenusStruct(iter).im_stack = im_stack;
        MCPVenusStruct(iter).im_med = im_med;
        MCPVenusStruct(iter).im_mean = im_mean;
        MCPVenusStruct(iter).ch = ch;
        MCPVenusStruct(iter).pw = pw;
        MCPVenusStruct(iter).condition = bi2de(c_vec);
        MCPVenusStruct(iter).fname = fname;
        % generate stack to store images
        im_stack = zeros(imDims(1),imDims(2),nSlicesPerStack);
        % increment
        iter = iter + 1               
    end
end
%% Compare average trends
condition_vec = [MCPVenusStruct.condition];
power_vec = [MCPVenusStruct.pw];
power_index = unique(power_vec);
channel_vec = [MCPVenusStruct.ch];

mcp_array = NaN(imDims(1),numel(power_index));
venus_array = NaN(imDims(1),numel(power_index));

for p = 1:numel(power_index)
    % get index
    v_ft = power_vec == power_index(p) & condition_vec == 1 & channel_vec == 0;
    m_ft = power_vec == power_index(p) & condition_vec == 1 & channel_vec == 1;
    mcp_array(:,p) = nanmean(MCPVenusStruct(m_ft).im_mean,2);
    venus_array(:,p) = nanmean(MCPVenusStruct(v_ft).im_mean,2);
end

% make scatter plot
sc_figure = figure;
hold on
plot(venus_array, mcp_array,'o')
xlabel('Dorsal-Venus intensity')
ylabel('MCP-mCherry intensity')
grid on 
box on
legend('7uW','14uW','21uW','Location','southeast')
set(gca,'Fontsize',14)
ylim([0, .5])
saveas(sc_figure,[FigPath 'Dl-venus vs MCP-mCherry.png'])

%% Use data sets with only Venus laser on to build model
fit_sigma_vec = linspace(.1,100,25);
coeff_array = NaN(numel(fit_sigma_vec),4);
rg_struct = struct;
N_stack = numel(MCPVenusStruct(1).im_stack);
index_vec = 1:N_stack;
n_regress = 1e5;
%%
for f = 1:numel(fit_sigma_vec)
    tic
    ft_sigma = fit_sigma_vec(f);
    bleedthrough_array = NaN(3*n_regress,3);    
    for p = 1:numel(power_index)
        % get index
        v_ft = power_vec == power_index(p) & condition_vec == 1 & channel_vec == 0;
        m_ft = power_vec == power_index(p) & condition_vec == 1 & channel_vec == 1;
        % record
        rg_ids = randsample(index_vec,n_regress,true);
        bleedthrough_array((p-1)*n_regress+1:p*n_regress,2) = power_index(p);
        v_stack = imgaussfilt(MCPVenusStruct(v_ft).im_stack,ft_sigma);
        m_stack1 = imgaussfilt(MCPVenusStruct(m_ft).im_stack,ft_sigma);
        bleedthrough_array((p-1)*n_regress+1:p*n_regress,1) = v_stack(rg_ids);
        bleedthrough_array((p-1)*n_regress+1:p*n_regress,3) = m_stack1(rg_ids);        
    end
    % convert to table array
    bt_table = array2table(bleedthrough_array,'VariableNames',{'venus','power','mcp'});
    bt_table.power = categorical(bt_table.power);
    lm_full = fitlm(bt_table,'mcp~venus+power');
    % record
    coeff_array(f,:) = lm_full.Coefficients.Estimate;
    rg_struct(f).mdl = lm_full;
    toc
end

%% make plot showing avg kernel dependence

sigma_fig = figure;
plot(fit_sigma_vec,coeff_array(:,1),'-o','Color','black')
xlabel('averaging kernel size (pixels)')
ylabel('inferred intercept')
grid on
box on
set(gca,'Fontsize',14)
saveas(sc_figure,[FigPath 'averaging effects.png'])

%% confirm that model fits data well
ft_sigma = 25; % inflection point in plot
bleedthrough_array = NaN(3*n_regress,3);    
for p = 1:numel(power_index)
    % get index
    v_ft = power_vec == power_index(p) & condition_vec == 1 & channel_vec == 0;
    m_ft = power_vec == power_index(p) & condition_vec == 1 & channel_vec == 1;
    % record
    rg_ids = randsample(index_vec,n_regress,true);
    bleedthrough_array((p-1)*n_regress+1:p*n_regress,2) = power_index(p);
    v_stack = imgaussfilt(MCPVenusStruct(v_ft).im_stack,ft_sigma);
    m_stack1 = imgaussfilt(MCPVenusStruct(m_ft).im_stack,ft_sigma);
    bleedthrough_array((p-1)*n_regress+1:p*n_regress,1) = v_stack(rg_ids);
    bleedthrough_array((p-1)*n_regress+1:p*n_regress,3) = m_stack1(rg_ids);        
end
% convert to table array
bt_table = array2table(bleedthrough_array,'VariableNames',{'venus','power','mcp'});
bt_table.power = categorical(bt_table.power);
lm_full = fitlm(bt_table,'mcp~venus+power');


% make prediction
tp = predict(lm_full,bt_table);

% plot
sc_figure = figure;
hold on
plot(venus_array, mcp_array,'d')
scatter(bt_table.venus,tp,5,'black','filled','MarkerEdgeColor','black')
xlabel('Dorsal-Venus intensity')
ylabel('MCP-mCherry intensity')
grid on 
box on
legend('data (7uW)','data (14uW)','data (21uW)','model predictions','Location','southeast')
set(gca,'Fontsize',14)
ylim([0, .5])
saveas(sc_figure,[FigPath 'Dl-venus vs MCP-mCherry_predictions.png'])

% save model 
save([OutPath 'bleedthrough_model.mat'],'lm_full')
%% Attempt to predict both laser scenario from individual lasers

bleedthrough_array = NaN(3*n_regress,4);    
ft_sigma = 25;
for p = 1%:numel(power_index)
    % get index
    v_ft = power_vec == power_index(p) & condition_vec == 3 & channel_vec == 0;
    m_ft1 = power_vec == power_index(p) & condition_vec == 2 & channel_vec == 1;
    m_ft2 = power_vec == power_index(p) & condition_vec == 3 & channel_vec == 1;
    % record
    rg_ids = randsample(index_vec,n_regress,true);
    bleedthrough_array((p-1)*n_regress+1:p*n_regress,2) = power_index(p);
    v_stack = imgaussfilt(MCPVenusStruct(v_ft).im_stack,ft_sigma);
    m_stack1 = imgaussfilt(MCPVenusStruct(m_ft1).im_stack,ft_sigma);
    m_stack2 = imgaussfilt(MCPVenusStruct(m_ft2).im_stack,ft_sigma);
    bleedthrough_array((p-1)*n_regress+1:p*n_regress,1) = v_stack(rg_ids);
    bleedthrough_array((p-1)*n_regress+1:p*n_regress,3) = m_stack1(rg_ids);        
    bleedthrough_array((p-1)*n_regress+1:p*n_regress,4) = m_stack2(rg_ids);        
end
bt_table = array2table(bleedthrough_array(:,1:3),'VariableNames',{'venus','power','mcp'});
bt_table.power = categorical(bt_table.power);
%%
alt_table = bt_table;
alt_table.mcp_diff = bleedthrough_array(:,4) - bleedthrough_array(:,3);
lm_diff = fitlm(alt_table,'mcp_diff~venus')
%% check predictions

venus_index = linspace(prctile(bt_table.venus,1),prctile(bt_table.venus,99));
pd_tbl = struct;
pd_tbl.venus = venus_index';
pd_tbl.power = categorical(repelem(7,100))';
pd_tbl.mcp = NaN(size(pd_tbl.power));
pd_tbl = struct2table(pd_tbl);
venus_bt_pd = predict(lm_full,pd_tbl);


close all
pd_figure = figure;
hold on
scatter(bt_table.mcp, bleedthrough_array(:,end)-bleedthrough_array(:,end-1))
scatter(pd_tbl.venus,venus_bt_pd,5,'black','filled','MarkerEdgeColor','black')
xlabel('Dorsal-Venus intensity')
ylabel('MCP-mCherry intensity')
grid on 
box on
% legend('data (7uW)','data (14uW)','data (21uW)','model predictions','Location','southeast')
set(gca,'Fontsize',14)
% saveas(sc_figure,[FigPath 'Dl-venus vs MCP-mCherry_predictions.png'])

