clear
close all

addpath('../utilities/')
% define basic path variables
DropboxFolder =  'E:\Nick\LivemRNA\Dropbox (Personal)\';
Bcd_GFP_project = 'Bcd-GFP_hbP2P-mCh';
Bcd_Venus_project = 'Bcd-Venus';


[RawResultsRoot, ~, ~] =   header_function(DropboxFolder, Bcd_GFP_project);
DataPath =  [DropboxFolder 'ProcessedEnrichmentData\absolute_calibration\'];
mkdir(DataPath);
FigPath = [DropboxFolder 'LocalEnrichmentFigures/PipelineOutput/absolute_calibration/venus_gfp_xcal/'];
mkdir(FigPath);

% load absolute Bcd concentration data from Gregor 2007
bcd_abs_path = 'E:\Nick\LivemRNA\Dropbox (Personal)\processedenrichmentdata\absolute_calibration\GregorData2007\';
bkg_data = readtable([bcd_abs_path 'Gregor2007Black.csv']);
bkg_data_clean.AP = bkg_data.Var1;
bkg_data_clean.nM = bkg_data.Var2;
bcd_data01 = readtable([bcd_abs_path 'Gregor2007BcdRed.csv']);
bcd_data02 = readtable([bcd_abs_path 'Gregor2007BcdBlue.csv']);
bcd_data_clean.AP = [bcd_data01.Var1 ; bcd_data02.Var1];
bcd_data_clean.nM = [bcd_data01.Var2 ; bcd_data02.Var2];
bcd_data_clean.setID = [repelem(1,numel(bcd_data01.Var1)), repelem(2,numel(bcd_data02.Var1))];

% find sheets
sheet_path = [RawResultsRoot 'DataStatus.xlsx'];
[~,sheet_names]=xlsfinfo(sheet_path);

sheet_index_gfp = find(ismember(sheet_names,Bcd_GFP_project));
sheet_index_venus = find(ismember(sheet_names,Bcd_Venus_project));

% get prefix names
% GFP
[~,~,sheet_cell] = xlsread(sheet_path,sheet_index_gfp);
name_col = sheet_cell(1:33,1); % hard coded for now
ready_ft = contains(name_col,'ReadyForEnrichment');
ready_cols = 1 + find([sheet_cell{ready_ft,2:end}]==1);
sheet_cell = sheet_cell(:,[1 ready_cols]);
% get list of project names
prefix_ft = contains(name_col,'Prefix');
prefix_cell_raw = sheet_cell(prefix_ft,2:end);
prefix_cell_gfp = {};
for i = 1:numel(prefix_cell_raw)
    if ~isempty(prefix_cell_raw{i})
        eval([prefix_cell_raw{i} ';'])
        prefix_cell_gfp = [prefix_cell_gfp{:} {Prefix}];
    end
end

% Venus
[~,~,sheet_cell] = xlsread(sheet_path,sheet_index_venus);
name_col = sheet_cell(1:33,1); % hard coded for now
ready_ft = contains(name_col,'ReadyForEnrichment');
ready_cols = 1 + find([sheet_cell{ready_ft,2:end}]==1);
sheet_cell = sheet_cell(:,[1 ready_cols]);
% get list of project names
prefix_ft = contains(name_col,'Prefix');
prefix_cell_raw = sheet_cell(prefix_ft,2:end);
prefix_cell_venus = {};
for i = 1:numel(prefix_cell_raw)
    if ~isempty(prefix_cell_raw{i})
        eval([prefix_cell_raw{i} ';'])
        prefix_cell_venus = [prefix_cell_venus{:} {Prefix}];
    end
end
% get pixel size (should be consistent across sets
load([RawResultsRoot prefix_cell_venus{1} '/FrameInfo.mat']) 
PixelSize = FrameInfo(1).PixelSize;
IntegrationRadius=2;       %Radius of the integration region in um
IntegrationRadius=floor(IntegrationRadius/FrameInfo(1).PixelSize); %Radius of the integration in pixels
if ~mod(IntegrationRadius,2)
    IntegrationRadius=IntegrationRadius+1;
end
Circle=false(3*IntegrationRadius,3*IntegrationRadius);
Circle=MidpointCircle(Circle,IntegrationRadius,1.5*IntegrationRadius+0.5,...
    1.5*IntegrationRadius+0.5,1);
n_pixels = sum(Circle(:));

% load nuclear fluorescence data

% Venus
time_vec_venus = [];
ap_vec_venus = [];
pt_vec_venus = [];
set_vec_venus = [];

for i = 1:numel(prefix_cell_venus)
    % load
    load([RawResultsRoot prefix_cell_venus{i} '/CompiledNuclei.mat']) 
    % extract data
    nc14_time = ElapsedTime(nc14);
    TimeMat = repmat(ElapsedTime',1,size(AllTracesVector,2))-nc14_time;
    APMat = repmat(AllTracesAP',size(AllTracesVector,1),1);
    % generate vectors
    qc_filter = ~isnan(AllTracesVector)&TimeMat>=0;
    time_vec = TimeMat(qc_filter)';
    ap_vec = APMat(qc_filter)';
    pt_vec = AllTracesVector(qc_filter)';
    % add to master vectors
    time_vec_venus = [time_vec_venus time_vec];
    ap_vec_venus = [ap_vec_venus ap_vec];
    pt_vec_venus = [pt_vec_venus pt_vec/n_pixels];
    set_vec_venus = [set_vec_venus repelem(i,numel(time_vec))];
end

% Bcd
time_vec_gfp = [];
ap_vec_gfp = [];
pt_vec_gfp = [];
set_vec_gfp = [];

for i = 1:numel(prefix_cell_gfp)
    % load
    load([RawResultsRoot prefix_cell_gfp{i} '/CompiledNuclei.mat']) 
    % extract data
    nc14_time = ElapsedTime(nc14);
    TimeMat = repmat(ElapsedTime',1,size(AllTracesVector,2))-nc14_time;
    APMat = repmat(AllTracesAP',size(AllTracesVector,1),1);
    % generate vectors
    qc_filter = ~isnan(AllTracesVector)&TimeMat>=0;
    time_vec = TimeMat(qc_filter)';
    ap_vec = APMat(qc_filter)';
    pt_vec = AllTracesVector(qc_filter)';
    % add to master vectors
    time_vec_gfp = [time_vec_gfp time_vec];
    ap_vec_gfp = [ap_vec_gfp ap_vec];
    pt_vec_gfp = [pt_vec_gfp pt_vec/n_pixels];
    set_vec_gfp = [set_vec_gfp repelem(i,numel(time_vec))];
end

% Make basic plots for each fluor
symbol_cell = {'o','s','x','*'};
close all

% Venus 
bcd_venus_fig = figure;
cmap = flipud(brewermap([],'Spectral'));
colormap(cmap);
hold on
% for i = 1:numel(prefix_cell_venus)
%     sft = set_vec_venus == i;
    scatter(ap_vec_venus,pt_vec_venus,15,time_vec_venus,'o','filled','MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2);
% end
h = colorbar;
xlabel('% embryo length')
ylabel('Bcd-Venus intensity (au)')
ylabel(h,'minutes into nc14')
xlim([.17 .38])
set(gca,'FontSize',14)
saveas(bcd_venus_fig,[FigPath 'bcd_venus_scatter.png'])

% GFP
bcd_gfp_fig = figure;
colormap(cmap);
hold on
% for i = 1:numel(prefix_cell_gfp)
%     sft = set_vec_gfp == i;
    scatter(ap_vec_gfp,pt_vec_gfp,15,time_vec_gfp,'filled','MarkerEdgeAlpha',0,'MarkerFaceAlpha',.2);
% end
h = colorbar;
xlabel('% embryo length')
ylabel('Bcd-GFP intensity (au)')
ylabel(h,'minutes into nc14')
% xlim([.17 .38])
set(gca,'FontSize',14)
saveas(bcd_gfp_fig,[FigPath 'bcd_gfp_scatter.png'])

% Perform simple cross-calibration using average values
ap_res = 1;
time_res = 1;
ap_range = 0:50;
time_range = 5:25;
time_mat = repmat(time_range',1,numel(ap_range));
nBoots = 100;

venus_fluo_array = NaN(numel(time_range),numel(ap_range),nBoots);
venus_ct_mat = NaN(numel(time_range),numel(ap_range));
gfp_fluo_array = NaN(numel(time_range),numel(ap_range),nBoots);
gfp_ct_mat = NaN(numel(time_range),numel(ap_range));

for  a = 1:numel(ap_range)
    for t = 1:numel(time_range)
        % Venus
        v_filter = round(100*ap_vec_venus)==ap_range(a) & round(time_vec_venus)==time_range(t);
        venus_ids = find(v_filter);
        % GFP
        g_filter = round(100*ap_vec_gfp)==ap_range(a) & round(time_vec_gfp)==time_range(t);
        g_ids = find(g_filter);
        for n = 1:nBoots
            if numel(venus_ids) >= 10
                boot_ids_v = randsample(venus_ids,numel(venus_ids),true);
                venus_fluo_array(t,a,n) = nanmean(pt_vec_venus(boot_ids_v));            
            end                
            if numel(g_ids) >= 10      
                boot_ids_g = randsample(g_ids,numel(g_ids),true);
                gfp_fluo_array(t,a,n) = nanmean(pt_vec_gfp(g_filter));      
            end
        end
        venus_ct_mat(t,a) = sum(v_filter);
        gfp_ct_mat(t,a) = sum(g_filter);
    end
end
% mean and se
venus_fluo_mean = nanmean(venus_fluo_array,3);
venus_fluo_ste = nanstd(venus_fluo_array,[],3);

gfp_fluo_mean = nanmean(gfp_fluo_array,3);
gfp_fluo_ste = nanstd(gfp_fluo_array,[],3);

% compare Bcd and Venus intensities
close all
xcal_fig = figure;
colormap(cmap)
scatter(venus_fluo_mean(:),gfp_fluo_mean(:),20,time_mat(:),'filled')
grid on 
h = colorbar;
xlabel('Bcd-Venus intensity (au)')
ylabel('Bcd-GFP intensity (au)')
ylabel(h,'minutes into nc14')
set(gca,'FontSize',14)
axis([0 .4 0 .4])
saveas(xcal_fig,[FigPath 'venus_gfp_xcal_scatter.png'])

%% Compare to Gregor Data
close all
% obtain estimate as function of AP
ap_vec_bcd = 100*bcd_data_clean.AP;
ap_vec_bkg = 100*bkg_data_clean.AP;
bcd_set_vec = bcd_data_clean.setID;
bcd_nM_p = polyfit(ap_vec_bcd,bcd_data_clean.nM,3);
bcd_nM_raw = polyval(bcd_nM_p,ap_range);
bkg_nM_p = polyfit(ap_vec_bkg,bkg_data_clean.nM,3);
bkg_nM_raw = polyval(bkg_nM_p,ap_range);

bcd_nM_vec = bcd_nM_raw - bkg_nM_raw;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%     Venus     %%%%%%%%%%%%%%%%
ven_factor = 2/3; % account for fact that Bcd-Venus transgene drives lower-than-endogenous-levels
bcd_nM_long = linspace(0,max(bcd_nM_vec));
% pull Bcd-Venus data at t=15
bcd_venus_mean_vec = venus_fluo_mean(time_range==15,:);
bcd_venus_ste_vec = venus_fluo_ste(time_range==15,:);
% perform simple linear regression
md_venus = fitlm(bcd_nM_vec, ven_factor*bcd_venus_mean_vec,'Intercept',false);
venus_au_per_nM = md_venus.Coefficients.Estimate(1);
venus_au_offset = 0;%md_venus.Coefficients.Estimate(1);
md_pd_trend_venus = venus_au_offset + venus_au_per_nM*bcd_nM_long;

% make fit figure
venus_cal_fig = figure;
hold on
plot(ven_factor*bcd_nM_long,md_pd_trend_venus,'Color','black');
errorbar(ven_factor*bcd_nM_vec,ven_factor*bcd_venus_mean_vec,bcd_venus_ste_vec,'.','Color','black','CapSize',0);
scatter(ven_factor*bcd_nM_vec,ven_factor*bcd_venus_mean_vec,20,'MarkerFaceColor',cmap(10,:),'MarkerEdgeAlpha',0);
% set axis labels etc
xlabel('[Bcd] (nM)')
ylabel('Bcd-Venus intensity (au per pixel)')
grid on
set(gca,'Fontsize',14)
xlim([0 35])
ylim([0 .35])
saveas(venus_cal_fig,[FigPath 'venus_calibration_scatter.png'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%     GFP     %%%%%%%%%%%%%%%%%%%

% pull Bcd-GFP data at t=15
bcd_gfp_mean_vec = gfp_fluo_mean(time_range==15,:);
bcd_gfp_ste_vec = gfp_fluo_ste(time_range==15,:);

% perform simple linear regression
md_gfp = fitlm(bcd_nM_vec, bcd_gfp_mean_vec);%,'Intercept',false);
gfp_au_per_nM = md_gfp.Coefficients.Estimate(2);
gfp_au_offset = md_gfp.Coefficients.Estimate(1);
md_pd_trend_gfp = gfp_au_offset + gfp_au_per_nM*bcd_nM_long;

% make fit figure
gfp_cal_fig = figure;
hold on
plot(bcd_nM_long,md_pd_trend_gfp,'Color','black');
errorbar(bcd_nM_vec,bcd_gfp_mean_vec,bcd_gfp_ste_vec,'.','Color','black','CapSize',0);
scatter(bcd_nM_vec,bcd_gfp_mean_vec,20,'MarkerFaceColor',cmap(end,:),'MarkerEdgeAlpha',0);
% set axis labels etc
xlabel('[Bcd] (nM)')
ylabel('Bcd-GFP intensity (au per pixel)')
xlim([0 35])
ylim([0 .35])
grid on
set(gca,'Fontsize',14)
saveas(gfp_cal_fig,[FigPath 'gfp_calibration_scatter.png'])

% calculate conversion factor
% kappa = 2.33 * 1.5 / 1.4;
VoxelSize = PixelSize^2 * .5;
micronLiter = 1e-15;
Nmol = 6.022e23;

nM_to_mol = Nmol * 1e-9 * VoxelSize / 1e15;

% replot with N_mol per pixel
venus_cal_fig = figure;
hold on
plot(ven_factor*bcd_nM_long*nM_to_mol,md_pd_trend_venus,'Color','black');
errorbar(ven_factor*bcd_nM_vec*nM_to_mol,ven_factor*bcd_venus_mean_vec,bcd_venus_ste_vec,'.','Color','black','CapSize',0);
scatter(ven_factor*bcd_nM_vec*nM_to_mol,ven_factor*bcd_venus_mean_vec,20,'MarkerFaceColor',cmap(10,:),'MarkerEdgeAlpha',0);
% set axis labels etc
xlabel('Bcd-Venus molecules per pixel')
ylabel('Bcd-Venus intensity (au per pixel)')
grid on
set(gca,'Fontsize',14)
% xlim([0 35])
ylim([0 .35])
saveas(venus_cal_fig,[FigPath 'venus_cal_nmol_scatter.png'])


gfp_cal_fig = figure;
hold on
plot(nM_to_mol*bcd_nM_long,md_pd_trend_gfp,'Color','black');
errorbar(nM_to_mol*bcd_nM_vec,bcd_gfp_mean_vec,bcd_gfp_ste_vec,'.','Color','black','CapSize',0);
scatter(nM_to_mol*bcd_nM_vec,bcd_gfp_mean_vec,20,'MarkerFaceColor',cmap(end,:),'MarkerEdgeAlpha',0);
% set axis labels etc
xlabel('Bcd-GFP molecules per pixel')
ylabel('Bcd-GFP intensity (au per pixel)')
% xlim([0 35])
ylim([0 .35])
grid on
set(gca,'Fontsize',14)
saveas(gfp_cal_fig,[FigPath 'gfp_calnmol_scatter.png'])


gregor_data_fig = figure;
hold on
scatter(100*bcd_data_clean.AP(bcd_set_vec==1),bcd_data_clean.nM(bcd_set_vec==1),30,'MarkerFaceColor',cmap(end,:),'MarkerEdgeAlpha',0);
scatter(100*bcd_data_clean.AP(bcd_set_vec==2),bcd_data_clean.nM(bcd_set_vec==2),30,'MarkerFaceColor',cmap(1,:),'MarkerEdgeAlpha',0);
scatter(100*bkg_data_clean.AP,bkg_data_clean.nM,30,'MarkerFaceColor','black','MarkerEdgeAlpha',0);
% set axis labels etc
xlabel('% embryo length')
ylabel('[Bcd]_{nuc} (nM)')
ylim([0 65])
grid on
set(gca,'Fontsize',14)
saveas(gregor_data_fig,[FigPath 'gregor_data_scatter.png'])

%% save calibration info
% basic parameters
calibration_info = struct;
calibration_info.PixelSize = PixelSize;
calibration_info.zStep = 0.5;
calibration_info.VoxelSize = VoxelSize;
calibration_info.Power = 7;
% cal vales
calibration_info.venus_au_per_nM = venus_au_per_nM;
calibration_info.venus_au_per_molecule = venus_au_per_nM/nM_to_mol;
calibration_info.venus_au_offset = venus_au_offset;

calibration_info.gfp_au_per_nM = gfp_au_per_nM;
calibration_info.gfp_au_per_molecule = gfp_au_per_nM/nM_to_mol;
calibration_info.gfp_au_offset = gfp_au_offset;

save([DataPath 'calibration_info.mat'],'calibration_info')
