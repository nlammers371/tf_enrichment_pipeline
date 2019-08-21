% Script to run regressions to test simple concentration-dependent
% threshold model for Kr Repression of Eve
clear 
close all
%%% id variables
project = 'kr_eve2_reporter';
%%% Load regression data
ReadPath = ['../../dat/' project '/'];
regression_table = readtable([ReadPath project '_analysis_data.csv']);
%%% set save path 
FigPath = ['../../fig/' project '/'];
mkdir(FigPath)
%%% Define spatiotemporal region to use for analyes
t_bounds = [15 45]*60;
ap_bounds = [40 50]/100;
%%% Extract useful vectors
ap_vec = regression_table.AP;
t_vec = regression_table.time;
on_vec = regression_table.on;
kr_vec = regression_table.pt_mf;
set_vec = regression_table.setID;
fluo_vec = regression_table.fluo_int;
set_index = unique(set_vec);
%%% define filter vec
dp_filter = ap_vec <=ap_bounds(2) & ap_vec >=ap_bounds(1) & ...
            t_vec <=t_bounds(2) & t_vec >=t_bounds(1)& ...
            ~isnan(on_vec);        
%%% conduct set-specific regression analyses
reg_struct = struct;
for s = 1:numel(set_index)
    setID = set_index(s);
    set_struct = struct;
    Y = categorical(on_vec(dp_filter&set_vec==setID));
    KR = kr_vec(dp_filter&set_vec==setID);
    AP = ap_vec(dp_filter&set_vec==setID);
    T = t_vec(dp_filter&set_vec==setID);    
    % simple KR repression
    [B1,dev1,stats1] = mnrfit(KR,Y);    
    set_struct(1).ID = 'kr_only';
    set_struct(1).beta = B1;
    set_struct(1).dev = dev1;
    set_struct(1).p = stats1.p;
    % KR w/ interaction terms
    [B2,dev2,stats2] = mnrfit([KR KR.*AP KR.*T],Y);
    set_struct(2).ID = 'kr_inter';
    set_struct(2).beta = B2;
    set_struct(2).dev = dev2;
    set_struct(2).p = stats2.p;
    % Only Time and AP
    [B3,dev3,stats3] = mnrfit([AP T AP.*T],Y);
    set_struct(3).ID = 'ap_time';
    set_struct(3).beta = B3;
    set_struct(3).dev = dev3;
    set_struct(3).p = stats3.p;
    % record
    reg_struct(s).setID = setID;
    reg_struct(s).results = set_struct;
end
%%% Make some exploratory plots
%%
% compare Kr slope and bkg values across sets
kr_slope_vec = [];
kr_intercept_vec = [];
for i = 1:numel(reg_struct)
    if ~ismember(set_index(i),[2,3,9])
        kr_intercept_vec = [kr_intercept_vec reg_struct(i).results(1).beta(1)];
        kr_slope_vec = [kr_slope_vec reg_struct(i).results(1).beta(2)];
    end
end

kr_fit_fig = figure;
cm = jet(128);
colormap(cm(1:12:120,:));
scatter(kr_intercept_vec, kr_slope_vec,40,'MarkerFaceColor',cm(30,:),'MarkerEdgeColor','black')
grid on
xlabel('kr offset')
ylabel('kr slope')
ylim([0 max(kr_slope_vec)*1.1])
saveas(kr_fit_fig,[FigPath 'simple_fit_fig.png'])

%% Make heatmaps comparing predicted and actual silencing trends
ap_ref = 100*ap_bounds(1):55;
t_ref = (t_bounds(1)/60:t_bounds(2)/60)*60;
hm_struct = struct;

for s = 1:numel(reg_struct)
    close all
    B = reg_struct(s).results(1).beta;
    pd_array = NaN(numel(t_ref),numel(ap_ref));
    obs_array = NaN(numel(t_ref),numel(ap_ref));
    kr_array = NaN(numel(t_ref),numel(ap_ref));
    % get vectors
    setID = set_index(s);    
    Y = on_vec(dp_filter&set_vec==setID);
    KR = kr_vec(dp_filter&set_vec==setID);
    AP = round(ap_vec(dp_filter&set_vec==setID)*100);
    T = round(t_vec(dp_filter&set_vec==setID)/60)*60;    
    pi = mnrval(B,KR);
    
    for a = 1:numel(ap_ref)
        for t = 1:numel(t_ref)
            pd_array(t,a) = nanmean(pi(AP==ap_ref(a)&T==t_ref(t),2));
            obs_array(t,a) = nanmean(Y(AP==ap_ref(a)&T==t_ref(t)));
            kr_array(t,a) = nanmean(double(KR(AP==ap_ref(a)&T==t_ref(t))));
        end
    end
    hm_struct(s).pd_grid = pd_array;
    hm_struct(s).obs_grid = obs_array;
    hm_struct(s).kr_grid = kr_array;
    hm_fig = figure('Position',[0 0 1024 512]);
    colormap(jet(128));
    subplot(1,3,1);
    title('Kr Concentration')
    imagesc(kr_array);
    subplot(1,3,2);
    title('Predicted Faction On')
    imagesc(pd_array);
    caxis([0,1])
    subplot(1,3,3);
    title('Observed Fraction On')
    imagesc(obs_array);
    caxis([0,1])    
    saveas(hm_fig,[FigPath 'set_' sprintf('%02d',set_index(s)) '.png'])
end
%% data qc checks
AP = round(ap_vec(dp_filter&set_vec==2)*100);
T = round(t_vec(dp_filter&set_vec==2)/60)*60;    
FLUO = round(regression_table.fluo_int(dp_filter&set_vec==2)/60)*60;  
ON = on_vec(set_vec==2);
f_ap = [];
ap_index = unique(AP);
for i = 1:numel(unique(AP))    
    f_ap = [f_ap nanmean(ON(AP==ap_index(i)))];
end