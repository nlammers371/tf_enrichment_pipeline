% script to explore trends in opto data
clear
close all

% set basic paths
DataRoot = 'C:\Users\nlamm\Dropbox (Personal)\InductionLogic\';
% DataRoot = 'S:\Nick\Dropbox\InductionLogic\';
project = '20201117_v2';
% project = '20200807_opto_chronic';
DataPath = [DataRoot project filesep];

% get list of data sets
SpotsFile = dir([DataPath '*spots_only.xlsx']);
ProteinFile = dir([DataPath '*protein_only.xlsx']);
SpotsFileWT = dir([DataPath '*WT.xlsx']);

% specify time res
dT = 30;

% load spots dataset (WT)
raw_spots_table_wt = readtable([DataPath SpotsFileWT(1).name]);
raw_spots_array_wt = raw_spots_table_wt{:,:};

% load spots dataset
raw_spots_table = readtable([DataPath SpotsFile(1).name]);
raw_spots_array = raw_spots_table{:,:};
time_vec = 0:dT:dT*size(raw_spots_array,1)-1;

% load protein dataset
raw_protein_table = readtable([DataPath ProteinFile(1).name]);
raw_protein_array = raw_protein_table{:,:};

% It looks like there are 10 spot frames for every protein frame
proteinFrames = 1:10:241;
raw_protein_array_interp = interp1(proteinFrames,raw_protein_array,proteinFrames(1):proteinFrames(end));

%%
pegTime = 150;
window = 45;
beforeIndices = pegTime-window:pegTime-5;
afterIndices = pegTime+5:pegTime+window;

before_fluo = nanmean(raw_spots_array(beforeIndices,:));
after_fluo = nanmean(raw_spots_array(afterIndices,:));

before_fluo_wt = nanmean(raw_spots_array_wt(beforeIndices,:));
after_fluo_wt = nanmean(raw_spots_array_wt(afterIndices,:));
% 
before_protein = nanmean(raw_protein_array_interp(beforeIndices,:));
after_protein = nanmean(raw_protein_array_interp(afterIndices,:));

testFilter = after_fluo == 0 & before_fluo ~= 0; 
testFilterWT = after_fluo_wt == 0 & before_fluo_wt ~= 0; 
sum(testFilter)
sum(testFilterWT)

close all
figure;
scatter(before_protein(~testFilter),after_protein(~testFilter))
hold on
scatter(before_protein(testFilter),after_protein(testFilter),'filled')

figure;
scatter(before_fluo,after_fluo)
