% Script to call Jonathan's MCMC function to estimate elongation time for
% cell culture traces. For now I think I'm going to generate a foux
% dynamics-results folder for this project to avoid making changes to
% Jonathan's code as much as possible
clear
close all
addpath(genpath('utilities'))

% get current directory
origDir = pwd;
makeNewSet = 1;

% move to livemRNA folder
livemRNAPath = 'C:\Users\nlamm\projects\LivemRNA\';
addpath(genpath(livemRNAPath));
cd(livemRNAPath)

% load WT data from Kirstin 
Prefix = '2020-08-07-wt_test_data';
mvPrefix = '2020-08-07\wt_test_data';
RawDataPath = 'C:\Users\nlamm\projects\induction_project\dat\20200807\';
wt_raw_table = readtable([RawDataPath 'WT.xlsx']);
wt_raw_array = wt_raw_table{:,:};
dT = 0.5; % in minutes
time_vec = 0:dT:dT*size(wt_raw_array,1);


%%%%%
% generate DynamicsResults entries expected by Jonathan's code
% this includes: 
% CompiledParticles
% nc13
% nc14
% ElapsedTime

% First add entry to MovieDatabase if necessary
[~, ~, defaultDropboxFolder, ~, ...
            ~, ~, movieDatabasePath]= DetermineLocalFolders;
movieDatabase = readtable(movieDatabasePath,'PreserveVariableNames',1);

if ~any(contains(movieDatabase.DataFolder,mvPrefix))
  movieDatabase.DataFolder(end+1) = {mvPrefix};
  movieDatabase.RootFolder(end) = {'Default'};
  movieDatabase.DropboxFolder(end) = {'Default'};
  writetable(movieDatabase,movieDatabasePath);
end
  
  
% now generate dummy folder and datasets as needed
PrefixPath = [defaultDropboxFolder '\' Prefix];
if ~exist(PrefixPath)|| makeNewSet
  mkdir(PrefixPath)
  CompiledParticles = struct;
  for i = 1:5%size(wt_raw_array,2)
    rawTrace = wt_raw_array(:,i)';
    frameVec = find(rawTrace~=0,1):find(rawTrace~=0,1,'last');
    rawTrace = rawTrace(frameVec);
    CompiledParticles(i).Fluo = rawTrace;
    CompiledParticles(i).Frame = frameVec;
    CompiledParticles(i).MeanAP = frameVec;
  end
  nc13 = 1;
  nc14 = 1;
  ElapsedTime = time_vec;
  save([PrefixPath '\CompiledParticles.mat'],'CompiledParticles','nc13','nc14','ElapsedTime')
end

% now let's see what happens when I call the MCMC script...
FitSingleParticleMCMC('Prefix',Prefix)