% script to write csv data to matlab data structures that can be parsed by
% the cpHMM analysis pipeline
clear
close all

% set basic paths
DataRoot = 'C:\Users\nlamm\Dropbox (Personal)\InductionLogic\';

% project = '20200807_WT';
project = '20200807_opto_chronic';
DataPath = [DataRoot project filesep];

% get list of data sets
FileList = dir([DataPath '*.xlsx']);
expStrings = {FileList.name};
expStrings = expStrings(~contains(expStrings,'~$'));
expIDs = 1:length(expStrings);

% specify time res
dT = 30;

% define key model architecture parameters 
alpha_frac = 1.302/(1.302 +  1.085/2);

% initialize structure
spot_struct = struct;
i_pass = 1;
for e = expIDs
  raw_table = readtable([DataPath expStrings{e}]);
  raw_array = raw_table{:,:};
  time_vec = 0:dT:dT*size(raw_array,1)-1;
  divFactor = 10^ceil(log10(size(raw_array,2)));
  for i = 1:size(raw_array,2)
    if any(raw_array(:,i))
        spot_struct(i_pass).fluoFull = raw_array(:,i)';
        spot_struct(i_pass).timeFull = time_vec;

        % generate find first and last nonzero entry
        first_i = find(raw_array(:,i)'~=0,1);
        last_i = find(raw_array(:,i)'~=0,1,'last');

        spot_struct(i_pass).fluoInterp = raw_array(first_i:last_i,i)';
        spot_struct(i_pass).fluo = spot_struct(i_pass).fluoInterp;
        spot_struct(i_pass).fluo(spot_struct(i_pass).fluo==0) = NaN; %NL: this mimics our formatting
        spot_struct(i_pass).timeInterp = time_vec(first_i:last_i);
        spot_struct(i_pass).time = spot_struct(i_pass).timeInterp;

        % ID variables
        spot_struct(i_pass).setID = e;
        spot_struct(i_pass).setName = expStrings{e};
        spot_struct(i_pass).particleID = e + i/divFactor;
        spot_struct(i_pass).alpha_frac = alpha_frac;
        spot_struct(i_pass).Tres = dT;
        spot_struct(i_pass).qcFlag = sum(spot_struct(i_pass).fluoFull>0)>20;
        i_pass = i_pass + 1;
    end    
  end
end

% save
save([DataPath 'spot_struct.mat'],'spot_struct')