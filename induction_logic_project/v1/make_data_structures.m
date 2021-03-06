% script to write csv data to matlab data structures that can be parsed by
% the cpHMM analysis pipeline
clear
close all

% set basic paths
% DataRoot = 'C:\Users\nlamm\Dropbox (Personal)\InductionLogic\';
DataRoot = 'S:\Nick\Dropbox\InductionLogic\';
project  = '20201117';
% project = '20200807_opto_chronic';
DataPath = [DataRoot project filesep];

% get list of data sets
TraceFileList = dir([DataPath '*traces_only.xlsx']);
ProteinFileList = dir([DataPath '*protein_only.xlsx']);
expStrings = {TraceFileList.name};
expStrings = expStrings(~contains(expStrings,'~$'));
expIDs = 1:length(expStrings);

% specify time res
dT = 30;

% define key model architecture parameters 
alpha_frac = 1.302/(1.302 +  1.085/2);


% initialize structure
spot_struct = struct;
if ~isempty(ProteinFileList)
    spot_struct_protein = struct;
end
i_pass = 1;
for e = 1
  raw_table = readtable([DataPath expStrings{e}]);
  raw_array = raw_table{:,:};
  time_vec = 0:dT:dT*size(raw_array,1)-1;
  divFactor = 10^ceil(log10(size(raw_array,2)));
  
  if length(ProteinFileList)>=e
      raw_protein_table = readtable([DataPath ProteinFileList(e).name]);
      raw_protein_array = raw_protein_table{:,:};
      % It looks like there are 10 spot frames for every protein frame
      proteinFrames = 1:10:length(time_vec);
      raw_protein_array_interp = interp1(proteinFrames,raw_protein_array,proteinFrames(1):proteinFrames(end));
  end
    
  for i = 1:size(raw_array,2)
    if any(raw_array(:,i))
%         spot_struct(i_pass).fluoFull = raw_array(:,i)';
%         spot_struct(i_pass).timeFull = time_vec;

        % generate find first and last nonzero entry
        first_i = 1;%find(raw_array(:,i)'~=0,1);
        last_i = size(raw_array,1);%find(raw_array(:,i)'~=0,1,'last');

        spot_struct(i_pass).fluoInterp = raw_array(first_i:last_i,i)';
        spot_struct(i_pass).fluo = spot_struct(i_pass).fluoInterp;
        spot_struct(i_pass).fluo(spot_struct(i_pass).fluo==0) = NaN; %NL: this mimics our formatting
        
        if length(ProteinFileList)>=e
          spot_struct_protein(i_pass).nuclear_protein_vecInterp = raw_protein_array_interp(first_i:last_i,i)';
          spot_struct_protein(i_pass).nuclear_protein_vec = spot_struct_protein(i_pass).nuclear_protein_vecInterp;
          spot_struct_protein(i_pass).nuclear_protein_vec(spot_struct_protein(i_pass).nuclear_protein_vec==0) = NaN; %NL: this mimics our formatting          
        end
        spot_struct(i_pass).timeInterp = time_vec(first_i:last_i);
        spot_struct(i_pass).time = spot_struct(i_pass).timeInterp;

        % ID variables
        spot_struct(i_pass).setID = e;
        spot_struct(i_pass).setName = expStrings{e};
        spot_struct(i_pass).particleID = e + i/divFactor;
        spot_struct(i_pass).alpha_frac = alpha_frac;
        spot_struct(i_pass).Tres = dT;
        spot_struct(i_pass).TraceQCFlag = sum(spot_struct(i_pass).fluoInterp>0)>20;
        spot_struct(i_pass).FrameQCFlags = repelem(spot_struct(i_pass).TraceQCFlag,length(spot_struct(i_pass).fluo));
        
        % add fields to protein structure
        if length(ProteinFileList)>=e
            fNames = fieldnames(spot_struct);
            for f = 1:length(fNames)
                spot_struct_protein(i_pass).(fNames{f}) = spot_struct(i_pass).(fNames{f});
            end
        end
        i_pass = i_pass + 1;
    end    
  end
end

% save
save([DataPath 'spot_struct.mat'],'spot_struct')
if ~isempty(ProteinFileList)
    save([DataPath 'spot_struct_protein.mat'],'spot_struct_protein')
end    