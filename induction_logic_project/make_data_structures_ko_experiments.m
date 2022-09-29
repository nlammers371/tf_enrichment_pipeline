% script to write csv data to matlab data structures that can be parsed by
% the cpHMM analysis pipeline
clear
close all

% set basic paths
DataRoot = [filesep 'Users' filesep 'nick' filesep 'Dropbox (Personal)' filesep];

if ~exist(DataRoot)
  DataRoot = 'S:\Nick\Dropbox (Personal)\';
end

% set basic paths
ReadRoot = [DataRoot 'InductionLogic' filesep 'raw_data' filesep];
WriteRoot = [DataRoot 'ProcessedEnrichmentData' filesep];
% DataRoot = 'S:\Nick\Dropbox\InductionLogic\';
% project  = '20220701_Oct4_opto';
project = '20220912_KO_experiments';
ReadPath = [ReadRoot project filesep];
WritePath = [WriteRoot project filesep];

% get list of data sets
TraceFileList = dir([ReadPath '*.xlsx']);
if isempty(TraceFileList)
    TraceFileList = dir([ReadPath '*.csv']);
end    
ProteinFileList = dir([ReadPath '*protein_only.xlsx']);
expStrings = {TraceFileList.name};
expStrings = expStrings(~contains(expStrings,'~$'));
expIDs = 1:length(expStrings);

% specify time res
dT = 30;
nz_flag = 1;

% define key model architecture parameters 

ms2_len = 1246;
gene_len_vec = [2398 2398 6671 6671];%[2398  2381];
elongation_rate = 2000 * dT/60;
mem_vec = ceil(gene_len_vec ./ elongation_rate);
alpha_frac_vec = ms2_len ./ gene_len_vec;

gen_id_vec = [1 1 2 2];

for e = 1:length(expStrings)
  i_iter = 1;
  i_pass = 1;
  % initialize structure
  spot_struct = struct;
  if ~isempty(ProteinFileList)
      spot_struct_protein = struct;
  end
  
  % get list of sheets
  f_name = TraceFileList(e).name;
  load_string = [TraceFileList(e).folder filesep f_name];
  dash_list = strfind(f_name,'_');
  gene_name = f_name(1:dash_list(1)-1);
  
  % load sheets
  is_csv_file = strcmp(f_name(end-2:end),'csv');
  if ~is_csv_file 
      for k = 1:numel(sheet_names)
          [~,sheet_names]=xlsfinfo(load_string);
          data_sheets = cell(1,length(sheet_names));
          data_sheets{k} = xlsread(load_string,sheet_names{k});
      end
  else
      data_sheets = cell(1,1);
      sheet_names = cell(1,1);
      table_temp = readtable(load_string);
      data_sheets{1} = table_temp{2:end,:};
      sheet_names{1} = f_name(1:end-3);
  end
  % iterate through sheets
  for k = 1:length(data_sheets)
    
    % extract table
    raw_array = data_sheets{k}(3:end,:);
    dark_pt = data_sheets{k}(1,:);
    light_pt = data_sheets{k}(2,:);
%      = raw_table{:,:};
    time_vec = 0:dT:dT*size(raw_array,1)-1;
    divFactor = 10^ceil(log10(size(raw_array,2)));

    if length(ProteinFileList)>=e
        raw_protein_table = readtable([ReadPath ProteinFileList(e).name]);
        raw_protein_array = raw_protein_table{:,:};
        % It looks like there are 10 spot frames for every protein frame
        proteinFrames = 1:10:length(time_vec);
        raw_protein_array_interp = interp1(proteinFrames,raw_protein_array,proteinFrames(1):proteinFrames(end));
    end

    for i = 1:size(raw_array,2)
      if mean(raw_array(:,i)~=0)>0 || ~nz_flag

          % generate find first and last nonzero entry
          first_i = 1;%find(raw_array(:,i)'~=0,1);
          last_i = size(raw_array,1);%find(raw_array(:,i)'~=0,1,'last');

          spot_struct(i_pass).fluoInterp = raw_array(first_i:last_i,i)';
          spot_struct(i_pass).fluo = spot_struct(i_pass).fluoInterp;
          spot_struct(i_pass).fluo(spot_struct(i_pass).fluo==0) = NaN; %NL: this mimics HG Lab formatting
        
          spot_struct(i_pass).dark_yap = dark_pt(i);
          spot_struct(i_pass).light_yap = light_pt(i);
          
          if length(ProteinFileList)>=e
            spot_struct_protein(i_pass).nuclear_protein_vecInterp = raw_protein_array_interp(first_i:last_i,i)';
            spot_struct_protein(i_pass).nuclear_protein_vec = spot_struct_protein(i_pass).nuclear_protein_vecInterp;
            spot_struct_protein(i_pass).nuclear_protein_vec(spot_struct_protein(i_pass).nuclear_protein_vec==0) = NaN; %NL: this mimics our formatting          
          end
          spot_struct(i_pass).timeInterp = time_vec(first_i:last_i);
          spot_struct(i_pass).time = spot_struct(i_pass).timeInterp;

          % ID variables
          spot_struct(i_pass).setID = i_iter;          
          spot_struct(i_pass).geneID = gen_id_vec(e);
          spot_struct(i_pass).expID = e;
          spot_struct(i_pass).geneName = gene_name;
          spot_struct(i_pass).repName = sheet_names{k};
          if any(raw_array(:,i))
              spot_struct(i_pass).particleID = i_iter + i/divFactor;
          else
              spot_struct(i_pass).particleID = NaN;
          end
          spot_struct(i_pass).alpha_frac = alpha_frac_vec(e);
          spot_struct(i_pass).nSteps = mem_vec(e);
          spot_struct(i_pass).Tres = dT;
          spot_struct(i_pass).KO_flag = contains(sheet_names{k},'KO');
          spot_struct(i_pass).off_flag = all(spot_struct(i_pass).fluoInterp==0);
          spot_struct(i_pass).diff_flag = contains(sheet_names{k},'diff');
          spot_struct(i_pass).TraceQCFlag = sum(spot_struct(i_pass).fluoInterp>2e3)>15;
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
    i_iter = i_iter + 1;
  end
  
  % save
  if nz_flag
      projectPath = [WriteRoot filesep project '_nz_' gene_name  filesep];
  else
      projectPath = [WriteRoot filesep project '_' gene_name  filesep];
  end
  mkdir(projectPath);
  save([projectPath  'spot_struct.mat'],'spot_struct')
  if ~isempty(ProteinFileList)
      save([projectPath 'spot_struct_protein.mat'],'spot_struct_protein')
  end    
end

