function assembleSnips(liveProject)
  disp('Assembling protein snip files into single data structure...')
  % get list of snip files
  snipPathTemp = [liveProject.dataPath '/snip_fragments_temp/'];
  filePattern = fullfile(snipPathTemp, 'snip_data*.mat'); % select all .mat files
  snipFiles = dir(filePattern);
  
  % initialize fields
  snip_data_master = struct;
  
  parfor k = 1 : length(snipFiles)
      baseFileName = snipFiles(k).name;
      snipFileName = [snipFiles(k).folder,filesep, baseFileName];      
      temp = load(snipFileName,'snip_data');
      snip_data = temp.snip_data;
      varNames = fieldnames(snip_data);
      for vn = 1:length(varNames)
        snip_data_master(k).(varNames{vn}) = snip_data.(varNames{vn});
      end
  end
  snip_data = snip_data_master;
  % save
  save([liveProject.dataPath 'snip_data.mat'],'snip_data', '-v7.3');
  % remove snip fragments
  snipCleanup(snipPathTemp)
  
  clear snip_data snip_data_master
  
  disp('Done.')
