function snipCleanup(snipPath)
    
    filePattern = fullfile(snipPath, 'snip_data*.mat'); % select all .mat files
    snipFiles = dir(filePattern);

    parfor k = 1 : length(snipFiles)
      baseFileName = snipFiles(k).name;
      snipFileName = [snipFiles(k).folder,filesep, baseFileName];     
      delete(snipFileName);
    end