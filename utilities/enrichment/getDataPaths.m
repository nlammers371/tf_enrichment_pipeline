function [InputDataPath, OutputDataPath] = getDataPaths(savioFlag,projectName)

if ~savioFlag
  [liveProject, ~, ~, ~, ~] = headerFunction(projectName);
  InputDataPath = liveProject.dataPath;
  OutputDataPath = InputDataPath;%[InputDataPath filesep 'cpHMM_results' filesep];
else
  currentDir = pwd;
  user_i = strfind(currentDir,'users/');
  slashes = strfind(currentDir,'/');
  indices = slashes(find(slashes>user_i,2));  
  if length(indices)==2
    stop_i = indices(2)-1;
  else
    stop_i = length(currentDir);
  end
  InputDataPath = [currentDir(1:stop_i) filesep 'dat' filesep 'tf_enrichment' filesep projectName filesep];
  OutputDataPath = InputDataPath;%[currentDir(1:stop_i) filesep 'tf_enrichment' filesep 'dat' filesep projectName filesep 'cpHMM_results' filesep];
end
disp(OutputDataPath)
mkdir(OutputDataPath);