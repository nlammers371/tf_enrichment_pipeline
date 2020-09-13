function livemRNAPath = FindLivemRNA
pathCell = regexp(path, pathsep, 'split');
pathIndices = find(contains(pathCell,'LivemRNA'));
% check to see if mRNADynamics folder is on the working path
if isempty(pathIndices)
  % search for livemRNAPath if it was not specified  
  disp('Searching for livemRNA...')
  currDir = pwd;
  slashes = strfind(currDir,'\');
  rootDir = currDir(1:slashes(end-1));
  livemRNAContents = dir([rootDir '**\livemRNA']);
  if isempty(livemRNAContents)
    disp('Could not find LivemRNAFolder. Please select the LivemRNA Directory')
    livemRNAPath = uigetdir;
  else
    livemRNAPath = livemRNAContents(end).folder; 
  end  

  % add mRNADynamics to filepath
  addpath(genpath([livemRNAPath '\mRNADynamics']));
else
  livemRNAPath = pathCell{pathIndices(1)}(1:strfind(pathCell{pathIndices(1)},'LivemRNA')+8);
end