function cleanUpmRNADynamics
pathCell = regexp(path, pathsep, 'split');
pathIndices = find(contains(pathCell,'LivemRNA\mRNADynamics'));
% check to see if mRNADynamics folder is on the working path
for i = pathIndices
  rmpath(pathCell{i});
end