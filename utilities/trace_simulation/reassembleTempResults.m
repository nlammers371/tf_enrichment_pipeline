function sweepTempFull = reassembleTempResults(tempSavePath,nIterations)

sweepTempFull = struct;
wb = waitbar(0,'Reassembling temporary files...');
for n = 1:nIterations    
    load([tempSavePath 'sweepTemp' sprintf('%08d',n) '.mat'],'sweepTemp')
    fnames = fieldnames(sweepTemp);
    for f = 1:length(fnames)
        sweepTempFull(n).(fnames{f}) = sweepTemp.(fnames{f});
    end
    waitbar(n/nIterations,wb)
end
