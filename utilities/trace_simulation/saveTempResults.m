function saveTempResults(sweep_step,sweepTemp,tempSavePath)

sweepTemp.sweep_step = sweep_step;
save([tempSavePath 'sweepTemp' sprintf('%08d',sweep_step) '.mat'],'sweepTemp')