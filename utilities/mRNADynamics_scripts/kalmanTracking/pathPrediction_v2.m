function ParticlesTemp = pathPrediction_v2(ParticlesTemp, kalmanOptions)
 
    % generate sorted position info
    FramesFull = ParticlesTemp.frames;    
    posData = [ParticlesTemp.xPosHR' ParticlesTemp.yPosHR' ParticlesTemp.zPosHR'];          

    % perform forward-backward kalman filtering
    KFTrack1 = kalmanFilterFwd(posData,kalmanOptions);
    KFTrack1 = kalmanFilterBkd(KFTrack1);        

    % flip the filters around to predict early spot positions in frames
    % prior to detection    
    posDataFlipped = flipud(posData);
    
    % 1:find(framesFull==frameVec(lastInd),1)
    KFTrack2 = kalmanFilterFwd(posDataFlipped,kalmanOptions);
    KFTrack2 = kalmanFilterBkd(KFTrack2);  
        
    % Add inferred position info to structure
    ParticlesTemp.framesFull = FramesFull;
    ParticlesTemp.logL = nanmean([KFTrack1.logL flipud(KFTrack2.logL)],2); 
    ParticlesTemp.logLMean = nanmean(ParticlesTemp.logL);
    ParticlesTemp.logLArray = nanmean(cat(3,KFTrack1.logLArray, flipud(KFTrack2.logLArray)),3); 

    % add predictions
    ParticlesTemp.smoothedPredictions = nanmean(cat(3,KFTrack1.smoothedTrack ,flipud(KFTrack2.smoothedTrack)),3);
    ParticlesTemp.smoothedPredictionsSE = nanmean(cat(3,KFTrack1.smoothedTrackSE, flipud(KFTrack2.smoothedTrackSE)),3);
    ParticlesTemp.xPosInf = ParticlesTemp.smoothedPredictions(:,1);
    ParticlesTemp.yPosInf = ParticlesTemp.smoothedPredictions(:,2);
    ParticlesTemp.zPosInf = ParticlesTemp.smoothedPredictions(:,3);
    ParticlesTemp.xPosSEInf = ParticlesTemp.smoothedPredictionsSE(:,1);
    ParticlesTemp.yPosSEInf = ParticlesTemp.smoothedPredictionsSE(:,2);
    ParticlesTemp.zPosSEInf = ParticlesTemp.smoothedPredictionsSE(:,3);
    
%     % create updated particle vectors
%     ParticlesTemp.xPosParticle = NaN(size(ParticlesTemp.xPosParticleOrig));
%     ParticlesTemp.yPosParticle = NaN(size(ParticlesTemp.yPosParticleOrig));
%     ParticlesTemp.zPosParticle = NaN(size(ParticlesTemp.zPosParticleOrig));
%     
%     ParticlesTemp.xPosParticle(~isnan(ParticlesTemp.xPosParticleOrig)) = ParticlesTemp.xPosInf(~isnan(ParticlesTemp.xPosParticleOrig));
%     ParticlesTemp.yPosParticle(~isnan(ParticlesTemp.yPosParticleOrig)) = ParticlesTemp.yPosInf(~isnan(ParticlesTemp.yPosParticleOrig));
%     ParticlesTemp.zPosParticle(~isnan(ParticlesTemp.zPosParticleOrig)) = ParticlesTemp.zPosInf(~isnan(ParticlesTemp.zPosParticleOrig));
%   
%         % supplement with early point predictions
%         ParticlesTemp.zPosInf = ParticlesTemp.zPosDetrendedInf' + trackingInfo.zPosStage(framesFull);           

    % make filter for convenience
%     ParticlesTemp.obsFrameFilter = ismember(FramesFull,frameVec);
        
