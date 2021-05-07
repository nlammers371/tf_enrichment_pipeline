function kalmanOptions = determineKalmanOptions(liveExperiment)

%     if true%strcmp(trackingInfo.kfType,'ConstantVelocity')
      
%     elseif strcmp(trackingInfo.kfType,'ConstantAcceleration')
%       nDims = 3;
%     end
    kalmanOptions.type = 'ConstantVelocity';
    nDims = 2;
    
    %% Set noise parameters
    kalmanOptions.MeasurementNoise = 0.1/liveExperiment.pixelSize_um; 
    kalmanOptions.MotionNoise = repelem(kalmanOptions.MeasurementNoise,nDims);
    kalmanOptions.InitialError = repelem(kalmanOptions.MeasurementNoise,nDims);
            
    kalmanOptions.measurementFields = {'xPos', 'yPos'};
    