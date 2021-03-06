% Pull data from representative and progressively downsample to estimate
% spot path prediction accuracies for different time resolutions
clear
close all

% set basic parameters
min_length = 20;
master_struct = struct;
PrefixCell = {'2019-08-12-A135P_eveWT_30uW_550V','2017-07-05-P190umA_eve_5uW','2020-08-21-WT_new_embryo2'};
tres_vec = [16.8430 10.1890];
% preload relevant structures
for p = 1:length(PrefixCell)
    Prefix = PrefixCell{p};
    liveExperiment = LiveEnrichmentExperiment(Prefix);
    % load particles
    master_struct(p).Particles = getParticles(liveExperiment);
    % get spots structure
    master_struct(p).Spots = getSpots(liveExperiment);
    % save info
    master_struct(p).Prefix = Prefix;
    master_struct(p).liveExperiment = liveExperiment;
end    
%% Set paths to example sets to load
% ~16 second time resolution 

for p = 1:length(PrefixCell)
    Prefix = PrefixCell{p};
%     liveExperiment = LiveEnrichmentExperiment(Prefix);
    % load particles
    Particles = master_struct(p).Particles;%getParticles(liveExperiment);
    % get spots structure
    Spots = master_struct(p).Spots;%getSpots(liveExperiment);
    % for each particle we need to loop through spots and extract the
    % fitted x/y positions, which have greater precision than what is
    % reported in Particles
    useFlags = false(size(Particles));
    iter = 1;
    for i = 1:length(Particles)
        IndexVec = Particles(i).Index;
        FrameVec = Particles(i).Frame;
        
        xFitVec = NaN(size(FrameVec));
        yFitVec = NaN(size(FrameVec));        
        if length(FrameVec) >= min_length
            useFlags(i) = 1;
            for f = 1:length(FrameVec)
                spotFit = Spots(FrameVec(f)).Fits(IndexVec(f));
                bzIndex = ismember(spotFit.z,spotFit.brightestZ);
                xFitVec(f) = spotFit.xFit(bzIndex);
                yFitVec(f) = spotFit.yFit(bzIndex);
            end
        end
        
        Particles(i).xFit = xFitVec;
        Particles(i).yFit = yFitVec;
    end
    master_struct(p).Particles = Particles(useFlags);
end    
    
%% Now let's use these datasets to generate testing data sets
% initialize longform array to store fitting results
spotPathPDOne = [];
spotPathPDTwo = [];

tID = 1;

for p = 2 % let's just progressivley downsample a single project for conistency
    % get basic kalman filter options
    liveExperiment = master_struct(p).liveExperiment;
    
    % initialize kalman options structure
%     kalmanOptions.type = 'ConstantVelocity';
%     nDims = 2;
    kalmanOptions.type = 'ConstantAcceleration';
    nDims = 3;
        
    kalmanOptions.MeasurementNoise = 1/liveExperiment.pixelSize_um; 
    kalmanOptions.MotionNoise = repelem(1,nDims);
    kalmanOptions.InitialError = repelem(kalmanOptions.MeasurementNoise,nDims);
            
    kalmanOptions.measurementFields = {'xPos', 'yPos'};
%     kalmanOptions = determineKalmanOptions(liveExperiment);
    
    % extract particles
    Particles = master_struct(p).Particles;
    
    % initialize struct to hold basic project info
    trackingInfo = struct;
    trackingInfo.nFrames = max([Particles.Frame]);
    trackingInfo.nExtrapFrames = 0;
    
    % get predictions for single gap case  
    for i = 1:length(Particles)
 
        % designate which frames to 
        FrameVec = Particles(i).Frame;
        FrameVecSkip = min(FrameVec):2:max(FrameVec);
        FrameVecFull = min(FrameVec):1:max(FrameVec);
        adjacent_frame_flags = conv([1 1 1],ismember(FrameVecFull,FrameVec),'Full');
        adjacent_frame_flags = adjacent_frame_flags(2:end-1);
        skipFrameFlags = adjacent_frame_flags(ismember(FrameVecFull,FrameVec))==3 & ismember(FrameVec,FrameVecSkip);
%         skipFrameFlags = ~ismember(FrameVec,FrameVecSkip) & skip_options;
        refFrameFlags = ~skipFrameFlags;
        % get kalman predictions
        ParticlesTemp = pathPrediction(Particles(i), trackingInfo, kalmanOptions,refFrameFlags);
        
        % output to use flags
        useFlags = ismember(FrameVecFull,FrameVec);
        
        % generate array to concatenate with master        
        xFitVec = ParticlesTemp.xFit;
        yFitVec = ParticlesTemp.yFit;
        pIDVec = double(repelem(i,length(FrameVec))');
        tIDVec = double(repelem(tID,length(FrameVec))');
        xInfVec = ParticlesTemp.xPosInf;
        yInfVec = ParticlesTemp.yPosInf;
  
        tempArray = [tIDVec pIDVec FrameVec' xFitVec' yFitVec' refFrameFlags' xInfVec(useFlags) yInfVec(useFlags)];
        
        % now generarte "dumb" predictions using (dumber) linear
        % interpolation...
        xLin = NaN(size(xFitVec));
        xLin(refFrameFlags) = xFitVec(refFrameFlags);
        xLin(~refFrameFlags) = interp1(FrameVec(refFrameFlags),xFitVec(refFrameFlags),FrameVec(~refFrameFlags));
        yLin = NaN(size(xFitVec));
        yLin(refFrameFlags) = yFitVec(refFrameFlags);
        yLin(~refFrameFlags) = interp1(FrameVec(refFrameFlags),yFitVec(refFrameFlags),FrameVec(~refFrameFlags));
        
        % ... and  (dumbest) the position of the last measurement
        xLast = NaN(size(xFitVec));
        xLast(refFrameFlags) = xFitVec(refFrameFlags);
        xLast(~refFrameFlags) = interp1(FrameVec(refFrameFlags),xFitVec(refFrameFlags),FrameVec(~refFrameFlags),'previous');
        yLast = NaN(size(xFitVec));
        yLast(refFrameFlags) = yFitVec(refFrameFlags);
        yLast(~refFrameFlags) = interp1(FrameVec(refFrameFlags),yFitVec(refFrameFlags),FrameVec(~refFrameFlags),'previous');
        
        % add to array
        tempArray(:,end+1:end+4) = [xLin' yLin' xLast' yLast'];
        
        if isempty(spotPathPDOne)
            spotPathPDOne = tempArray;
        else
            spotPathPDOne = cat(1,spotPathPDOne,tempArray);
        end
    end
    
    % get predictions for dobule gap case
    for i = 1:length(Particles)
 
        % designate which frames to 
        FrameVec = Particles(i).Frame;
        FrameVecSkip = min(FrameVec):4:max(FrameVec);        
        FrameVecFull = min(FrameVec):1:max(FrameVec);
        
        adjacent_frame_flags = conv([1 1 1 1 1],ismember(FrameVecFull,FrameVec(1:2:end)),'Full');
        adjacent_frame_flags = adjacent_frame_flags(3:end-2);
        skipFrameFlags = adjacent_frame_flags(ismember(FrameVecFull,FrameVec))==3 & ismember(FrameVec,FrameVecSkip);
        refFrameFlags =  bwdist(skipFrameFlags)==2;
        
%         skipFrameFlags = ~ismember(FrameVec,FrameVecSkip) & skip_options;
%         refFrameFlags = ~skipFrameFlags;
        if sum(refFrameFlags) >= 2
            % get kalman predictions
            ParticlesTemp = pathPrediction(Particles(i), trackingInfo, kalmanOptions, refFrameFlags);

            % output to use flags
            useFlags = ismember(FrameVecFull,FrameVec);

            % generate array to concatenate with master        
            xFitVec = ParticlesTemp.xFit;
            yFitVec = ParticlesTemp.yFit;
            pIDVec = double(repelem(i,length(FrameVec))');
            tIDVec = double(repelem(tID,length(FrameVec))');
            xInfVec = ParticlesTemp.xPosInf;
            yInfVec = ParticlesTemp.yPosInf;

            tempArray4 = [tIDVec pIDVec FrameVec' xFitVec' yFitVec' skipFrameFlags' refFrameFlags' xInfVec(useFlags) yInfVec(useFlags)];

            % now generarte "dumb" predictions using (dumber) linear
            % interpolation...
            xLin = NaN(size(xFitVec));
            xLin(refFrameFlags) = xFitVec(refFrameFlags);
            xLin(~refFrameFlags) = interp1(FrameVec(refFrameFlags),xFitVec(refFrameFlags),FrameVec(~refFrameFlags));
            yLin = NaN(size(xFitVec));
            yLin(refFrameFlags) = yFitVec(refFrameFlags);
            yLin(~refFrameFlags) = interp1(FrameVec(refFrameFlags),yFitVec(refFrameFlags),FrameVec(~refFrameFlags));

            % ... and  (dumbest) the position of the last measurement
            xLast = NaN(size(xFitVec));
            xLast(refFrameFlags) = xFitVec(refFrameFlags);
            xLast(~refFrameFlags) = interp1(FrameVec(refFrameFlags),xFitVec(refFrameFlags),FrameVec(~refFrameFlags),'previous');
            yLast = NaN(size(xFitVec));
            yLast(refFrameFlags) = yFitVec(refFrameFlags);
            yLast(~refFrameFlags) = interp1(FrameVec(refFrameFlags),yFitVec(refFrameFlags),FrameVec(~refFrameFlags),'previous');

            % add to array
            tempArray4(:,end+1:end+4) = [xLin' yLin' xLast' yLast'];

            if isempty(spotPathPDTwo)
                spotPathPDTwo = tempArray4;
            else
                spotPathPDTwo = cat(1,spotPathPDTwo,tempArray4);
            end
        end
    end

end    
    
%% Assess average errors for each method
gapFlags = spotPathPDOne(:,6)==0;
kalmanErrors2 = sqrt(sum((spotPathPDOne(gapFlags,4:5)-spotPathPDOne(gapFlags,7:8)).^2,2));
interpErrors2 = sqrt(sum((spotPathPDOne(gapFlags,4:5)-spotPathPDOne(gapFlags,9:10)).^2,2));
prevErrors2 = sqrt(sum((spotPathPDOne(gapFlags,4:5)-spotPathPDOne(gapFlags,11:12)).^2,2));

gapFlags4 = spotPathPDTwo(:,6)==1;
kalmanErrors4 = sqrt(sum((spotPathPDTwo(gapFlags4,4:5)-spotPathPDTwo(gapFlags4,8:9)).^2,2));
interpErrors4 = sqrt(sum((spotPathPDTwo(gapFlags4,4:5)-spotPathPDTwo(gapFlags4,10:11)).^2,2));
prevErrors4 = sqrt(sum((spotPathPDTwo(gapFlags4,4:5)-spotPathPDTwo(gapFlags4,12:13)).^2,2));

close all
err_bins = linspace(0,25);
figure('Position',[100 100 512 512]);
hold on
histogram(interpErrors2,err_bins,'Normalization','probability')
histogram(interpErrors4,err_bins,'Normalization','probability')
% histogram(kalmanErrors2,err_bins,'Normalization','probability')


kalErr = nanmedian(kalmanErrors2)
interpErr = nanmedian(interpErrors2)
prevErr = nanmedian(prevErrors2)

figure('Position',[100 100 512 512]);
scatter(tempArray4(tempArray4(:,7)==1,4),tempArray4(tempArray4(:,7)==1,5))
hold on
scatter(tempArray4(tempArray4(:,6)==1,4),tempArray4(tempArray4(:,6)==1,5),'x')
% plot(tempArray4(:,8),tempArray4(:,9),'-^')
plot(tempArray4(:,10),tempArray4(:,11),'-.')
% plot(tempArray(:,4),tempArray(:,5),'-')

% particles_load_path1 = "S:\Nick\Dropbox\EveMutantResults\2019-08-12-A135P_eveWT_30uW_550V\CompiledParticles_2019-08-12-A135P_eveWT_30uW_550V.mat";
% % ~10 second resolution
% particles_load_path2 = "S:\Nick\Dropbox\eveProject\eve7stripes\2017-07-05-P190umA_eve_5uW\CompiledParticles_2017-07-05-P190umA_eve_5uW.mat";