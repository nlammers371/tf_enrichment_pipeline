% Pull data from representative and progressively downsample to estimate
% spot path prediction accuracies for different time resolutions
clear
close all

FigurePath = 'C:\Users\nlamm\Dropbox (Personal)\LocalEnrichmentFigures\PipelineOutput\PathPredictionTests\';
mkdir(FigurePath)

% set basic parameters
min_length = 40;
master_struct = struct;
PrefixCell = {'2017-07-06-A200umP_eve_5uW','2017-07-17-A200umP_eve_5uW','2017-09-11-P60umA_200-50V2_5uW'};
Tres = 10.1890;
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
%% Obtain fitted spot positions for each particle
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
    
%% Combine Particles into a single master set and generate down-sampled 
%%%Traces
downsample_struct = struct;
ds_factor_vec = [4 3 2 1]; % keep every second, every third, or every fourht point

% concatenate particles
ParticlesFull = [master_struct.Particles];
ParticlesFullLong = struct;
iter = 1;
% progressively downsample each particle trace
for p = 1:length(ParticlesFull)
    % Frame-related stuff
    FramesRaw = ParticlesFull(p).Frame;
    FramesFull = min(FramesRaw):max(FramesRaw);
    ObsVec = ismember(FramesFull,FramesRaw);
    QueryPointVec = false(size(FramesFull));
    
    % position stuff
    xFitVec = ParticlesFull(p).xFit;
    yFitVec = ParticlesFull(p).yFit;
    xFitVecFull = NaN(size(FramesFull));
    xFitVecFull(ObsVec) = xFitVec;
    yFitVecFull = NaN(size(FramesFull));
    yFitVecFull(ObsVec) = yFitVec;
    
    % do the downsampling
    for d = 1:length(ds_factor_vec)     
        
        % ds factor
        ds_factor = ds_factor_vec(d);
        
        % identify valid options
        conv_kernel = [1 zeros(1,ds_factor-1) 1 zeros(1,ds_factor-1) 1];
        nb_flag_vec = conv(conv_kernel,ObsVec);
        nb_flag_vec = nb_flag_vec(ds_factor+1:end-ds_factor);
        valid_option_flags = nb_flag_vec==3;
        
        % pick series with appropriate spacing
        skip_flag_vec = false(size(valid_option_flags));
        v_indices = find(valid_option_flags);
        skip_flag_vec(v_indices(1)) = 1;
        ind_last = 1;
        ind = 2;
        while ind <= length(v_indices)        
            if v_indices(ind)-v_indices(ind_last)>=2*ds_factor
                skip_flag_vec(v_indices(ind)) = 1;    
                ind_last = ind;
            end
            ind = ind + 1;
        end
        ref_flag_vec = bwdist(skip_flag_vec)==ds_factor;
        fi = 1;%find(ref_flag_vec,1);
        li = length(ref_flag_vec);%find(ref_flag_vec,1,'last');
        % add the fields
        ParticlesFullLong(iter).OriginalParticle = p;
        ParticlesFullLong(iter).dsFactor = ds_factor;
        ParticlesFullLong(iter).FramesFull = FramesFull(fi:li);
        ParticlesFullLong(iter).SkipFlags = skip_flag_vec(fi:li);
        ParticlesFullLong(iter).RefFlags = ref_flag_vec(fi:li);
        ParticlesFullLong(iter).xPosFull = xFitVecFull(fi:li);
        ParticlesFullLong(iter).yPosFull = yFitVecFull(fi:li);
        ParticlesFullLong(iter).xPos = xFitVecFull(fi:li);
        ParticlesFullLong(iter).xPos(~ref_flag_vec(fi:li)) = NaN;
        ParticlesFullLong(iter).yPos = yFitVecFull;
        ParticlesFullLong(iter).yPos(~ref_flag_vec(fi:li)) = NaN;
        
        % increment
        iter = iter + 1;
    end
end

%% Now let's use these datasets to generate testing data sets
% initialize array to store results
spotPathArray = NaN(length([ParticlesFullLong.xPos]),13);

% initialize kalman options
kalmanOptions.type = 'ConstantAcceleration';
nDims = 3;
kalmanOptions.MeasurementNoise = 1/liveExperiment.pixelSize_um; 
kalmanOptions.MotionNoise = repelem(1,nDims);
kalmanOptions.InitialError = repelem(kalmanOptions.MeasurementNoise,nDims);

kalmanOptions.measurementFields = {'xPos', 'yPos'};

% initialize variables to track position in array
start_i = 1;

for p = 1:length(ParticlesFullLong)
    
    FramesFull = ParticlesFullLong(p).FramesFull';
    RefFlags = ParticlesFullLong(p).RefFlags';
    % get kalman predictions
    ParticlesTemp = pathPrediction(ParticlesFullLong(p), kalmanOptions); 

    % generate array to concatenate with master        
    xFitVec = ParticlesTemp.xPosFull;
    yFitVec = ParticlesTemp.yPosFull;
    pIDVec = double(repelem(i,length(yFitVec))');
    dsIDVec = double(repelem(ParticlesFullLong(p).dsFactor,length(yFitVec))');
    opIDVec = double(repelem(ParticlesFullLong(p).OriginalParticle,length(yFitVec))');
    xInfVec = ParticlesTemp.xPosInf;
    yInfVec = ParticlesTemp.yPosInf;

    tempArray = [dsIDVec pIDVec FramesFull xFitVec' yFitVec' ...
                ParticlesFullLong(p).RefFlags' ParticlesFullLong(p).SkipFlags' xInfVec yInfVec];

    % now generarte "dumb" predictions using (dumber) linear
    % interpolation...
    xLin = NaN(size(xFitVec));
    xLin(RefFlags) = xFitVec(RefFlags);
    xLin(~RefFlags) = interp1(FramesFull(RefFlags),xFitVec(RefFlags),FramesFull(~RefFlags));
    
    yLin = NaN(size(xFitVec));
    yLin(RefFlags) = yFitVec(RefFlags);
    yLin(~RefFlags) = interp1(FramesFull(RefFlags),yFitVec(RefFlags),FramesFull(~RefFlags));

    % ... and  (dumbest) the position of the last measurement
    xLast = NaN(size(xFitVec));
    xLast(RefFlags) = xFitVec(RefFlags);
    xLast(~RefFlags) = interp1(FramesFull(RefFlags),xFitVec(RefFlags),FramesFull(~RefFlags),'previous');
    
    yLast = NaN(size(xFitVec));
    yLast(RefFlags) = yFitVec(RefFlags);
    yLast(~RefFlags) = interp1(FramesFull(RefFlags),yFitVec(RefFlags),FramesFull(~RefFlags),'previous');

    % add to array
    tempArray(:,end+1:end+4) = [xLin' yLin' xLast' yLast'];
      
    last_i = start_i + size(tempArray,1) - 1;
    spotPathArray(start_i:last_i,:) = tempArray;
    
    start_i = last_i + 1;    
end    
    
%% Assess average errors for each method

gapFlags = spotPathArray(:,7)==1;
dsIDVec = spotPathArray(gapFlags,1);
kalmanErrors = sqrt(sum((spotPathArray(gapFlags,4:5)-spotPathArray(gapFlags,8:9)).^2,2));
interpErrors = sqrt(sum((spotPathArray(gapFlags,4:5)-spotPathArray(gapFlags,10:11)).^2,2));
prevErrors = sqrt(sum((spotPathArray(gapFlags,4:5)-spotPathArray(gapFlags,12:13)).^2,2));


% close all
% err_bins = linspace(0,25);
% figure('Position',[100 100 512 512]);
% hold on
% histogram(interpErrors(dsIDVec==3),err_bins,'Normalization','probability')
% histogram(interpErrors(dsIDVec==2),err_bins,'Normalization','probability')
% histogram(interpErrors(dsIDVec==1),err_bins,'Normalization','probability')
% histogram(kalmanErrors2,err_bins,'Normalization','probability')

outlier_flags = interpErrors > 10 | prevErrors > 10 | kalmanErrors > 10;

% kalErr = nanmedian(kalmanErrors)
px_size = liveExperiment.pixelSize_um;
kal_err_mean = NaN(1,length(ds_factor_vec));
interp_err_mean = NaN(1,length(ds_factor_vec));
prev_err_mean = NaN(1,length(ds_factor_vec));

for d = 1:length(ds_factor_vec)
    kal_err_mean(d) = px_size*nanmean(kalmanErrors(dsIDVec==ds_factor_vec(d)));
    interp_err_mean(d) = px_size*nanmean(interpErrors(dsIDVec==ds_factor_vec(d)));
    prev_err_mean(d) = px_size*nanmean(prevErrors(dsIDVec==ds_factor_vec(d)));
end
    
x_axis = Tres*(ds_factor_vec);


err_fig = figure;
cmap1 = brewermap([],'Set2');
hold on

plot(x_axis,prev_err_mean,'Color','k')
s1 = scatter(x_axis,prev_err_mean,'MarkerFaceColor',cmap1(1,:),'MarkerEdgeColor','k');

plot(x_axis,interp_err_mean,'Color','k')
s2 = scatter(x_axis,interp_err_mean,'MarkerFaceColor',cmap1(2,:),'MarkerEdgeColor','k');

plot(x_axis,kal_err_mean,'Color','k')
s3 = scatter(x_axis,kal_err_mean,'MarkerFaceColor',cmap1(3,:),'MarkerEdgeColor','k');

xlabel('time gap (seconds)')
ylabel('prediction error (\mu m)')

grid on
set(gca,'Fontsize',14)
legend([s1 s2 s3],'previous position','linear interpolation','kalman filter','Location','southeast')
xlim([0 50])
ylim([0 0.7])
saveas(err_fig,[FigurePath 'prediction_error_v_time_gap.png'])


close all

pd_ex = figure;
hold on
% plot(tempArray(:,4),tempArray(:,5),'--','Color','k');
s1 = scatter(tempArray(tempArray(:,6)==1,4),tempArray(tempArray(:,6)==1,5));
s2 = scatter(tempArray(tempArray(:,7)==1,4),tempArray(tempArray(:,7)==1,5),'s');
p1 = plot(tempArray(:,8),tempArray(:,9),'-^');
% p2 = plot(tempArray(:,10),tempArray(:,11),'-s');

grid on
set(gca,'Fontsize',14)

legend([s1 s2 p1],'observed spots','missed spots','kalman filter prediction','Location','southeast')
xlabel('x position (pixels)')
ylabel('y position (pixels)')
saveas(pd_ex,[FigurePath 'example_plot.png'])
