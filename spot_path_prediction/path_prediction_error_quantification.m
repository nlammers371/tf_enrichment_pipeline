% Pull data from representative and progressively downsample to estimate
% spot path prediction accuracies for different time resolutions
clear
close all

% set basic parameters
min_length = 20;
master_struct = struct;
PrefixCell = {'2019-08-12-A135P_eveWT_30uW_550V','2017-07-05-P190umA_eve_5uW'};
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
    
% particles_load_path1 = "S:\Nick\Dropbox\EveMutantResults\2019-08-12-A135P_eveWT_30uW_550V\CompiledParticles_2019-08-12-A135P_eveWT_30uW_550V.mat";
% % ~10 second resolution
% particles_load_path2 = "S:\Nick\Dropbox\eveProject\eve7stripes\2017-07-05-P190umA_eve_5uW\CompiledParticles_2017-07-05-P190umA_eve_5uW.mat";