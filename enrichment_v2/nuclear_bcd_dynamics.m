clear
close all

projectName = 'Bcd-GFP_hbMS2-mCh_Airy_fast_int';

liveProject = LiveEnrichmentProject(projectName);
resultsRoot = [liveProject.dataPath filesep];
FigurePath = [liveProject.figurePath filesep 'nucleus_movies'];
RefPath = [liveProject.dataPath 'refFrames' filesep];

load([resultsRoot 'nucleusMovieData.mat'],'movieData') 

%% first, we need to reconfigure stuff to be in particle-centric

nucleus_snip_radius_um = 2.25; % in microns
PixelSize = liveProject.includedExperiments{1}.pixelSize_um;
nucleus_snip_radius_px = round(nucleus_snip_radius_um/PixelSize);
snip_size = 2*nucleus_snip_radius_px+1;
int_radius = 0.25 / PixelSize;
nMovies = 100; % number of nucleus movies


% arrangements
longSliceArray = cat(3,movieData.sliceArray);

% initialize arrays
mean_bcd_vec = NaN(size(longSliceArray,3),1);
ctrl_bcd_vec = NaN(size(longSliceArray,3),1);
sym_cell = cell(1,size(longSliceArray,3));
mean_ent_vec = NaN(size(longSliceArray,3),1);
ref_ent_vec = NaN(size(longSliceArray,3),1);
ctrl_ent_vec = NaN(size(longSliceArray,3),1);

gini_vec = NaN(size(longSliceArray,3),1);

rand_gini_vec = NaN(size(longSliceArray,3),1);


time_vec = [movieData.time_vec];
fluo_vec = [movieData.fluo_vec];


% calculate fluorescence scale

try
  parpool(18)
catch
end

parfor s = 1:size(longSliceArray,3)
         

    nuc_slice = longSliceArray(:,:,s);
    
    % generate mask
    nc_mask = false(size(nuc_slice));
    nc_mask(ceil(size(nc_mask,1)/2),ceil(size(nc_mask,2)/2)) = true;
    nc_mask = bwdist(nc_mask)<=nucleus_snip_radius_px;
    
    % calculate mean level and entropy
    px_vec = nuc_slice(nc_mask==1);
    mean_bcd_vec(s) = nanmean(px_vec);    
    px_vec_norm = px_vec / sum(px_vec) + 1e-10;
    mean_ent_vec(s) = -sum(px_vec_norm.*log2(px_vec_norm));
    
    px_uni = ones(size(px_vec))/length(px_vec);
    px_rand = rand(size(px_vec))*rand();
    px_rand_norm = px_rand / sum(px_rand);
    ref_ent_vec(s) = -sum(px_uni.*log2(px_uni));
    ctrl_ent_vec(s) = -sum(px_rand_norm.*log2(px_rand_norm));
    ctrl_bcd_vec(s) = mean(px_rand);
    
    sym_cell{s} = px_vec_norm;
    
    % calculate gini coefficient
    diffs = triu(px_vec-px_vec');
    gini_vec(s) = sum(sum(abs(diffs(diffs~=0)))) / (2*length(px_vec)^2*nanmean(px_vec));
    
    diffs_rand = triu(px_rand-px_rand');
    rand_gini_vec(s) = sum(sum(abs(diffs_rand(diffs_rand~=0)))) / (2*length(px_rand)^2*nanmean(px_rand));
end


%%
close all

info_gain_vec = gini_vec;
info_gain_vec_ctrl = ctrl_ent_vec;

figure;
hold on
scatter(time_vec/60,info_gain_vec)
xlabel('minutes into nc14')
ylabel('intra-nuclear Bcd info content (bits)')


figure;
hold on
scatter(fluo_vec,info_gain_vec)
xlabel('spot intensity (au)')
ylabel('intra-nuclear Bcd info content (bits)')

figure;
hold on
scatter(mean_bcd_vec,info_gain_vec)
xlabel('[Bcd]')
ylabel('intra-nuclear Bcd info content (bits)')

figure;
hold on
scatter(ctrl_bcd_vec,info_gain_vec_ctrl)
xlabel('[Bcd] (cotnrol)')
ylabel('rand info content (bits)')

%%
high_info_indices = find(info_gain_vec >= 10.8 );
low_info_indices = find(info_gain_vec <= 10.4 & info_gain_vec >=10.3);
bins = linspace(0,5e-3);
ind = 10;
figure;
hold on
histogram(sym_cell{high_info_indices(ind)},bins)
histogram(sym_cell{low_info_indices(ind)},bins)

%%
close all
figure;
imagesc(longSliceArray(:,:,high_info_indices(ind))/mean_bcd_vec(high_info_indices(ind)))
colorbar
% caxis([0 2e3])

figure;
imagesc(longSliceArray(:,:,low_info_indices(ind))/mean_bcd_vec(low_info_indices(ind)))
colorbar
% caxis([0 2e3])