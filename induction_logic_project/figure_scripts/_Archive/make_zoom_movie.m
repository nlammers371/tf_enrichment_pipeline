clear
close all
file_dir = 'E:\LocalEnrichment\Data\RawDynamicsData\2019-11-26\2xDl_Venus_snaBAC_MCPmCherry_Leica_Zoom7_21uW14uW_01\2xDl_Venus_snaBAC_MCPmCherry_Leica_Zoom7_21uW14uW_01.lif';


project = 'zoom7';
DropboxFolder =  'E:\Nick\LivemRNA\Dropbox (Personal)\';
[~, DataPath, FigRoot] =   header_function(DropboxFolder, project);
FigPath = [FigRoot project '\zoom_movies\'];
mkdir(FigPath)

data = bfopen(file_dir);
%% collect protein and transcription channels into separate stacks

nChannels = 3;
nSlices = 12;
series = 2;

T = size(data{series,1},1)/nChannels/nSlices;
im1 = data{series,1}{1,1};

% initialize stacks
mcp_stack = zeros(size(im1,1),size(im1,2),nSlices,T);
protein_stack = zeros(size(im1,1),size(im1,2),nSlices,T);

for i = 1:T*nSlices
    ind = (i-1)*nChannels+1;
    t = ceil(i/nSlices);
    z = mod(i-1,nSlices)+1;
    protein_stack(:,:,z,t) = data{series,1}{ind,1};
    mcp_stack(:,:,z,t) = data{series,1}{ind+1,1};
end

%% specify slices of interest and generate z projections
zRange = 5:7;
protein_stack_pj = imgaussfilt(squeeze(nanmax(protein_stack(:,:,zRange,1:7),[],3)),1.5);
mcp_stack_pj = imgaussfilt(squeeze(nanmean(mcp_stack(:,:,zRange,1:7),3)),1.5);

% apply bleaching correction to protein channel
mean_venus_vec = NaN(1,size(protein_stack_pj,3));
for i = 1:size(protein_stack_pj,3)
    mean_venus_vec(i) = nanmean(nanmean(protein_stack_pj(:,:,i)));
end

protein_stack_pj_bc = protein_stack_pj;
for i = 1:size(protein_stack_pj,3)
    protein_stack_pj_bc(:,:,i) = protein_stack_pj_bc(:,:,i) / mean_venus_vec(i);
end

% convert to grayscale and overlay
rgb_series = NaN(size(protein_stack_pj,1),size(protein_stack_pj,2),size(protein_stack_pj,3),3);

for i = 1:size(protein_stack_pj,3)
    pt_gs = mat2gray(protein_stack_pj_bc(:,:,i));
    mcp_gs = mat2gray(mcp_stack_pj(:,:,i));
    rgb = cat(3,mcp_gs,pt_gs,zeros(size(mcp_gs)));
    imwrite(rgb,[FigPath 'frame_' sprintf('%02d',i) '.tif'])
end
