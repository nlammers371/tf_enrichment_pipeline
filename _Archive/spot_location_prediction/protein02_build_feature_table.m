% Script to investigate spatio-temporal dynamics of protein distributions
% within nuclei. Ultimate goal is to find a way to infer position of active
% loci using protein channel alone
clear
close all

% define core ID variables
project = 'Dl-Ven_snaBAC-mCh';
% project = 'Dl-Ven_hbP2P-mCh';
dropboxFolder =  'E:\Nick\LivemRNA\Dropbox\';
dataPath = [dropboxFolder 'ProcessedEnrichmentData\' project '/'];
% load training set
load([dataPath 'spot_loc_train_set.mat'],'training_struct')
% set spot radius used for pixel labels
spot_rad = 2; % in pixels

% reshape protein array
nc_protein_array = training_struct.nc_protein_array;
nc_protein_array = reshape(nc_protein_array,size(nc_protein_array,1),size(nc_protein_array,2),[]);
nc_mcp_array = training_struct.nc_mcp_array;
nc_mcp_array = reshape(nc_mcp_array,size(nc_protein_array,1),size(nc_protein_array,2),[]);
% nc_protein_array = nc_protein_array(:,:,1:50);
% extract nucleus mask
nc_mask = training_struct.nc_mask;
nc_mask_array = repmat(nc_mask,1,1,size(nc_protein_array,3));
nc_protein_normed = nc_protein_array - nansum(nansum(nc_protein_array.*nc_mask_array))/sum(nc_mask(:));
nc_protein_normed = nc_protein_normed ./ nanstd(nanstd(nc_protein_normed));
disp('generating region labels...')
% generate array containing region labels
region_label_array = false(size(nc_protein_array));
[x_ref, y_ref] = meshgrid(1:size(region_label_array,2),1:size(region_label_array,1));
for i = 1:size(region_label_array,3)
    xp = training_struct.spot_x_vec(i);
    yp = training_struct.spot_y_vec(i);
    rp = sqrt((x_ref-xp).^2+(y_ref-yp).^2);
    % record
    slice = region_label_array(:,:,i);
    slice(rp<=spot_rad) = true;
    region_label_array(:,:,i) = slice;
end
id_mat = reshape(1:size(region_label_array,3),1,1,[]).*ones(size(nc_protein_array));

disp('generating gaussian blur features...')
tic
% Gaussian blur
fg1 = imgaussfilt(nc_protein_normed,1);
fg3 = imgaussfilt(nc_protein_normed,3);
fg5 = imgaussfilt(nc_protein_normed,5);
fg7 = imgaussfilt(nc_protein_normed,7);
toc
disp('generating gaussian DoG features...')
tic
% DoG
DoG31 = fg1 - fg3;
DoG51 = fg1 - fg5;
DoG53 = fg3 - fg5;
DoG75 = fg5 - fg7;
toc
disp('generating gradient features...')
tic
% Generate gradient images of blurs
d_raw = NaN(size(fg1));
d_g1 = NaN(size(fg1));
d_g3 = NaN(size(fg1));
d_g5 = NaN(size(fg1));
d_g7 = NaN(size(fg1));

for i = 1:size(fg1,3)
    d_raw(:,:,i) = imgradient(nc_protein_normed(:,:,i));
    d_g1(:,:,i) = imgradient(fg1(:,:,i));
    d_g3(:,:,i) = imgradient(fg3(:,:,i));
    d_g5(:,:,i) = imgradient(fg5(:,:,i));
    d_g7(:,:,i) = imgradient(fg7(:,:,i));
end
toc
disp('compileing feature array...')
nc_mask_ft = nc_mask_array(:)' & ~isnan(nc_protein_normed(:))';
% generate feature array
feature_array = [region_label_array(:) nc_protein_normed(:) fg1(:) fg3(:) fg5(:) fg7(:)...
    DoG31(:) DoG51(:) DoG53(:) DoG75(:) d_raw(:) d_g1(:) d_g3(:) d_g5(:) d_g7(:)];
feature_array = feature_array(nc_mask_ft==1,:);


% array of other relevant features
id_array = cat(2,repelem((1:size(nc_mask_array,3)),numel(nc_mask))',...
    repmat((1:numel(nc_mask))',size(nc_mask_array,3),1),...
    repelem(training_struct.spot_particle_vec,numel(nc_mask))',...
    repelem(training_struct.spot_frame_vec,numel(nc_mask))', ...
    repelem(training_struct.spot_time_vec,numel(nc_mask))', ...
    nc_mcp_array(:), ...
... %     repelem(training_struct.mf_protein_vec,numel(nc_mask))' ...
    nc_protein_array(:));
id_array = id_array(nc_mask_ft==1,:);

%%%
feature_table = array2table(feature_array,'VariableNames',{'class','raw_normed','g1','g3','g5','g7'...
    ,'DoG31','DoG51','DoG53','DoG75','d_raw','dg1','dg3','dg5','dg7'});
id_table = array2table(id_array,'VariableNames',{'image_num','pixel_num','ptID',...
    'frame','time','raw_mcp','raw_protein'});
%%% generate training, testing, and full image tables
wt_vec = feature_table.class * .99 + .01;
id_index = 1:size(feature_table,1);
train_ids = randsample(id_index,100000,true,wt_vec);
test_ids = randsample(id_index(~ismember(id_index,train_ids)),100000,true,...
    wt_vec(~ismember(id_index,train_ids)));
image_ids = randsample(1:size(nc_protein_array,3),100,false);
% training data
feature_table_train = feature_table(train_ids,:);
id_table_train = id_table(train_ids,:);
writetable(feature_table_train,[dataPath 'spot_loc_feature_train.csv'])
writetable(id_table_train,[dataPath 'spot_loc_id_train.csv'])
% testing data
feature_table_test = feature_table(test_ids,:);
id_table_test = id_table(test_ids,:);
writetable(feature_table_test,[dataPath 'spot_loc_feature_test.csv'])
writetable(id_table_test,[dataPath 'spot_loc_id_test.csv'])
% image data
im_ft = ismember(id_table.image_num,image_ids);
feature_table_im = feature_table(im_ft,:);
id_table_im = id_table(im_ft,:);
writetable(feature_table_im,[dataPath 'spot_loc_feature_im.csv'])
writetable(id_table_im,[dataPath 'spot_loc_id_im.csv'])


