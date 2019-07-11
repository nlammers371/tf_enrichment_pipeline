clear
close all
addpath('utilities')
% define core ID variables
project = 'Dl-Ven_snaBAC-mCh';
% project = 'Dl-Ven_hbP2P-mCh';
dropboxFolder =  'E:\Nick\LivemRNA\Dropbox\';
dataPath = [dropboxFolder 'ProcessedEnrichmentData\' project '/'];

% load training data 
feature_table_train = readtable([dataPath 'spot_loc_feature_train.csv']);
id_table_train = readtable([dataPath 'spot_loc_id_train.csv']);

% train RF classifier
Mdl = TreeBagger(100,feature_table_train,'class', 'OOBPrediction','on');

figure;
oobErrorBaggedEnsemble = oobError(Mdl);
plot(oobErrorBaggedEnsemble)
xlabel 'Number of grown trees';
ylabel 'Out-of-bag classification error';
%%
% test performance on second chunck of training data
[~, score_array] = predict(Mdl,feature_table_train(test_ids,:));

%% make figures examining model performance
class_vec = feature_table_train(test_ids,:).class;
hist_fig = figure;
hm_cm = flipud(brewermap([],'RdYlBu'));
colormap(hm_cm);
histogram(score_array(class_vec==0,2),linspace(0,1),'Normalization','probability')
hold on
histogram(score_array(class_vec==1,2),linspace(0,1),'Normalization','probability')

%% info gain
H_init = log2(sum(~isnan(class_vec)));
p_vec_norm = score_array(:,2) / nansum(score_array(:,2));
H_new = -nansum(p_vec_norm.*log2(p_vec_norm));
Information = H_init - H_new

%%
% load image data 
feature_table_im = readtable([dataPath 'spot_loc_feature_im.csv']);
id_table_im = readtable([dataPath 'spot_loc_id_im.csv']);

%%
im_options = unique(id_table_im.image_num);
ind1 = im_options(10);
ex_ft = id_table_im.image_num == ind1;
pixel_id_vec = id_table_im(ex_ft,:).pixel_num;
pixel_class_vec = feature_table_im(ex_ft,:).class;
[~, score_array_pd1] = predict(Mdl,feature_table_im(ex_ft,:));
dim = 57;
nc1_mat = zeros(dim,dim);
spot_prob_vec = score_array_pd1(:,2) / nansum(score_array_pd1(:,2));
nc1_mat(pixel_id_vec) = spot_prob_vec;
id1_mat = zeros(dim,dim);
id1_mat(pixel_id_vec) = pixel_class_vec;
[y, x] = find(id1_mat);
yc = nanmean(y);
xc = nanmean(x);
%%
pd_fig1 = figure;
hm_cm = flipud(brewermap([],'RdYlBu'));
colormap(hm_cm);
p = imagesc(nc1_mat);
hold on
scatter(xc,yc,'filled','MarkerFaceColor','black','MarkerEdgeAlpha',0)
colorbar
% caxis([0 1])