% Build up data set that can be used for inference. Conduct a few
% exploratory analyses
clear
close all
%%% id variables
Tres = 20;
project = 'kr_eve2_reporter';
%%% Load trace data
ReadPath = ['../../dat/' project '/'];
eve_table = readtable([ReadPath project '_viterbi_fit_data.csv']);
fluo = eve_table.fluo;
% eve_table = eve_table(fluo>=prctile(fluo,80),:);
eve_table = eve_table(~isnan(eve_table.viterbi_state),:);
set_key = readtable([ReadPath project '_set_key.csv'],'Delimiter', ',');
set_index = unique(eve_table.setID);

%%% Load nucleus data
load(['../../dat/' project '/inference_nuclei_' project '_dT' num2str(Tres) '.mat'])

%%% set parameters for protein analysis
snippet_size = 7;
protein_channel = 1;
%%% Read from raw tifs
RawPath = 'E:\Sean\LivemRNA\LivemRNAFISH\Data\PreProcessedData\';
% keyword = 'BacKr';
keyword = 'Rpb1_GFP';
file_list = dir([RawPath '*' keyword '*']);
name_key = {file_list.name};

%%% initialize variable vectors
protein_spot_snippets = NaN(2*snippet_size+1,2*snippet_size+1,size(eve_table,1),1);
protein_spot_vec = NaN(size(eve_table,1),1);
protein_null_snippets = NaN(2*snippet_size+1,2*snippet_size+1,size(eve_table,1),1);
protein_null_vec = NaN(size(eve_table,1),1);
distance_vec_spot = NaN(size(eve_table,1),1);
distance_vec_null = NaN(size(eve_table,1),1);
z_vec = NaN(size(eve_table,1),1);
% binary_vec = NaN(size(eve_table,1),1);
viterbi_vec = NaN(size(eve_table,1),1);
fluo_vec = NaN(size(eve_table,1),1);
fluo_raw = NaN(size(eve_table,1),1);
fluo_raw_snippets = NaN(2*snippet_size+1,2*snippet_size+1,size(eve_table,1),1);
% fluo_raw = NaN(2*snippet_size+1,2*snippet_size+1,size(eve_table,1)*2);
ncID_vec = NaN(size(eve_table,1),1);
time_vec = NaN(size(eve_table,1),1);
ap_vec = NaN(size(eve_table,1),1);
r_vec = NaN(size(eve_table,1),1);
%%% first build up distribution of distances from nuclear center for making
%%% null set
deltas = [eve_table.x-eve_table.xNC eve_table.y-eve_table.yNC];
samp_vec = 1:size(deltas,1);
set_vec = eve_table.setID;

%%% iterate through projects
parfor i = 1:size(eve_table,1)
    setID = eve_table(i,:).setID;
    set_table = eve_table(eve_table.setID==setID,:);
    src = set_key{set_key.setID==setID,2};
    
    % now iterate through observations in set table    
    frame = eve_table(i,:).frame;
    ncID = eve_table(i,:).nucleus_id;
    xp = eve_table(i,:).x-1;
    yp = eve_table(i,:).y-1;
    zp = eve_table(i,:).z;
    v = eve_table(i,:).viterbi_state;
    xn = eve_table(i,:).xNC;
    yn = eve_table(i,:).yNC;
    fluo = eve_table(i,:).fluo;
    time = eve_table(i,:).time;
    AP = eve_table(i,:).AP;
    % load particle frame        
    mcp_name = [src{:} '_' sprintf('%03d',frame) '_z' sprintf('%02d',zp) '_ch0' num2str(1*(2~=protein_channel) + 1) '.tif'];
    mcp_frame = imread([RawPath src{:} '/' mcp_name]);
    mcp_snippet = mcp_frame(max(1,yp-snippet_size):min(size(mcp_frame,1),yp...            
        +snippet_size),max(1,xp-snippet_size):min(size(mcp_frame,2),xp+snippet_size));    
    % load protein frame 
    pt_name = [src{:} '_' sprintf('%03d',frame) '_z' sprintf('%02d',zp) '_ch0' num2str(protein_channel) '.tif'];
    pt_frame = imread([RawPath src{:} '/' pt_name]);
    pt_snippet = pt_frame(max(1,yp-snippet_size):min(size(mcp_frame,1),yp...            
        +snippet_size),max(1,xp-snippet_size):min(size(mcp_frame,2),xp+snippet_size));    
    % now take sample from opposite side of nucleus
    s_id = randsample(samp_vec(set_vec==setID),1);
    x_opp = xn-sign(xp-xn)*abs(deltas(s_id,1));
    y_opp = yn-sign(yp-yn)*abs(deltas(s_id,2));  
        
    r_samp = sqrt((xp-x_opp)^2+(yp-y_opp)^2);
    
    pt_snippet_null = pt_frame(max(1,y_opp-snippet_size):min(size(mcp_frame,1),y_opp...            
        +snippet_size),max(1,x_opp-snippet_size):min(size(mcp_frame,2),x_opp+snippet_size));
    full_size = snippet_size*2 + 1;    
    % only take samples where neither spot nor control was truncated
    if size(pt_snippet,1) ==  full_size && size(pt_snippet,2) == full_size&& ...
        size(pt_snippet_null,1) == full_size && size(pt_snippet_null,2) == full_size        
        
        protein_null_snippets(:,:,i) = pt_snippet_null;
        protein_spot_snippets(:,:,i) = pt_snippet;
        fluo_raw_snippets(:,:,i) = mcp_snippet;
        % save info    
        fluo_raw(i) = mcp_frame(yp,xp);        
        fluo_vec(i) = fluo;        
        protein_spot_vec(i) = pt_frame(yp,xp);        
        protein_null_vec(i) = pt_frame(y_opp,x_opp);        
        distance_vec_spot(i) = sqrt(deltas(i,1)^2+deltas(i,2)^2);
        distance_vec_null(i) = sqrt((x_opp-xn)^2+(y_opp-yn)^2);
        z_vec(i) = zp;    
        viterbi_vec(i) = v;
        time_vec(i) = time;
        ap_vec(i) = AP;
        ncID_vec(i) = ncID;         
        r_vec(i) = r_samp;
    end
end
%%% save
regression_table = array2table([double(protein_spot_vec) double(protein_null_vec) distance_vec_spot distance_vec_null...
    double(r_vec) double(z_vec) double(viterbi_vec)  double(fluo_vec) double(fluo_raw) ncID_vec  time_vec ap_vec],...
    'VariableNames',{'protein_spot','protein_null','distance_spot','distance_null','r_centers',...
    'z','viterbi','fluo','fluo_raw','ncID','time','AP'});

regression_table = regression_table(~isnan(regression_table.ncID),:);
writetable(regression_table,[ReadPath '/' project '_regression_tbl.csv'])

snippet_struct = struct;
snippet_struct.raw_fluo = fluo_raw_snippets;
snippet_struct.protein_spot = protein_spot_snippets;
snippet_struct.protein_null = protein_null_snippets;
save([ReadPath '/' project '_snippets.mat'],'snippet_struct')

