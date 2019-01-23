% main04_check_control_selection(project, RawPath)
%
% DESCRIPTION
% Generates figures overlaying segmentation results onto histone channel
% and interface in which user can accept or reject individual frames for
% use in subsequent analyses
%
% ARGUMENTS
% project: master ID variable 
%
% RawPath: Full or relative path to PreProcessed folder (pr
%               equivalent)

%
% OUTPUT: ref_frame_struct: compiled data set

function main04_check_control_selection(project, RawPath)
% Script to check efficacy of control selection procedure

%%% Load trace data
DataPath = ['../../dat/' project '/'];
nucleus_struct_pt = load([DataPath  'nucleus_struct_pt.mat'],'nucleus_struct_pt');
nucleus_struct_pt = nucleus_struct_pt.nucleus_struct_pt;

snip_struct = load([DataPath  'snip_struct.mat'],'snip_struct');
snip_struct = snip_struct.snip_struct;

ref_frame_struct = load([DataPath  'ref_frame_struct.mat'],'ref_frame_struct');
ref_frame_struct = ref_frame_struct.ref_frame_struct;

% Get path to ref images
RefPath = [DataPath '/mf_images/'];
% Define path to write control figures to
WritePath = [DataPath '/ctrl_check_images/'];
mkdir(WritePath)
% get list of control and spot points
ctrl_vec = [nucleus_struct_pt.ctrl_flags];
frame_vec = [nucleus_struct_pt.frames];
xSpot_vec = [nucleus_struct_pt.xPosParticle];
ySpot_vec = [nucleus_struct_pt.yPosParticle];
xNull_vec = [nucleus_struct_pt.xPosNull_edge];
yNull_vec = [nucleus_struct_pt.yPosNull_edge];
set_vec = [];
for i = 1:numel(nucleus_struct_pt)
    set_vec = [set_vec repelem(nucleus_struct_pt(i).setID,numel(nucleus_struct_pt(i).ctrl_flags))];
end

set_vec = set_vec(ctrl_vec==1);
frame_vec = frame_vec(ctrl_vec==1);
xSpot_vec = xSpot_vec(ctrl_vec==1);
ySpot_vec = ySpot_vec(ctrl_vec==1);
xNull_vec = xNull_vec(ctrl_vec==1);
yNull_vec = yNull_vec(ctrl_vec==1);
set_index = unique(set_vec);
set_frame_array = unique([set_vec' frame_vec'],'row');
% check to see if there are already files in the write directory
write_files = dir([WritePath '*.tif']);
overwrite = 1;
fnames = {};
if numel(write_files) > 0
    ip = input('Files detected in output folder. Do you want to use existing files? (1=Yes/0=No)');
    if ip
        overwrite = 0;
        fnames = {write_files.name};
    end
end
   
% generate source index 
src_cell = {1,numel(set_index)};
pt_id_vec = [nucleus_struct_pt.ParticleID];
for s = 1:numel(set_index)
    ind = find([nucleus_struct_pt.setID]==set_index(s)&~isnan(pt_id_vec),1);   
    src = nucleus_struct_pt(ind).source_path;    
    src_cell{s} = src;     
end
cm = jet(128);
parfor i = 1:size(set_frame_array,1)
    setID = set_frame_array(i,1);
    frame = set_frame_array(i,2);
    % generate write name
    out_name = ['ctrl_check_frame_' sprintf('%03d',frame) '_set' sprintf('%03d',setID) '.tif'];        
    % if using preexisting files, check to see if this set-frame combo
    % exists already
    if ~overwrite
        if sum(contains(fnames,out_name)) > 0
            continue
        end
    end
    % get nucleus and spot position info
    xs = xSpot_vec(set_vec==setID&frame_vec==frame);
    ys = ySpot_vec(set_vec==setID&frame_vec==frame);    
    xn = xNull_vec(set_vec==setID&frame_vec==frame);
    yn = yNull_vec(set_vec==setID&frame_vec==frame);
    cm_vec = mod(5*(1:numel(yn)),128)+1;    
        
    % generate ref slice name and load
    ref_name = [RefPath '/ref_frame_' sprintf('%03d',frame) '_set' sprintf('%03d',setID) '.mat'];
    rf = load(ref_name,'rf');
    rf = ~isnan(rf.rf);
    % generate histone slice name and load
    src = src_cell{set_index==setID};    
    his_slice = imread([RawPath src '/' src '-His_' sprintf('%03d',frame) '.tif']);
    ct_fig = figure('Visible','off');
    imshow(imadjust(his_slice));
    hold on
    % generate composite  image
    boundaries = bwboundaries(rf);
    for k=1:numel(boundaries)
        b = boundaries{k};
        plot(b(:,2),b(:,1),'g','LineWidth',1);
    end
    scatter(xs,ys,30,cm(cm_vec,:),'filled','MarkerEdgeColor','black')
    hold on
    scatter(xn,yn,30,cm(cm_vec,:),'filled','MarkerEdgeColor','black','Marker','^')    
    saveas(ct_fig, [WritePath '/' out_name])
end

%% Final QC control for frames. Accepted frames will be used for enrichment calculations
check = input('Perform manual review of nucleus segmentation and control assignment? (1=Yes/0=No)');
% Check for pre-existing qc array
set_ref = set_frame_array(:,1);
frame_ref = set_frame_array(:,2);
if check
%     if isfile([DataPath 'qc_ref_array.mat'])
%         ip = input('Preexisting QC file found. Do you wish to use this file? (y/n)');
%         if strcmp(ip,'n')
%             disp('Overwriting Old File')
%             ref_qc_final = false(1,size(set_frame_array,1));        
%         else
%             disp('Using Old File')
%             load([DataPath 'qc_ref_array.mat'])        
%             ref_qc_final = ref_array(:,end);
%         end
%     else
%         ref_qc_final = false(1,size(set_frame_array,1));        
%     end
    ref_qc_final = false(1,size(set_frame_array,1));        
    i = 1;
    iter = 1;
    while i <= numel(ref_qc_final)
        a_str = 'N';
        if  ref_qc_final(i)
            a_str = 'Y';
        end
        setID = set_frame_array(i,1);
        frame = set_frame_array(i,2);        
        
        next_set_ind = find(set_frame_array(:,1)>setID,1);
        if isempty(next_set_ind)
            next_set_ind = numel(ref_qc_final)+1;
        end
        ref_name =  [WritePath '/ctrl_check_frame_' sprintf('%03d',frame) '_set' sprintf('%03d',setID) '.tif'];    
        im = imread(ref_name);
        qc_fig = figure;
        imshow(im)
        title(['Set ' num2str(setID) ' Frame ' num2str(frame) ' (Current Status: ' a_str ')'])
        ct = waitforbuttonpress;
        cc = get(qc_fig,'currentcharacter');
        disp('Press 1 to approve, 0 to disapprove, 2 to approve current frame and all following')        
    %     value = double(get(gcf,'CurrentCharacter'));     
        if strcmp(cc,'x')
            disp('saving and exiting')
            i = numel(ref_qc_final)+1;
        elseif strcmp(cc,'2')            
            ref_qc_final(i:next_set_ind-1) = 1;
            i = next_set_ind;
        else
            ref_qc_final(i) = str2num(cc);   
            i = i + 1;
        end
        iter = iter + 1;
        close all    
    end
else
    warning('Skipping manual review. All frames will be included in analysis set')
    ref_qc_final = true(1,size(set_frame_array,1));        
end
nucleus_struct_ctrl = nucleus_struct_pt;
snip_struct_ctrl = snip_struct;
% aply frame selection to nucleus struct
for i = 1:numel(nucleus_struct_ctrl)
    f_vec = nucleus_struct_ctrl(i).frames;
    s_vec = repelem(nucleus_struct_ctrl(i).setID,numel(f_vec));
    ctrl_edge = nucleus_struct_ctrl(i).ctrl_flags;    
    for j = 1:numel(f_vec)
        if ~isnan(ctrl_vec)
            if ctrl_edge(j)
                ctrl_edge(j) = ref_qc_final(f_vec(j)==frame_ref&s_vec(j)==set_ref);        
            end
        end
    end
    nucleus_struct_ctrl(i).ctrl_flags_final = ctrl_edge;    
end
% aply frame selection to snip struct
snip_frames = snip_struct_ctrl.ptFrame_vec;
snip_sets = snip_struct_ctrl.set_vec;
ctrl_flags_final = snip_struct.ctrl_flags_edge;
for i = 1:numel(snip_sets)    
    if ~isnan(ctrl_flags_final(i))
        if ctrl_flags_final(i)
            ctrl_flags_final(i) = ref_qc_final(snip_frames(i)==frame_ref&snip_sets(i)==set_ref);    
        end
    end
end
snip_struct_ctrl.ctrl_flags_final = ctrl_flags_final;
save([DataPath 'nucleus_struct_ctrl.mat'],'nucleus_struct_ctrl')
save([DataPath 'snip_struct_ctrl.mat'],'snip_struct_ctrl')
