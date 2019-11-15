clear
close all
oldPath = "E:\Nick\LivemRNA\Dropbox (Personal)\ProcessedEnrichmentData\Dl-Ven_snaBAC-mCh\nucleus_struct.mat";
load(oldPath);
nucleus_struct_orig = nucleus_struct;
newPath = "E:\Nick\LivemRNA\Dropbox (Personal)\ProcessedEnrichmentData\Dl-Ven_New_snaBAC-mCh\nucleus_struct.mat";
load(newPath);
nucleus_struct_new = nucleus_struct;

%%%
time_index = 0:60:50*60;
new_set_index = unique(floor([nucleus_struct_new.ParticleID]));
new_set_index = new_set_index(~isnan(new_set_index));
orig_set_index = unique(floor([nucleus_struct_orig.ParticleID]));
orig_set_index = orig_set_index(~isnan(orig_set_index));

% generate reference arrays
new_time_ref = round([nucleus_struct_new.time]/60)*60;
orig_time_ref = round([nucleus_struct_orig.time]/60)*60;
new_offset_ref = [nucleus_struct_new.fluoOffset];
orig_offset_ref = [nucleus_struct_orig.fluoOffset];

orig_set_ref = [];
for i = 1:numel(nucleus_struct_orig)
    orig_set_ref = [orig_set_ref repelem(nucleus_struct_orig(i).setID,numel(nucleus_struct_orig(i).time))];
end

new_set_ref = [];
for i = 1:numel(nucleus_struct_new)
    new_set_ref = [new_set_ref repelem(nucleus_struct_new(i).setID,numel(nucleus_struct_new(i).time))];
end
%
% calculate average offset for new and old
new_offset_array = NaN(numel(new_set_index),numel(time_index));
orig_offset_array = NaN(numel(orig_set_index),numel(time_index));
% 
for t = 1:numel(time_index)    
    for i = 1:numel(orig_set_index)
        orig_offset_array(i,t) = nanmean(orig_offset_ref(orig_time_ref==time_index(t) & orig_set_ref==orig_set_index(i)));
    end
    for i = 1:numel(new_set_index)
        new_offset_array(i,t) = nanmean(new_offset_ref(new_time_ref==time_index(t) & new_set_ref==new_set_index(i)));
    end
end

%%
close all
figure;
hold on


p1 = plot(nanmean(new_offset_array,1)','LineWidth',2);
p2 = plot(nanmean(orig_offset_array,1)','LineWidth',2);
plot(new_offset_array','Color',[p1.Color .5])
plot(orig_offset_array','Color',[p2.Color .5])


legend([p2 p1],'original line','new line')        
xlabel('minutes into nc14')
ylabel('offset (au)')