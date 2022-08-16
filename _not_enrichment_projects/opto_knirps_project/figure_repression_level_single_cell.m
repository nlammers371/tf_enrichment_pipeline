clc
clear
close all

addpath(genpath('./lib'))

%% create color map
%Green color map v1
cmap_green = uint8([[247,252,245];[229,245,224];[199,233,192];[161,217,155];[116,196,118];[65,171,93];[35,139,69];[0,109,44];[0,68,27]]);

old_stepNum = size(cmap_green,1);
new_stepNum = 256;

cmap_green_1 = interp1(linspace(0,1,old_stepNum),double(cmap_green(:,1)),linspace(0,1,new_stepNum));
cmap_green_2 = interp1(linspace(0,1,old_stepNum),double(cmap_green(:,2)),linspace(0,1,new_stepNum));
cmap_green_3 = interp1(linspace(0,1,old_stepNum),double(cmap_green(:,3)),linspace(0,1,new_stepNum));

cmap_green_new = [cmap_green_1' cmap_green_2' cmap_green_3']/256;

%Green color map v2
cmap_green_v2 = [[247,252,253];[229,245,249];[204,236,230];[153,216,201];[102,194,164];[65,174,118];[35,139,69];[0,109,44];[0,68,27]];

old_stepNum = size(cmap_green_v2,1);
new_stepNum = 256;

cmap_green_1_v2 = interp1(linspace(0,1,old_stepNum),double(cmap_green_v2(:,1)),linspace(0,1,new_stepNum));
cmap_green_2_v2 = interp1(linspace(0,1,old_stepNum),double(cmap_green_v2(:,2)),linspace(0,1,new_stepNum));
cmap_green_3_v2 = interp1(linspace(0,1,old_stepNum),double(cmap_green_v2(:,3)),linspace(0,1,new_stepNum));

cmap_green_new_v2 = [cmap_green_1_v2' cmap_green_2_v2' cmap_green_3_v2']/256;

%Red color map v1
cmap_red = [[255,255,229];[255,247,188];[254,227,145];[254,196,79];[254,153,41];[236,112,20];[204,76,2];[153,52,4];[102,37,6]];

old_stepNum = size(cmap_red,1);
new_stepNum = 256;

cmap_red_1 = interp1(linspace(0,1,old_stepNum),double(cmap_red(:,1)),linspace(0,1,new_stepNum));
cmap_red_2 = interp1(linspace(0,1,old_stepNum),double(cmap_red(:,2)),linspace(0,1,new_stepNum));
cmap_red_3 = interp1(linspace(0,1,old_stepNum),double(cmap_red(:,3)),linspace(0,1,new_stepNum));

cmap_red_new = [cmap_red_1' cmap_red_2' cmap_red_3']/256;

%Color for lines
color_green = [38 143 75]/256;
color_red = [209 82 5]/256;

left_color = color_green;
right_color = color_red;


% Define some other colors  
yw = [234 194 100]/255; % yellow
bl = [115 143 193]/255; % blue
gr = [191 213 151]/255; % green

%% initialization
ap_lim = 0.02;
eYFP_background = 375698.13; %prctile(double(knirps_vec_long),1);

projectName = {'optokni_eve4+6_WT','optokni_eve4+6_ON_CONST','optokni_eve4+6_ON_LOW_FULL'}; 

spot_traces_WT = [];
spot_traces_CONST = [];
spot_traces_LOW = [];

knirps_traces_WT = [];
knirps_traces_CONST = [];
knirps_traces_LOW = [];

for i = 1:length(projectName)

liveProject = LiveEnrichmentProject(projectName{i});
resultsRoot = [liveProject.dataPath filesep];    

% load data
load([resultsRoot 'spot_struct.mat'])
%FigurePath = [liveProject.figurePath 'WT_vs_CONST' filesep];
%mkdir(FigurePath)

ever_on_vec = zeros(size(spot_struct));
off_time_vec = NaN(size(spot_struct));
on_time_vec = NaN(size(spot_struct));
off_spot_fluo = NaN(size(spot_struct));
off_knirps_vec = NaN(size(spot_struct));
off_ap = NaN(size(spot_struct));
on_knirps_vec = NaN(size(spot_struct));
on_ap = NaN(size(spot_struct));
mean_ap = NaN(size(spot_struct));
mean_knirps = NaN(size(spot_struct));
mean_fluo = NaN(size(spot_struct));
% initialize longform vectors for regression
ap_vec_long = [];
time_vec_long = [];
knirps_vec_long_raw = [];
fluo_raw_long = [];
fluo_zeros_long = [];
mRNA_vec_long = [];

post_turn_on_flags = [];
post_turn_off_flags = [];
ever_on_flags = [];

for j = 1:length(spot_struct)
  
    % extract core vectors 
    fluo_vec = spot_struct(j).fluo;
    time_vec = spot_struct(j).time;
    knirps_vec = spot_struct(j).rawNCProtein;
    ap_vec = spot_struct(j).APPosNucleus;
    
    %if time_vec(end) - time_vec(1) >=20*60
        
        % make average vectors        
      
        % get off and on indices
        ever_on_vec(j) = any(~isnan(fluo_vec));
        mean_ap(j) = nanmean(ap_vec);
        mean_knirps(j) = nanmean(knirps_vec);
        
        post_on_vec = zeros(size(ap_vec));
        post_off_vec = zeros(size(ap_vec));
        if ever_on_vec(j)
            start_i = find(~isnan(fluo_vec),1);
            stop_i = find(~isnan(fluo_vec),1,'last');
            if true%stop_i < length(fluo_vec)-10
                off_time_vec(j) = time_vec(stop_i);
                off_knirps_vec(j) = knirps_vec(stop_i);
                off_spot_fluo(j) = fluo_vec(stop_i);
                off_ap(j) = ap_vec(stop_i);
                post_off_vec(stop_i+1:end) = 1;
            end
            
            if start_i > 1
                on_time_vec(j) = time_vec(start_i);
                on_knirps_vec(j) = knirps_vec(start_i);
                on_ap(j) = ap_vec(start_i);
                post_on_vec(start_i+1:end) = 1;
            end
        end
        
        % make regression vectors
        ever_on_flags = [ever_on_flags repelem(ever_on_vec(j),length(ap_vec))];
        post_turn_on_flags = [post_turn_on_flags post_on_vec];
        post_turn_off_flags = [post_turn_off_flags post_off_vec];
        fluo_raw_long = [fluo_raw_long fluo_vec];
        
        fluo_zeros = fluo_vec;
        all_zeros = fluo_vec;
        fluo_zeros(post_on_vec&~post_off_vec&isnan(fluo_vec)) = 0;
        fluo_zeros_long = [fluo_zeros_long fluo_zeros];
        
        all_zeros(isnan(all_zeros)) = 0;
        mRNA_vec_long = [mRNA_vec_long all_zeros];
        
        mean_fluo(j) = nanmean(fluo_zeros);
        
        knirps_vec_long_raw = [knirps_vec_long_raw knirps_vec];
        ap_vec_long = [ap_vec_long ap_vec];
        time_vec_long = [time_vec_long time_vec];
        
    %end
end


% Look at long vectors

%set parameters
%nBoots = 1000;
timeBins = 61;

%ap_bins = linspace(-0.2,0.2,26);
ap_bins = linspace(-0.2,0.2,31);
ap_bins_plot = (ap_bins(1:end-1) + ap_bins(2:end))/2;
time_bins = linspace(0,35*60,timeBins);
%time_bins_plot = (time_bins(1:end-1) + time_bins(2:end))/2/60;
time_bins_plot = time_bins(1:end-1)/60;
%knirps_bins = linspace(0,15e5,nLinBins);

% calculate mean vectors
knirps_vec_long = knirps_vec_long_raw - eYFP_background;

ap_groups = discretize(ap_vec_long,ap_bins); 
time_groups = discretize(time_vec_long,time_bins); 
%knirps_groups = discretize(knirps_vec_long,knirps_bins); 

frac_inst_on_time_array_mean = NaN(length(time_bins)-1,length(ap_bins)-1);
frac_on_time_array_mean = NaN(length(time_bins)-1,length(ap_bins)-1);
%frac_on_knirps_array_mean = NaN(length(time_bins)-1,length(ap_bins)-1);
eve_time_array_mean = NaN(length(time_bins)-1,length(ap_bins)-1);
%eve_time_array_full = NaN(length(time_bins)-1,length(ap_bins)-1,nBoots);
knirps_time_array_mean = NaN(length(time_bins)-1,length(ap_bins)-1);

frac_inst_on_time_array_std = NaN(length(time_bins)-1,length(ap_bins)-1);
frac_on_time_array_std = NaN(length(time_bins)-1,length(ap_bins)-1);
%frac_on_knirps_array_std = NaN(length(time_bins)-1,length(ap_bins)-1);
eve_time_array_std = NaN(length(time_bins)-1,length(ap_bins)-1);
knirps_time_array_std = NaN(length(time_bins)-1,length(ap_bins)-1);

frac_on_time_array_num = NaN(length(time_bins)-1,length(ap_bins)-1);
frac_inst_on_time_array_num = NaN(length(time_bins)-1,length(ap_bins)-1);
%frac_on_knirps_array_num = NaN(length(time_bins)-1,length(ap_bins)-1);
eve_time_array_num = NaN(length(time_bins)-1,length(ap_bins)-1);
knirps_time_array_num = NaN(length(time_bins)-1,length(ap_bins)-1);


frac_on_time_array_ste = NaN(length(time_bins)-1,length(ap_bins)-1);
frac_on_knirps_array_ste = NaN(length(time_bins)-1,length(ap_bins)-1);
eve_time_array_ste = NaN(length(time_bins)-1,length(ap_bins)-1);
knirps_time_array_ste = NaN(length(time_bins)-1,length(ap_bins)-1);

for t = 1:length(time_bins)
    t
    for a = 1:length(ap_bins)
        time_window_filter = time_groups==t & ap_groups==a & post_turn_on_flags;
        if sum(time_window_filter) > 10
            still_on_flags_time = ~post_turn_off_flags(time_window_filter);
            inst_on_flags_time = mRNA_vec_long(time_window_filter)>0;
            %boot_samples_time = bootstrp(nBoots,@mean,still_on_flags_time);
            %frac_on_time_array_mean(t,a) = mean(boot_samples_time);
            %frac_on_time_array_ste(t,a) = std(boot_samples_time);
            frac_on_time_array_mean(t,a) = nanmean(still_on_flags_time);
            frac_on_time_array_num(t,a) = sum(~isnan(still_on_flags_time));
            frac_inst_on_time_array_mean(t,a) = mean(inst_on_flags_time);
            frac_inst_on_time_array_num(t,a) = length(inst_on_flags_time);
            
            %frac_on_time_array_std(t,a) = std(boot_samples_time);
            %frac_on_time_array_ste(t,a) = frac_on_time_array_std(t,a)/sqrt(frac_on_time_array_num(t,a));
            
            %boot_samples_mRNA = bootstrp(nBoots,@nanmean,mRNA_vec_long(time_window_filter));
            %eve_time_array_mean(t,a) = mean(boot_samples_mRNA);
            %eve_time_array_ste(t,a) = std(boot_samples_mRNA);
            %eve_time_array_full(t,a,:) = boot_samples_mRNA; % NL: note that this is not quite right. Should really be sampling at the level of loci 
            eve_sample = mRNA_vec_long(time_window_filter);
            eve_sample(isnan(eve_sample)) = 0;
            eve_time_array_mean(t,a) = mean(eve_sample);
            eve_time_array_num(t,a) = sum(~isnan(eve_sample));
            eve_time_array_std(t,a) = std(eve_sample);
            eve_time_array_ste(t,a) = eve_time_array_std(t,a)/sqrt(eve_time_array_num(t,a));
            
            %boot_samples_knirps = bootstrp(nBoots,@nanmean,knirps_vec_long_raw(time_window_filter));
            %knirps_time_array_mean(t,a) = nanmean(boot_samples_knirps); % NL: note that this is not quite right. Should really be sampling at the level of loci 
            %knirps_time_array_ste(t,a) = nanstd(boot_samples_knirps);
            knirps_sample = knirps_vec_long(time_window_filter);
            knirps_sample(isnan(knirps_sample)) = 0;
            knirps_time_array_mean(t,a) = mean(knirps_sample);
            knirps_time_array_num(t,a) = sum(~isnan(knirps_sample));
            knirps_time_array_std(t,a) = std(knirps_sample);
            knirps_time_array_ste(t,a) = knirps_time_array_std(t,a)/sqrt(knirps_time_array_num(t,a));
            
        end
                
%         knirps_window_filter = knirps_groups==t & ap_groups==a & post_turn_on_flags;
%         if sum(knirps_window_filter)>10
%             still_on_flags_knirps = ~post_turn_off_flags(knirps_window_filter);
%             boot_samples_knirps = bootstrp(nBoots,@mean,still_on_flags_knirps);
%             frac_on_knirps_array_mean(t,a) = mean(boot_samples_knirps);
%             frac_on_knirps_array_ste(t,a) = std(boot_samples_knirps);
%         end
    end
end
knirps_time_array_mean(isnan(knirps_time_array_mean)) = 0;
eve_time_array_mean(isnan(eve_time_array_mean)) = 0;
frac_inst_on_time_array_mean(isnan(frac_inst_on_time_array_mean)) = 0;

knirps_time_array_mean = movmean(knirps_time_array_mean,3,1);
eve_time_array_mean = movmean(eve_time_array_mean,3,1);
frac_inst_on_time_array_mean = movmean(frac_inst_on_time_array_mean,3,1);

if i == 1
    WT.knirps_mean = knirps_time_array_mean;
    WT.eve_mean = eve_time_array_mean;
    WT.inst_on = frac_inst_on_time_array_mean;
else
    if i == 2
        CONST.knirps_mean = knirps_time_array_mean;
        CONST.eve_mean = eve_time_array_mean;
        CONST.inst_on = frac_inst_on_time_array_mean;
    else
        LOW.knirps_mean = knirps_time_array_mean;
        LOW.eve_mean = eve_time_array_mean;
        LOW.inst_on = frac_inst_on_time_array_mean;
    end
end

num = 0;

for j = 1:length(spot_struct)
    
    temp_trace = zeros(1,151);
    temp_knirps = NaN(1,151);
    
    ap_vec = spot_struct(j).APPosNucleus;
    ap_pos = mean(ap_vec);
    
    if (ap_pos>=-ap_lim) && (ap_pos<=ap_lim)
        num = num+1;
        spot_fluo = spot_struct(j).fluoInterp;
        knirps_protein = spot_struct(j).rawNCProteinInterp;
        
        frame_start = spot_struct(j).timeInterp(1)/20 + 1;
        frame_final = spot_struct(j).timeInterp(end)/20 + 1;

        if ~isnan(frame_start) && ~isnan(frame_final)
            temp_trace(frame_start:frame_final) = spot_fluo;
            temp_knirps(frame_start:frame_final) = knirps_protein;
            if i == 1
                spot_traces_WT = [spot_traces_WT;temp_trace];
                knirps_traces_WT = [knirps_traces_WT;temp_knirps];
            else
                if i == 2
                    spot_traces_CONST = [spot_traces_CONST;temp_trace];
                    knirps_traces_CONST = [knirps_traces_CONST;temp_knirps];
                else
                    spot_traces_LOW = [spot_traces_LOW;temp_trace];
                    knirps_traces_LOW = [knirps_traces_LOW;temp_knirps];
                end
            end
        end
    end
    
end

end

%% plot sample single traces;

sample_num = 100;

time_vec = (0:150)/3;

random_sample_traces_WT = spot_traces_WT(randperm(size(spot_traces_WT, 1)), :);
random_sample_traces_CONST = spot_traces_CONST(randperm(size(spot_traces_CONST, 1)), :);
random_sample_traces_LOW = spot_traces_LOW(randperm(size(spot_traces_LOW, 1)), :);

sample_traces_frac_WT = mean(spot_traces_WT>0,1);
sample_traces_frac_CONST = mean(spot_traces_CONST>0,1);
sample_traces_frac_LOW = mean(spot_traces_LOW>0,1);

sample_traces_mean_WT = mean(spot_traces_WT,1);
sample_traces_mean_CONST = mean(spot_traces_CONST,1);
sample_traces_mean_LOW = mean(spot_traces_LOW,1);

sample_knirps_mean_WT = nanmean(knirps_traces_WT,1)-eYFP_background;
sample_knirps_mean_CONST = convert_from_458(nanmean(knirps_traces_CONST,1))-eYFP_background;
sample_knirps_mean_LOW = nanmean(knirps_traces_LOW,1)-eYFP_background;

single_traces_WT = figure;
imagesc('XData',time_vec,'CData',random_sample_traces_WT(1:sample_num,:))
%imagesc('XData',time_vec,'CData',random_sample_traces_WT)
%imagesc('XData',time_vec,'CData',spot_traces_WT)
xlabel('time (min)')
xlim([0 30])
ylim([1 sample_num])
colorbar
colormap(plasma)
%colormap(turbo)
caxis([0 5E5])
pbaspect([2 1 1])

single_traces_CONST = figure;
imagesc('XData',time_vec,'CData',random_sample_traces_CONST(1:sample_num,:))
%imagesc('XData',time_vec,'CData',ramdom_sample_traces_CONST)
%imagesc('XData',time_vec,'CData',spot_traces_CONST)
xlabel('time (min)')
xlim([0 30])
ylim([1 sample_num])
colorbar
colormap(plasma)
%colormap(turbo)
caxis([0 5E5])
pbaspect([2 1 1])

single_traces_LOW = figure;
imagesc('XData',time_vec,'CData',random_sample_traces_LOW(1:sample_num,:))
%imagesc('XData',time_vec,'CData',random_sample_traces_LOW)
%imagesc('XData',time_vec,'CData',spot_traces_LOW)
xlabel('time (min)')
xlim([0 30])
ylim([1 sample_num])
colorbar
colormap(plasma)
%colormap(turbo)
caxis([0 5E5])
pbaspect([2 1 1])

single_traces_mean = figure;
plot(time_vec,  sample_traces_mean_WT)
hold on
plot(time_vec,  sample_traces_mean_CONST)
plot(time_vec,  sample_traces_mean_LOW)
xlim([0 30])
%ylim([0 1])
xlabel('time (min)')
ylabel('mean spot fluorescence (au)')
pbaspect([2 1 1])


knirps_traces_mean = figure;
plot(time_vec,  sample_knirps_mean_WT)
hold on
plot(time_vec,  sample_knirps_mean_CONST)
plot(time_vec,  sample_knirps_mean_LOW)
xlim([0 30])
xlabel('time (min)')
ylabel('mean knirps fluorescence (au)')
pbaspect([2 1 1])


single_traces_movmean = figure;
plot(time_vec,  movmean(sample_traces_mean_WT,3))
hold on
plot(time_vec,  movmean(sample_traces_mean_CONST,3))
plot(time_vec,  movmean(sample_traces_mean_LOW,3))
xlim([0 30])
%ylim([0 1])
xlabel('time (min)')
ylabel('mean spot fluorescence (au)')
pbaspect([2 1 1])


single_traces_frac_on = figure;
plot(time_vec, sample_traces_frac_WT)
hold on
plot(time_vec,  sample_traces_frac_CONST)
plot(time_vec,  sample_traces_frac_LOW)
xlim([0 30])
ylim([0 1])
xlabel('time (min)')
ylabel('fraction of active nuclei')
pbaspect([3 1 1])
legend('WT','HIGH','LOW')


single_traces_mean_on = figure;
plot(time_vec, sample_traces_mean_WT./sample_traces_frac_WT)
hold on
plot(time_vec,  sample_traces_mean_CONST./sample_traces_frac_CONST)
plot(time_vec,  sample_traces_mean_LOW./sample_traces_frac_LOW)
xlim([0 30])
ylim([0 3E5])
xlabel('time (min)')
%ylabel('fraction of active nuclei')
pbaspect([3 1 1])
legend('WT','HIGH','LOW')


%% plot calculated mean vector

% plot for wildtype
WT_mean_kni_fig = figure;
imagesc(ap_bins_plot,time_bins_plot,WT.knirps_mean)
colorbar
colormap(cmap_green_new_v2)
caxis([0 10E5])
pbaspect([3 2 1])
ylim([0 35])
xlabel('AP position (% embryo length)')
ylabel('time (min)')


WT_mean_eve_fig = figure;
imagesc(ap_bins_plot,time_bins_plot,WT.eve_mean)
colorbar
colormap(cmap_red_new)
caxis([0 3.15E5])
pbaspect([3 2 1])
ylim([5 30])
xlabel('AP position (% embryo length)')
ylabel('time (min)')

% % plot for const export
% CONST_mean_kni_fig = figure;
% imagesc(ap_bins_plot,time_bins_plot,CONST.knirps_mean)
% colorbar
% colormap(cmap_green_new_v2)
% caxis([0 15E5])
% pbaspect([3 2 1])
% ylim([0 35])
% xlabel('AP position (% embryo length)')
% ylabel('time (min)')

CONST_mean_eve_fig = figure;
imagesc(ap_bins_plot,time_bins_plot,CONST.eve_mean)
colorbar
colormap(cmap_red_new)
caxis([0 3.15E5])
pbaspect([3 2 1])
ylim([5 30])
xlabel('AP position (% embryo length)')
ylabel('time (min)')


LOW_mean_eve_fig = figure;
imagesc(ap_bins_plot,time_bins_plot,LOW.eve_mean)
colorbar
colormap(cmap_red_new)
caxis([0 3.15E5])
pbaspect([3 2 1])
ylim([5 30])
xlabel('AP position (% embryo length)')
ylabel('time (min)')


WT_frac_eve_fig = figure;
imagesc(ap_bins_plot,time_bins_plot,WT.inst_on)
colorbar
colormap(cmap_red_new)
%colormap(jet)
caxis([0 1])
pbaspect([3 2 1])
ylim([5 30])
xlabel('AP position (% embryo length)')
ylabel('time (min)')

CONST_frac_eve_fig = figure;
imagesc(ap_bins_plot,time_bins_plot,CONST.inst_on)
colorbar
colormap(cmap_red_new)
%colormap(jet)
caxis([0 1])
pbaspect([3 2 1])
ylim([5 30])
xlabel('AP position (% embryo length)')
ylabel('time (min)')

LOW_frac_eve_fig = figure;
imagesc(ap_bins_plot,time_bins_plot,LOW.inst_on)
colorbar
colormap(cmap_red_new)
%colormap(jet)
caxis([0 1])
pbaspect([3 2 1])
ylim([5 30])
xlabel('AP position (% embryo length)')
ylabel('time (min)')

%% plot sample trace together

sample_WT = 104;
% none?80, 91,95,99,101? 104! 106!
%sample_CONST = 98; 
sample_CONST = 119;
% 55,57,58,150,98!!! 171 173 101! 105! 108! 112!!! 119!!!
%sample_LOW = 150; 
sample_LOW = 65;
% 9? 16, 22, high; 44?; 47!!!, 139!!!  199 146? 65(with high level)!150! sparce: 155!157!


test_panel_fig = figure;
tiledlayout(3,1)
nexttile
plot(time_vec,spot_traces_CONST(sample_CONST,:),'- .')
xlim([14.5 30])
ylim([-0.25E5 5.75E5 ])
pbaspect([3 1 1])
nexttile
plot(time_vec,spot_traces_LOW(sample_LOW,:),'- .')
xlim([14.5 30])
ylim([-0.25E5 5.75E5 ])
pbaspect([3 1 1])
nexttile
plot(time_vec,spot_traces_WT(sample_WT,:),'- .')
xlim([14.5 30])
ylim([-0.25E5 5.75E5 ])
pbaspect([3 1 1])




test_fig = figure;

hold on
plot(time_vec,spot_traces_CONST(sample_CONST,:))
plot(time_vec,spot_traces_LOW(sample_LOW,:))
plot(time_vec,spot_traces_WT(sample_WT,:))
xlim([14.5 30])
ylim([-0.25E5 5.5E5 ])
pbaspect([3 2 1])





