clear
close all
clc

addpath(genpath('./lib'))

%% Initialization

projectNameCell = {'optokni_eve4+6_ON_CONST'}; 


%%
% color to be used
k_green = brighten([38 142 75]/256,.4);
color_green = [38 143 75]/256; % color from Jake
mRNA_red = brighten([212 100 39]/256,.2);

eYFP_background = 375698.13; %prctile(double(knirps_vec_long),1);

ap_lim = 0.02; % AP range for analysis, 0.02 seems to be a reasonable number

time_threshold = 2; %2min seems to be reasonable

% histogram parameters
binNum = 15;% best 14; ok, 15
binMax = 8;% best 8
%binNum = 10;%14
%binMax = 8;
edges = linspace(0,binMax,binNum);
%edges1 = linspace(0,10,11); % bin for histogram
%edges1 = linspace(0,15,16);
%edges = linspace(0,6,9); % bin for histogram

% temporary correction
%correction_factor = 1.346; % calibration factor
%correction_factor = 1.243;
%correction_factor = 1;

% timerange to analyze for response time
analysis_range = 8;

condition_filter = [];
knirps_int_last_on_all_cond = [];
response_time_all_cond = [];
fluo_mean_all_cond = [];
knirps_mean_all_cond = [];
%data_filter_all_cond = [];

%% Start analysis

for project = 1:length(projectNameCell)

projectName = projectNameCell{project};
liveProject = LiveEnrichmentProject(projectName);
resultsRoot = [liveProject.dataPath filesep];

% load data
load([resultsRoot 'spot_struct.mat'])

embryo = [];

%% For "optokni_eve4+6_ON_CONST"
% Embryo 21
embryo(1).expID = 1;
embryo(1).frame_on = 62;

% Embryo 22
embryo(2).expID = 2;
embryo(2).frame_on = 50;

% Embryo 24
embryo(3).expID = 3;
embryo(3).frame_on = 42;


%% Figure: plot mean fluorescence vs time (not aligned)

data_filter_full = [];
silence_time_full = [];
response_time_full = [];
silence_time_full_frac = [];
response_time_full_frac = [];
knirps_last_on_full = [];
knirps_int_last_on_full = [];
fluo_mean_full = [];
knirps_mean_full = [];

time_aligned_full_long = [];
fluo_full_long = [];
knirps_full_long = [];
frame_full_long = [];

for i = 1:length(embryo)

    expID = embryo(i).expID;
    frame_on = embryo(i).frame_on;
    
    nBoots = 100;

    ever_on_vec = [];
    mean_ap = [];
    time_orig_long = [];
    fluo_orig_long = [];
    frame_orig_long = [];
    off_time_long = [];
    knirps_orig_long = [];
    last_on_long = [];
    last_on_long_frac = [];
    first_on_long = [];
    first_on_long_frac = [];
    knirps_last_on_long = [];
    knirps_last_on_long_frac = [];
    knirps_int_last_on_long = [];
    fluo_mean_long = [];
    knirps_mean_long = [];

    count = 0;
    count_zeros = 0;

    for j = 1:length(spot_struct)

        if (spot_struct(j).TraceQCFlag == 1) && (spot_struct(j).setID == expID)
            % extract core vectors 

            % extract core vectors 
            fluo_vec_orig = spot_struct(j).fluo;
            fluo_vec_orig_zeros = fluo_vec_orig;
            fluo_vec_orig_zeros(isnan(fluo_vec_orig_zeros)) = 0;
            
            time_vec_orig = spot_struct(j).time;
            frame_vec_orig = spot_struct(j).frames;
            knirps_vec_orig = spot_struct(j).rawNCProtein;
            ap_vec_orig = spot_struct(j).APPosNucleus;

            % calculate mean
            ever_on_orig = any(~isnan(fluo_vec_orig));
            mean_ap_orig = nanmean(ap_vec_orig);
            mean_knirps_orig = nanmean(knirps_vec_orig);

            if ever_on_orig
                last_on_index = find(~isnan(fluo_vec_orig) & (frame_vec_orig <= frame_on),1,'last');
                last_on_frame = frame_vec_orig(last_on_index);
                first_on_index = find(~isnan(fluo_vec_orig) & (frame_vec_orig > frame_on),1);
                first_on_frame = frame_vec_orig(first_on_index);
                frame_on_index = find((frame_vec_orig == frame_on));
                if ~isnan(frame_on_index)
                    frame_delay_index = find((time_vec_orig-time_vec_orig(frame_on_index))>=5*60,1);
                    frame_prior_index = find((time_vec_orig-time_vec_orig(frame_on_index))<=-2*60,1,'last');
                else
                    frame_delay_index = NaN;
                end
            end
            

            if (mean_ap_orig > -ap_lim) && (mean_ap_orig < ap_lim)
               count_zeros = count_zeros + 1;
               time_orig_long = [time_orig_long time_vec_orig];
               frame_orig_long = [frame_orig_long frame_vec_orig];
               fluo_orig_long = [fluo_orig_long fluo_vec_orig];
               knirps_orig_long = [knirps_orig_long knirps_vec_orig];
               

               if ~isempty(last_on_frame) && ~isempty(first_on_frame)
                   last_on_long = [last_on_long last_on_frame];
                   first_on_long = [first_on_long first_on_frame];
                   knirps_last_on_long = [knirps_last_on_long knirps_orig_long(last_on_index)];
                   if last_on_frame<frame_on
                       knirps_int_last_on_long = [knirps_int_last_on_long sum(knirps_orig_long(last_on_index:frame_on_index-1)-eYFP_background)];
                   else
                       knirps_int_last_on_long = [knirps_int_last_on_long NaN];
                   end
                   
                   %if ~isnan(frame_delay_index) && (first_on_index < frame_delay_index)
                   %    fluo_mean_long = [fluo_mean_long mean(fluo_vec_orig_zeros(first_on_index:frame_delay_index))];
                   %else
                   %    fluo_mean_long = [fluo_mean_long NaN];
                   %end
                   
                   if ~isnan(frame_delay_index)
                       fluo_mean_long = [fluo_mean_long mean(fluo_vec_orig_zeros(frame_on_index:frame_delay_index))];
                   else
                       fluo_mean_long = [fluo_mean_long NaN];
                   end
                   
                   if ~isnan(frame_prior_index)
                       knirps_mean_long = [knirps_mean_long mean(knirps_vec_orig(frame_prior_index:frame_on_index-1))];
                   else
                       knirps_mean_long = [knirps_mean_long NaN];
                   end
                   

               end
               
               if ~isempty(last_on_frame)
                   knirps_last_on_long_frac = [knirps_last_on_long knirps_orig_long(last_on_index)];
                   last_on_long_frac = [last_on_long_frac last_on_frame];
                   if ~isempty(first_on_frame)
                       first_on_long_frac = [first_on_long_frac first_on_frame];
                   else
                       first_on_long_frac = [first_on_long_frac NaN];
                   end
               end
               

               %plot((time_vec_orig)/60,fluo_vec_orig,'Color', [175 175 175]/255);
               count = count + 1;
            end

    %         if (mean_ap(end) > -0.01) && (mean_ap(end) < 0.01)
    %             plot(time_vec-time_vec(end),fluo_vec);
    %             count = count + 1
    %         end

        end

    end
    
    
    % calculate mean knirps and fraction on (before/after perturbation)

    fluo_orig_long(isnan(fluo_orig_long)) = 0;

    frame_len = max(frame_orig_long);
    
    time_vec = zeros(1,frame_len);
    fluo_vec_mean = zeros(frame_len,1);
    fluo_vec_ste = zeros(frame_len,1);
    knirps_vec_mean = zeros(frame_len,1);
    knirps_vec_ste = zeros(frame_len,1);
    frac_on = zeros(frame_len,1);
    frac_on_ste = zeros(frame_len,1);

    fluo_orig_long_zero = fluo_orig_long;
    fluo_orig_long_zero(isnan(fluo_orig_long)) = 0;

    fluo_orig_long_binary = fluo_orig_long;
    fluo_orig_long_binary(isnan(fluo_orig_long)) = 0;
    fluo_orig_long_binary(fluo_orig_long>0) = 1;


    for j = 1:frame_len

        time_filter_long = frame_orig_long==j;

        %if sum(time_filter_long) > 10        
        %    boot_samples_fluo = bootstrp(nBoots,@nanmean,fluo_orig_long(time_filter_long));
        %    fluo_vec_mean(i) = nanmean(boot_samples_fluo);
        %    fluo_vec_ste(i) = std(boot_samples_fluo);

        time_vec(j) = mean(time_orig_long(time_filter_long))/60;

        fluo_vec_mean(j) = nanmean(fluo_orig_long_zero(time_filter_long));
        fluo_vec_ste(j) = std(fluo_orig_long_zero(time_filter_long),'omitnan');

        knirps_vec_mean(j) = nanmean(knirps_orig_long(time_filter_long));
        knirps_vec_ste(j) = std(knirps_orig_long(time_filter_long),'omitnan');

        frac_on(j) = nanmean(fluo_orig_long_binary(time_filter_long));
        frac_on_ste(j) = std(fluo_orig_long_binary(time_filter_long));

    end

    time_vec_on = time_vec-time_vec(frame_on);

    %knirps_vec_mean(time_vec_on>=0) = knirps_vec_mean(time_vec_on>=0)/correction_factor;
    knirps_vec_mean(time_vec_on>=0) = convert_from_458(knirps_vec_mean(time_vec_on>=0));    
    knirps_vec_mean = knirps_vec_mean-eYFP_background;
    
    knirps_last_on_long = knirps_last_on_long - eYFP_background;
    knirps_last_on_long_frac = knirps_last_on_long_frac - eYFP_background;
    knirps_mean_long = knirps_mean_long - eYFP_background;

    % record the result for this embryo
    data_filter = (time_vec(last_on_long) <= time_vec(frame_on)-time_threshold);
    data_filter_full = [data_filter_full double(data_filter)];
    
    silence_time_full = [silence_time_full time_vec(frame_on)-time_vec(last_on_long)];
    response_time_full = [response_time_full time_vec(first_on_long)-time_vec(frame_on)];
    knirps_last_on_full = [knirps_last_on_full knirps_last_on_long];
    knirps_int_last_on_full = [knirps_int_last_on_full knirps_int_last_on_long];
    fluo_mean_full = [fluo_mean_full fluo_mean_long];
    knirps_mean_full = [knirps_mean_full knirps_mean_long];
    
    silence_time_full_frac = [silence_time_full_frac time_vec(frame_on)-time_vec(last_on_long_frac)];
    for i = 1:length(first_on_long_frac)
        if (first_on_long_frac(i))>0
            response_time_full_frac = [response_time_full_frac time_vec(first_on_long_frac(i))-time_vec(frame_on)];
        else
            response_time_full_frac = [response_time_full_frac NaN];
        end
    end
    %response_time_full_frac = [response_time_full_frac time_vec(first_on_long_frac)-time_vec(frame_on)];
    
    
    % record results for combining the traces
    time_vec_aligned = time_orig_long/60 - time_vec(frame_on);
    time_aligned_full_long = [time_aligned_full_long time_vec_aligned];
    fluo_full_long = [fluo_full_long fluo_orig_long];
    knirps_full_long = [knirps_full_long knirps_orig_long];
    frame_full_long = [frame_full_long frame_orig_long];
    
    
    temp_traj_fig  = figure('Position',[10 10 800 800]);
    tiledlayout(2,1)
    nexttile
    hold on
    %time_interp = min(time_vec_on):0.1:max(time_vec_on);
    time_interp = -10:0.1:10;
    frac_on_mean = movmean(frac_on,5);
    frac_on_interp = interp1(time_vec_on(~isnan(frac_on_mean)),frac_on_mean(~isnan(frac_on_mean)),time_interp,'spline');
    frac_on_interp = movmean(frac_on_interp,5);
    %frac_on_interp = interp1(time_vec_on,frac_on,time_interp,'v5cubic');

    errorbar(time_vec_on,knirps_vec_ste,'Color','k','CapSize',0);
    plot(time_vec_on,knirps_vec_mean,'-k','LineWidth',1)
    scatter(time_vec_on,knirps_vec_mean,50,'MarkerFaceColor',k_green,'MarkerEdgeColor','k')
    xlim([-10 7])
    ylim([3.75E5 9E5])
    xlabel(['time relative to perturbation (min)'])
    ylabel(['Knirps concentration (AU)'])
    pbaspect([3 2 1])
    
    nexttile
    hold on
    plot(time_vec_on,frac_on,'.')
    plot(time_interp,frac_on_interp,'-','LineWidth',2);
    scatter(time_vec_on,frac_on,50,'MarkerFaceColor',mRNA_red,'MarkerEdgeColor','k')
    xlim([-10 7])
    ylim([0 1])
    xlabel(['time relative to perturbation (min)'])
    ylabel(['fraction of nuclei on'])
    pbaspect([3 2 1])
    
    %saveas(temp_traj_fig,[FigurePath 'figure_temporal_trajectory_' num2str(i) '.pdf'])
    
    
    mean_fluo_traj_fig  = figure('Position',[10 10 800 800]);
    tiledlayout(2,1)
    nexttile
    hold on
    %time_interp = min(time_vec_on):0.1:max(time_vec_on);
    time_interp = -10:0.1:10;
    frac_on_mean = movmean(frac_on,5);
    frac_on_interp = interp1(time_vec_on(~isnan(frac_on_mean)),frac_on_mean(~isnan(frac_on_mean)),time_interp,'spline');
    frac_on_interp = movmean(frac_on_interp,5);
    %frac_on_interp = interp1(time_vec_on,frac_on,time_interp,'v5cubic');

    errorbar(time_vec_on,knirps_vec_ste,'Color','k','CapSize',0);
    plot(time_vec_on,knirps_vec_mean,'-k','LineWidth',1)
    scatter(time_vec_on,knirps_vec_mean,50,'MarkerFaceColor',k_green,'MarkerEdgeColor','k')
    xlim([-10 7])
    ylim([3.75E5 9E5])
    xlabel(['time relative to perturbation (min)'])
    ylabel(['Knirps concentration (AU)'])
    pbaspect([3 2 1])
    
    nexttile
    hold on
    plot(time_vec_on,fluo_vec_mean,'.')
    %plot(time_interp,fluo_vec_mean,'-','LineWidth',2);
    scatter(time_vec_on,fluo_vec_mean,50,'MarkerFaceColor',mRNA_red,'MarkerEdgeColor','k')
    xlim([-10 7])
    %ylim([0 1])
    xlabel(['time relative to perturbation (min)'])
    ylabel(['mean fluorescence (AU)'])
    pbaspect([3 2 1])

end

%% plot combined results

% calculate mean knirps and fraction on (before/after perturbation)

fluo_full_long(isnan(fluo_full_long)) = 0;

%time_bin_full = double((0:max(frame_full_long)))+0.5;
time_bin_full = linspace(-15,15,58);% 58,56, best
time_vec_plot = (time_bin_full(1:end-1) + time_bin_full(2:end))/2;
time_groups_full = discretize(time_aligned_full_long,time_bin_full);
%time_vec_full = zeros(1,length(time_bin_full)-1);

fluo_vec_full_mean = zeros(length(time_bin_full)-1,1);
fluo_vec_full_ste = zeros(length(time_bin_full)-1,1);

knirps_vec_full_mean = zeros(length(time_bin_full)-1,1);
knirps_vec_full_ste = zeros(length(time_bin_full)-1,1);

frac_on_full = zeros(length(time_bin_full)-1,1);
frac_on_full_ste = zeros(length(time_bin_full)-1,1);

fluo_full_long_zero = fluo_full_long;
fluo_full_long_zero(isnan(fluo_full_long)) = 0;

fluo_full_long_binary = fluo_full_long;
fluo_full_long_binary(isnan(fluo_full_long)) = 0;
fluo_full_long_binary(fluo_full_long>0) = 1;


for j = 1:length(time_bin_full)-1

    time_filter_long_full = time_groups_full==j;

    time_vec_full(j) = mean(time_aligned_full_long(time_filter_long_full));

    fluo_vec_full_mean(j) = nanmean(fluo_full_long_zero(time_filter_long_full));
    fluo_vec_full_ste(j) = std(fluo_full_long_zero(time_filter_long_full),'omitnan');

    knirps_vec_full_mean(j) = nanmean(knirps_full_long(time_filter_long_full));
    knirps_vec_full_ste(j) = std(knirps_full_long(time_filter_long_full),'omitnan');

    frac_on_full(j) = nanmean(fluo_full_long_binary(time_filter_long_full));
    frac_on_full_ste(j) = std(fluo_full_long_binary(time_filter_long_full));

end

%knirps_vec_full_mean(time_vec_plot>=0) = knirps_vec_full_mean(time_vec_plot>=0)/correction_factor;
knirps_vec_full_mean(time_vec_plot>=0) = convert_from_458(knirps_vec_full_mean(time_vec_plot>=0));
knirps_vec_full_mean = knirps_vec_full_mean - eYFP_background;

combined_traj_fig  = figure('Position',[10 10 800 800]);
tiledlayout(2,1)
nexttile

%time_full_interp = -10:0.1:10;
%knirps_vec_full_movmean = movmean(knirps_vec_full_mean,3);
%knirps_vec_full_interp = interp1(time_vec_full(~isnan(knirps_vec_full_movmean)),knirps_vec_full_movmean(~isnan(knirps_vec_full_movmean)),time_full_interp,'spline');
%knirps_vec_full_interp = movmean(knirps_vec_full_interp,1);

hold on
%plot(time_full_interp,knirps_vec_full_interp,'-k','LineWidth',2)
plot(time_vec_plot,knirps_vec_full_mean,'-k','LineWidth',1)
scatter(time_vec_plot,knirps_vec_full_mean,50,'MarkerFaceColor',k_green,'MarkerEdgeColor','k')
%plot(time_vec_plot,knirps_vec_full_mean,'o');
xlim([-10 7])
ylim([3.75E5 10E5])
xlabel(['time relative to perturbation (min)'])
ylabel(['Knirps concentration (AU)'])
pbaspect([3 2 1])

nexttile

time_full_interp = -10:0.1:10;
frac_on_full_mean = movmean(frac_on_full,3);
frac_on_full_interp = interp1(time_vec_full(~isnan(frac_on_full_mean)),frac_on_full_mean(~isnan(frac_on_full_mean)),time_full_interp,'spline');
frac_on_full_interp = movmean(frac_on_full_interp,1);

hold on
plot(time_vec_plot,frac_on_full,'o')
plot(time_full_interp,frac_on_full_interp,'-','LineWidth',2);
scatter(time_vec_plot,frac_on_full,50,'MarkerFaceColor',mRNA_red,'MarkerEdgeColor','k')
xlabel(['time relative to perturbation (min)'])
ylabel(['fraction of nuclei on'])
xlim([-10 7])
ylim([0 0.9])
pbaspect([3 2 1])

%%
%data_filter_all_cond = [data_filter_all_cond data_filter_full];

knirps_int_last_on_final = knirps_int_last_on_full(data_filter_full==1);
response_time_final = response_time_full(data_filter_full==1);
fluo_mean_final = fluo_mean_full((data_filter_full==1));
knirps_mean_final = knirps_mean_full((data_filter_full==1));

response_time_all_cond = [response_time_all_cond response_time_final];
knirps_int_last_on_all_cond = [knirps_int_last_on_all_cond knirps_int_last_on_final];
fluo_mean_all_cond = [fluo_mean_all_cond fluo_mean_final];
condition_filter = [condition_filter project*ones(1,length(knirps_int_last_on_final))];
knirps_mean_all_cond = [knirps_mean_all_cond knirps_mean_final];


%% plot single traces

sample_traces = [];
time_vec_interp = -10:0.1:7;

% extract single-trace data
for j = 1:length(spot_struct)
    
    embryo_num = find(spot_struct(j).setID == [embryo(:).expID]);
    
    if ~isnan(embryo_num)
        frame_on = embryo(embryo_num).frame_on;

        temp_trace = zeros(1,151);

        ap_vec = spot_struct(j).APPosNucleus;
        ap_pos = mean(ap_vec);

        if (ap_pos>=-0.02) && (ap_pos<=0.02)
            spot_fluo = spot_struct(j).fluo;
            knirps_fluo = spot_struct(j).rawNCProtein;
            spot_fluo(isnan(spot_fluo)) = 0;
            
            time_vec = spot_struct(j).time/60;
            frame_vec = spot_struct(j).frames;

            frame_start = frame_vec(1);
            frame_final = frame_vec(end);

            if (frame_on>=frame_start) && (frame_on<=frame_final) && ~isnan(frame_start) && ~isnan(frame_final)
                time_on = time_vec(frame_vec == frame_on);
                time_vec_new = time_vec-time_on;

                spot_fluo_interp = interp1(time_vec_new,spot_fluo,time_vec_interp,'linear');

                sample_traces = [sample_traces;spot_fluo_interp]; 

            end
        end
    end
    
end

last_on_time_long = zeros(size(sample_traces,1),1);
first_on_time_long = zeros(size(sample_traces,1),1);

for i = 1:size(sample_traces,1)
    spot_vec_temp = sample_traces(i,:);
    
    last_on_time = time_vec_interp(find((spot_vec_temp>0) & (time_vec_interp<0),1,'last'));
    first_on_time = time_vec_interp(find((spot_vec_temp>0) & (time_vec_interp>0),1));
    
    if ~isempty(last_on_time) && ~isempty(first_on_time) && (last_on_time<-2)
        last_on_time_long(i) = last_on_time;
        first_on_time_long(i) = first_on_time;
    else
        last_on_time_long(i) = NaN;
        first_on_time_long(i) = NaN;
    end
end

[B,I] = sort(first_on_time_long,'descend');

%I = I((~isnan(B)) & (B<=5));
I = I((~isnan(B)));

sample_traces_fig = figure;
imagesc('XData',time_vec_interp,'CData',sample_traces(I,:))
xlim([-10 7])
ylim([1 length(I)])
xlabel('time relative to perturbation (min)')
%colormap(plasma)
caxis([0 4.5E5])
colorbar
pbaspect([3 1 1])

%saveas(sample_traces_fig,[FigurePath 'figure_ON_sample_traces.pdf'])

%{
%% plot single sample trace

trace_num = I(91);%74, 89, 91, 106,107,117???

single_traces_fig = figure;
plot(time_vec_interp,sample_traces(trace_num,:))
xlabel('time (min)')
ylabel('spot fluorescence (au)')
xlim([-10 7])
pbaspect([4 1 1])

%saveas(single_traces_fig,[FigurePath 'figure_ON_single_trace.pdf'])
%}
%% calculate silenced duration vs response time

analysis_range_sil_dur = 20;

response_time_final = response_time_full(data_filter_full==1);
silence_time_final = silence_time_full(data_filter_full==1);


x = silence_time_final;
y = response_time_final;
filter = (x<analysis_range_sil_dur) & (y<analysis_range);

xFit = x(filter);
yFit = y(filter);

%mdl = fitlm(xFit,yFit,'RobustOpts','on');
%mdl = fitlm(xFit,yFit);
mdl = fitlm(x,y)
p = coefTest(mdl)

b = mdl.Coefficients{1,1};
k = mdl.Coefficients{2,1};
%options = fitoptions('poly1');
%options.Upper = [Inf -3];
%mdl = fit(x(filter)',y(filter)','poly1',options);
%k = mdl.p1;
%b = mdl.p2;

regMdl = regARIMA(0,0,0);
regMdl.Distribution = struct('Name','t','DoF',3);
%regMdl is a regARIMA model object. It is a template for estimation.

%Estimate the regression model with ARIMA errors. Plot the regression line.
estRegMdl = estimate(regMdl,yFit','X',xFit');


xRange = linspace(0,20,1000);
yResult = k*xRange + b;

memory_fig = figure;
%scatter(x(filter),y(filter),50,'filled','MarkerFaceColor',[234 194 100]/255,'MarkerEdgeColor',[0 .5 .5],'LineWidth',0.5)
scatter(x(filter),y(filter),50,'filled','MarkerFaceColor',[115 142 193]/255,'MarkerEdgeColor',[0 .5 .5],'LineWidth',0.5)
hold on
plot(xRange,yResult,'-','LineWidth',2)
%plot(x,mdl)
xlabel('silenced duration (min) before illumination');
ylabel('response time (min)');
xlim([time_threshold analysis_range_sil_dur])
ylim([0 analysis_range])
pbaspect([3 2 1])
%saveas(memory_fig,[FigurePath 'figure_memory.pdf'])


%{
% plot the binned result
binNum_comp = 12;
binMax_comp = 15;
edges_comp = linspace(time_threshold,binMax_comp,binNum_comp);

[~,~,loc]=histcounts(xFit,edges_comp);
meany = accumarray(loc(:),yFit(:))./accumarray(loc(:),1);
stdy = accumarray(loc(:),yFit(:),[],@std)./sqrt(accumarray(loc(:),1));
xmid = 0.5*(edges_comp(1:end-1)+edges_comp(2:end));

test_fig = figure;
%scatter(x(filter),y(filter),30,'filled','MarkerFaceColor','#BFBF99','MarkerEdgeColor',[0 .5 .5],'LineWidth',0.5)
scatter(x(filter),y(filter),25,'filled','MarkerFaceColor',[200 200 200]/256,'LineWidth',0.5)

hold on
%errorbar(xmid, meany, stdy,'- .','CapSize',18,'MarkerSize',20,'Color','#D64D4D')
errorbar(xmid, meany, stdy,'- .','CapSize',18,'MarkerSize',20,'Color','#D64D4D','LineWidth',1)

xlabel('silenced duration (min) before illumination');
ylabel('response time (min)');
xlim([time_threshold analysis_range_sil_dur])
ylim([0 analysis_range])
pbaspect([3 2 1])
%}

%% fit gamma function

%a = response_time_final(response_time_final<=6);
a = response_time_final(response_time_final<=8);
%a = response_time_final;

% fit gamma function
[muhat,muci] = mle(a,'distribution','gamma'); % Generic function
%[muhat,muci] = gamfit(a); % Distribution specific function

x = 0:0.1:20;
y1 = gampdf(x,muhat(1),muhat(2))/(binNum/binMax);
%y2 = gamcdf(x,muhat(1),muhat(2));

response_time_fig = figure;
hold on
h = histogram(a,edges,'Normalization','probability');
h.FaceColor = mRNA_red;
%plot(x,y1,'LineWidth',3,'Color',mRNA_red)
mean(a)

%xlim([0 analysis_range])
xlim([0 8])
ylim([0 0.23])
xlabel('response time (min)')
ylabel('probability')
pbaspect([1 1 1])

%saveas(response_time_fig,[FigurePath 'figure_response_time_hist.pdf'])



%% calculate silenced knirps concentration vs response time

%analysis_range_sil_dur = 20;

knirps_last_on_final = knirps_last_on_full(data_filter_full==1);
response_time_final = response_time_full(data_filter_full==1);
%silence_time_final = silence_time_full(data_filter_full==1);


x = knirps_last_on_final;
y = response_time_final;
%filter = (x<analysis_range_sil_dur) & (y<analysis_range);

%xFit = x(filter);
%yFit = y(filter);

xFit = double(x);
yFit = y;

%mdl = fitlm(xFit,yFit,'RobustOpts','on');
%mdl = fitlm(xFit,yFit);
mdl = fitlm(x,y)
p = coefTest(mdl)

b = mdl.Coefficients{1,1};
k = mdl.Coefficients{2,1};
%options = fitoptions('poly1');
%options.Upper = [Inf -3];
%mdl = fit(x(filter)',y(filter)','poly1',options);
%k = mdl.p1;
%b = mdl.p2;

regMdl = regARIMA(0,0,0);
regMdl.Distribution = struct('Name','t','DoF',3);
%regMdl is a regARIMA model object. It is a template for estimation.

%Estimate the regression model with ARIMA errors. Plot the regression line.
estRegMdl = estimate(regMdl,yFit','X',xFit');


%xRange = linspace(0,20,1000);
%yResult = k*xRange + b;

memory_fig = figure;
%scatter(x(filter),y(filter),50,'filled','MarkerFaceColor',[234 194 100]/255,'MarkerEdgeColor',[0 .5 .5],'LineWidth',0.5)
scatter(x(filter),y(filter),50,'filled','MarkerFaceColor',[115 142 193]/255,'MarkerEdgeColor',[0 .5 .5],'LineWidth',0.5)
hold on
%plot(xRange,yResult,'-','LineWidth',2)
%plot(x,mdl)
xlabel('silenced knirps concentration (AU)');
ylabel('response time (min)');
%xlim([time_threshold analysis_range_sil_dur])
%ylim([0 analysis_range])
pbaspect([1 1 1])
%saveas(memory_fig,[FigurePath 'figure_memory.pdf'])


end






