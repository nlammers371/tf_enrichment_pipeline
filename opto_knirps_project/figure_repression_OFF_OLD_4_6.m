clear
close all
clc

addpath(genpath('./lib'))

%% Initialization

projectName = 'optokni_eve4+6_oFF_OLD'; 

liveProject = LiveEnrichmentProject(projectName);
resultsRoot = [liveProject.dataPath filesep];

% load data
load([resultsRoot 'spot_struct.mat'])
FigurePath = [liveProject.figurePath 'repression_dynamics' filesep];
mkdir(FigurePath)

% Embryo 13
embryo(1).expID = 1;
embryo(1).frame_on = 17;
embryo(1).time_on = 5.41981;

% Embryo 14
embryo(2).expID = 2;
embryo(2).frame_on = 49;
embryo(2).time_on = 16.9258;

% Embryo 16
embryo(3).expID = 3;
embryo(3).frame_on = 17;
embryo(3).time_on = 5.43871;

% color to be used
k_green = brighten([38 142 75]/256,.4);
color_green = [38 143 75]/256; % color from Jake
mRNA_red = brighten([212 100 39]/256,.2);

knirps_offset = 2.5e5;%prctile(double(knirps_vec_long),1);

ap_lim = 0.02; % AP range for analysis, 0.02 seems to be a reasonable number

time_threshold = 2; %min
%time_threshold = 1;

% histogram parameters
binNum = 25;
binMax = 12;
edges = linspace(0,binMax,binNum);
%edges1 = linspace(0,10,11); % bin for histogram
%edges1 = linspace(0,15,16);
%edges = linspace(0,6,9); % bin for histogram

% temporary correction
correction_factor = 1.5;
%correction_factor = 1;

% timerange to analyze for response time
analysis_range = 8;


%% Figure: plot mean fluorescence vs time (not aligned)

data_filter_full = [];
on_time_full = [];
response_time_full = [];

count1 = 0;

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
    %last_on_long = [];
    %first_on_long = [];
    last_off_long = [];
    first_off_long = [];

    count = 0;
    for j = 1:length(spot_struct)

        if (spot_struct(j).TraceQCFlag == 1) && (spot_struct(j).setID == expID)
            % extract core vectors 

            % extract core vectors 
            fluo_vec_orig = spot_struct(j).fluo;
            time_vec_orig = spot_struct(j).time;
            frame_vec_orig = spot_struct(j).frames;
            knirps_vec_orig = spot_struct(j).rawNCProtein;
            ap_vec_orig = spot_struct(j).apPosNucleus;

            % calculate mean
            ever_on_orig = any(~isnan(fluo_vec_orig));
            mean_ap_orig = nanmean(ap_vec_orig);
            mean_knirps_orig = nanmean(knirps_vec_orig);

            if ever_on_orig
                %last_on_frame = frame_vec_orig(find(~isnan(fluo_vec_orig) & (frame_vec_orig <= frame_on),1,'last'));
                %first_on_frame = frame_vec_orig(find(~isnan(fluo_vec_orig) & (frame_vec_orig > frame_on),1));
                last_off_frame = frame_vec_orig(find(isnan(fluo_vec_orig) & (frame_vec_orig <= frame_on),1,'last'));
                first_off_frame = frame_vec_orig(find(~isnan(fluo_vec_orig) & (frame_vec_orig > frame_on),1,'last'))+1;
            end

            if (mean_ap_orig > -ap_lim) && (mean_ap_orig < ap_lim)
               count1 = count1+1;
               time_orig_long = [time_orig_long time_vec_orig];
               frame_orig_long = [frame_orig_long frame_vec_orig];
               fluo_orig_long = [fluo_orig_long fluo_vec_orig];
               knirps_orig_long = [knirps_orig_long knirps_vec_orig-knirps_offset];

               %if ~isempty(last_on_frame) && ~isempty(first_on_frame)
               %    last_on_long = [last_on_long last_on_frame];
               %    first_on_long = [first_on_long first_on_frame];
               %end
               
               if ~isempty(last_off_frame) && ~isempty(first_off_frame)
                   last_off_long = [last_off_long last_off_frame];
                   first_off_long = [first_off_long first_off_frame];
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

    %time_bin = 2:max(frame_orig_long);
    %time_groups = discretize(frame_orig_long,time_bin);
    
    frame_len = max(frame_orig_long);
    
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

    %knirps_vec_mean = [];
    %knirps_vec_ste = [];
    time_vec = [];

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

    knirps_vec_mean(time_vec_on<0) = knirps_vec_mean(time_vec_on<0)/correction_factor;
    
    try
        on_time_full = [on_time_full time_vec(frame_on)-time_vec(last_off_long)];
        response_time_full = [response_time_full time_vec(first_off_long)-time_vec(frame_on)];
        
        % record the result for this embryo
        data_filter = (time_vec(last_off_long) <= time_vec(frame_on)-time_threshold);
        data_filter_full = [data_filter_full data_filter];
    end
        
    temp_traj_fig  = figure('Position',[10 10 800 800]);
    tiledlayout(3,1)
    nexttile
    hold on
    %time_interp = min(time_vec_on):0.1:max(time_vec_on);
    time_interp = -5:0.1:10;
    frac_on_mean = movmean(frac_on,5);
    frac_on_interp = interp1(time_vec_on(~isnan(frac_on_mean)),frac_on_mean(~isnan(frac_on_mean)),time_interp,'spline');
    frac_on_interp = movmean(frac_on_interp,5);
    %frac_on_interp = interp1(time_vec_on,frac_on,time_interp,'v5cubic');

    errorbar(time_vec_on,knirps_vec_ste,'Color','k','CapSize',0);
    plot(time_vec_on,knirps_vec_mean,'-k','LineWidth',1)
    scatter(time_vec_on,knirps_vec_mean,50,'MarkerFaceColor',k_green,'MarkerEdgeColor','k')
    xlim([-2.5 7])
    ylim([2E5 10E5])
    xlabel(['time relative to perturbation (min)'])
    ylabel(['Knirps concentration (AU)'])
    pbaspect([3 2 1])
    
    nexttile
    hold on
    plot(time_vec_on,frac_on,'.')
    plot(time_interp,frac_on_interp,'-','LineWidth',2);
    scatter(time_vec_on,frac_on,50,'MarkerFaceColor',mRNA_red,'MarkerEdgeColor','k')
    xlim([-2.5 7])
    ylim([0 1])
    xlabel(['time relative to perturbation (min)'])
    ylabel(['fraction of nuclei on'])
    pbaspect([3 2 1])
    
    nexttile
    hold on 
    plot(time_vec_on,fluo_vec_mean./frac_on,'.')
    %plot(time_interp,frac_on_interp,'-','LineWidth',2);
    scatter(time_vec_on,fluo_vec_mean./frac_on,50,'MarkerFaceColor',mRNA_red,'MarkerEdgeColor','k')
    xlim([-2.5 7])
    ylim([0 3E5])
    xlabel(['time relative to perturbation (min)'])
    ylabel(['mean transcription rate (au)'])
    pbaspect([3 2 1])
    
    %saveas(temp_traj_fig,[FigurePath 'figure_temporal_trajectory_' num2str(i) '.pdf'])

end

%% plot single traces

sample_traces = [];
time_vec_interp = -10:0.1:10;

% extract single-trace data
for j = 1:length(spot_struct)
    
    embryo_num = find(spot_struct(j).setID == [embryo(:).expID]);
    
    if ~isnan(embryo_num)
        frame_on = embryo(embryo_num).frame_on;

        temp_trace = zeros(1,151);

        ap_vec = spot_struct(j).apPosNucleus;
        ap_pos = mean(ap_vec);

        if (ap_pos>=-ap_lim) && (ap_pos<=ap_lim)
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

last_off_time_long = zeros(size(sample_traces,1),1);
first_off_time_long = zeros(size(sample_traces,1),1);

for i = 1:size(sample_traces,1)
    spot_vec_temp = sample_traces(i,:);
    
    last_off_time = time_vec_interp(find((spot_vec_temp==0) & (time_vec_interp<0),1,'last'));
    first_off_time = time_vec_interp(find((spot_vec_temp>0) & (time_vec_interp>0),1,'last'));
    
    if ~isempty(last_off_time) && ~isempty(first_off_time) && (last_off_time<-time_threshold)
        last_off_time_long(i) = last_off_time;
        first_off_time_long(i) = first_off_time;
    else
        last_off_time_long(i) = NaN;
        first_off_time_long(i) = NaN;
    end
end

[B,I] = sort(first_off_time_long,'descend');

%I = I((~isnan(B)) & (B<=5));
I = I((~isnan(B)));

sample_traces_fig = figure;
imagesc('XData',time_vec_interp,'CData',sample_traces(I,:))
xlim([-10 10])
ylim([1 length(I)])
xlabel('time relative to perturbation (min)')
colormap(plasma)
caxis([0 5E5])
pbaspect([3 1 1])

saveas(sample_traces_fig,[FigurePath 'figure_ON_sample_traces.pdf'])

%% plot single sample trace

trace_num = I(10);

single_traces_fig = figure;
plot(time_vec_interp,sample_traces(trace_num,:))
xlabel('time (min)')
ylabel('spot fluorescence (au)')
pbaspect([4 1 1])

%saveas(single_traces_fig,[FigurePath 'figure_ON_single_trace.pdf'])

% %% calculate silenced duration vs response time
% 
% response_time_final = response_time_full(data_filter_full==1);
% silence_time_final = silence_time_full(data_filter_full==1);
% 
% x = silence_time_final;
% y = response_time_final;
% filter = (x<analysis_range) & (y<analysis_range);
% 
% xFit = x(filter);
% yFit = y(filter);
% 
% %mdl = fitlm(xFit,yFit,'RobustOpts','on');
% %mdl = fitlm(xFit,yFit);
% mdl = fitlm(x,y)
% p = coefTest(mdl)
% 
% b = mdl.Coefficients{1,1};
% k = mdl.Coefficients{2,1};
% %options = fitoptions('poly1');
% %options.Upper = [Inf -3];
% %mdl = fit(x(filter)',y(filter)','poly1',options);
% %k = mdl.p1;
% %b = mdl.p2;
% 
% regMdl = regARIMA(0,0,0);
% regMdl.Distribution = struct('Name','t','DoF',3);
% %regMdl is a regARIMA model object. It is a template for estimation.
% 
% %Estimate the regression model with ARIMA errors. Plot the regression line.
% estRegMdl = estimate(regMdl,yFit','X',xFit');
% 
% 
% xRange = linspace(0,20,1000);
% yResult = k*xRange + b;
% 
% memory_fig = figure;
% %scatter(x(filter),y(filter),50,'filled','MarkerFaceColor',[234 194 100]/255,'MarkerEdgeColor',[0 .5 .5],'LineWidth',0.5)
% scatter(x(filter),y(filter),50,'filled','MarkerFaceColor',[115 142 193]/255,'MarkerEdgeColor',[0 .5 .5],'LineWidth',0.5)
% hold on
% plot(xRange,yResult,'-','LineWidth',2)
% %plot(x,mdl)
% xlabel('silenced duration (min) before illumination');
% ylabel('response time (min)');
% xlim([0 analysis_range])
% ylim([0 analysis_range])
% pbaspect([3 2 1])
% %saveas(memory_fig,[FigurePath 'figure_memory.pdf'])


%% fit gamma function

%a = response_time_final(response_time_final<=6);
response_time_final = response_time_full(data_filter_full==1);
%silence_time_final = silence_time_full(data_filter_full==1);
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
ylim([0 0.275])
xlabel('response time (min)')
ylabel('probability')
pbaspect([3 2 1])

%saveas(response_time_fig,[FigurePath 'figure_response_time_hist.pdf'])




%{
%% Bursting data analysis

% load inference results
load([resultsRoot filesep 'cpHMM_results' filesep 'hmm_input_output.mat']);

% color to be used
k_green = brighten([38 142 75]/256,.4);
color_green = [38 143 75]/256; % color from Jake
mRNA_red = brighten([212 100 39]/256,.2);

% set parameters
ap_lim = 0.02;
% histogram
binNum = 17;
binMax = 12;
edges = linspace(0,binMax,binNum);
analysis_range = 8;

response_time = [];

for i = 1:length(hmm_input_output)

    time = hmm_input_output(i).time/60;
    pState = hmm_input_output(i).promoter_state;
    fluo = hmm_input_output(i).fluo;
    predicted_fluo = hmm_input_output(i).predicted_fluo;
    embryo_num = fix(hmm_input_output(i).ncID(1));
    
    ap_vec = hmm_input_output(i).APPosParticle;
    mean_ap = nanmean(ap_vec);
    
    frame_before = find(time <= embryo(embryo_num).time_on);
    frame_after = find(time>embryo(embryo_num).time_on);
    
    if (mean_ap > -ap_lim) && (mean_ap < ap_lim)
    
        if ~isempty(frame_before) && (pState(max(frame_before))>1)
            last_on = intersect(find(pState>1,1,'last'),frame_after);
            response_time_temp = time(last_on)-embryo(embryo_num).time_on + 1/3;
            response_time = [response_time response_time_temp];
        end
        
    end
end


% fit gamma function

a = response_time(response_time<=8);

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
plot(x,y1,'LineWidth',3,'Color',mRNA_red)
mean(a)
mean(response_time)

xlim([0 analysis_range])
xlabel('response time (min)')
ylabel('probability')
pbaspect([3 2 1])

saveas(response_time_fig,[FigurePath 'figure_repression_response_time_hist.pdf'])
%}

% TraceFig = figure;
% yyaxis left
% stairs(time,pState-1);
% ylim([-0.1 2.1])
% yyaxis right
% plot(time,Fluo)
% hold on
% plot(time,predFluo)



%% calculate silenced duration vs response time
% 
% response_time_final = response_time_full(data_filter_full==1);
% silence_time_final = silence_time_full(data_filter_full==1);
% 
% 
% x = silence_time_final;
% y = response_time_final;
% filter = (x<analysis_range) & (y<analysis_range);
% 
% xFit = x(filter);
% yFit = y(filter);
% 
% mdl = fitlm(xFit,yFit,'RobustOpts','on');
% %mdl = fitlm(x(filter),y(filter));
% b = mdl.Coefficients{1,1};
% k = mdl.Coefficients{2,1};
% %options = fitoptions('poly1');
% %options.Upper = [Inf -3];
% %mdl = fit(x(filter)',y(filter)','poly1',options);
% %k = mdl.p1;
% %b = mdl.p2;
% 
% regMdl = regARIMA(0,0,0);
% regMdl.Distribution = struct('Name','t','DoF',3);
% %regMdl is a regARIMA model object. It is a template for estimation.
% 
% %Estimate the regression model with ARIMA errors. Plot the regression line.
% estRegMdl = estimate(regMdl,yFit','X',xFit');
% 
% 
% xRange = linspace(0,20,1000);
% yResult = k*xRange + b;
% 
% memory_fig = figure;
% %scatter(x(filter),y(filter),50,'filled','MarkerFaceColor',[234 194 100]/255,'MarkerEdgeColor',[0 .5 .5],'LineWidth',0.5)
% scatter(x(filter),y(filter),50,'filled','MarkerFaceColor',[115 142 193]/255,'MarkerEdgeColor',[0 .5 .5],'LineWidth',0.5)
% hold on
% plot(xRange,yResult,'-','LineWidth',2)
% %plot(x,mdl)
% xlabel('silenced duration (min) before illumination');
% ylabel('response time (min)');
% xlim([0 analysis_range])
% ylim([0 analysis_range])
% pbaspect([2 3 1])
% saveas(memory_fig,[FigurePath 'figure_memory.pdf'])
% 
% 
% %% fit gamma function
% 
% a = response_time_final(response_time_final<=8);
% 
% % fit gamma function
% [muhat,muci] = mle(a,'distribution','gamma'); % Generic function
% %[muhat,muci] = gamfit(a); % Distribution specific function
% 
% x = 0:0.1:20;
% y1 = gampdf(x,muhat(1),muhat(2))/(binNum/binMax);
% %y2 = gamcdf(x,muhat(1),muhat(2));
% 
% response_time_fig = figure;
% hold on
% h = histogram(a,edges,'Normalization','probability');
% h.FaceColor = mRNA_red;
% plot(x,y1,'LineWidth',3,'Color',mRNA_red)
% mean(a)
% 
% xlim([0 analysis_range])
% xlabel('response time (min)')
% ylabel('probability')
% pbaspect([3 2 1])
% 
% saveas(response_time_fig,[FigurePath 'figure_response_time_hist.pdf'])