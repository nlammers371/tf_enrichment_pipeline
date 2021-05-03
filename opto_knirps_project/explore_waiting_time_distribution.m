clear
close all
clc

addpath(genpath('./lib'))

%% Initialization

projectName = 'optokni_eve4+6_MCP-GFP_Homo'; 

liveProject = LiveEnrichmentProject(projectName);
resultsRoot = [liveProject.dataPath filesep 'cpHMM_results' filesep];

% load data
load([resultsRoot 'hmm_input_output.mat'])
FigurePath = [liveProject.figurePath 'waiting_time' filesep];
mkdir(FigurePath)

% color to be used
k_green = brighten([38 142 75]/256,.4);
color_green = [38 143 75]/256; % color from Jake
mRNA_red = brighten([212 100 39]/256,.2);

%% Calculate waiting time distribution

on_off_waiting_time_vec = [];
off_on_waiting_time_vec = [];

test = 0;

for i = 1:length(hmm_input_output)
    
    QCflag_temp = hmm_input_output(i).TraceQCFlag;
    promoter_state_temp = hmm_input_output(i).promoter_state;
    promoter_state_binary = promoter_state_temp>1;
    
    if (sum(QCflag_temp) == length(QCflag_temp))
        
        flag = promoter_state_temp(1)>1;
        count = 1;
        
        for j = 2:length(promoter_state_binary)
            
            if (promoter_state_binary(j) ~= promoter_state_binary(j-1))
                if promoter_state_binary(j) == 1
                    off_on_waiting_time_vec = [off_on_waiting_time_vec count];
                    count = 1;
                else
                    if promoter_state_binary(j) == 0
                        on_off_waiting_time_vec = [on_off_waiting_time_vec count];
                        count = 1;
                    end
                end
            else
                count = count + 1;
            end
            
        end
        
    end
   
end

%% Draw histograms

% histogram parameters
binNum = 40;
binMax = 12;
edges = linspace(0,binMax,binNum);
%edges = 0:1/3:8;

x = off_on_waiting_time_vec*1/3;
%a = a(a<2.6);
x_fit = x(x<3.1);

% fit gamma function
[muhat,muci] = mle(x_fit,'distribution','gamma'); % Generic function
%[muhat,muci] = gamfit(a_fit); % Distribution specific function

xRange = 0:0.1:20;
y1 = gampdf(xRange,muhat(1),muhat(2))/(binNum/binMax);
%y2 = gamcdf(x,muhat(1),muhat(2));

fig = figure(1);
hold on
histogram(x,edges, 'Normalization','probability')
h.FaceColor = mRNA_red;
plot(xRange,y1,'LineWidth',3,'Color',mRNA_red)
mean(x_fit)

xlim([0 8])
xlabel('response time (min)')
ylabel('probability')
pbaspect([3 2 1])

%fig = figure(2);
%hold on
%histogram(on_off_waiting_time_vec*1/3,edges, 'Normalization','probability')
%h.FaceColor = mRNA_red;
%xlim([0 8])

%% Try two distribution fit

%ft = fittype('a*gampdf(x,b,c) + (1-a)*gampdf(x,d,e)');
%myfit = fit(x',y',myfittype)

%f = @(x,lambda,a1,b1,a2,b2) lambda*gampdf(x,a1,b1) + (1-lambda)*gampdf(x,a2,b2);

%[a b] = mle(x,'pdf',@(x,lambda,a1,b1,a2,b2) lambda*gampdf(x,a1,b1) + (1-lambda)*gampdf(x,a2,b2),'start',[1 1 1 1 1],'LowerBound',[0.9 1 1 1 1],'UpperBound',[1 1.1 1.1 1.1 1.1]);

