ResultsPaths = {'hbBAC-MS2-27_5C\cpHMM_results\compiledResults_w7_K3_p0_ap9_t8_f2D.mat',...
    'hbBAC-MS2-25C\cpHMM_results\compiledResults_w7_K3_p0_ap9_t8_f2D.mat',...
    'hbBAC-MS2-22_5C\cpHMM_results\compiledResults_w7_K3_p0_ap9_t8_f2D.mat',...
    'hbBAC-MS2-20C\cpHMM_results\compiledResults_w7_K3_p0_ap9_t8_f2D.mat',...
    'hbBAC-MS2-17_5C\cpHMM_results\compiledResults_w7_K3_p0_ap9_t8_f2D.mat'};
outdir = 'S:\Gabriella\Dropbox\ProcessedEnrichmentData\GMPlots20210304\hbBAC-MS2-AllTemps\V1\';

GenerateStandardPlots(ResultsPaths)
GenerateAPSubplots(ResultsPaths, outdir)
GenerateProjectSubplots(ResultsPaths, outdir)
%%

ResultsPaths2 = {'hbBAC-MS2-27_5C\cpHMM_results\compiledResults_w7_K3_p0_ap9_t8_f2D.mat',...
    'hbBAC-MS2-25C\cpHMM_results\compiledResults_w7_K3_p0_ap9_t8_f2D.mat',...
    'hbBAC-MS2-22_5C\cpHMM_results\compiledResults_w7_K3_p0_ap9_t10_f2D.mat',...
    'hbBAC-MS2-20C\cpHMM_results\compiledResults_w7_K3_p0_ap9_t13_f2D.mat',...
    'hbBAC-MS2-17_5C\cpHMM_results\compiledResults_w7_K3_p0_ap9_t18_f2D.mat'};
outdir2 = 'S:\Gabriella\Dropbox\ProcessedEnrichmentData\GMPlots20210304\hbBAC-MS2-AllTemps\V2\';

GenerateStandardPlots(ResultsPaths2)
GenerateAPSubplots(ResultsPaths2, outdir2)
GenerateProjectSubplots(ResultsPaths2, outdir2)

%%

ResultsPaths3 = { 'hbBAC-MS2-25C\cpHMM_results\compiledResults_w7_K3_p0_ap9_t8_f2D.mat',...
    'HbMS2JB3\cpHMM_results\compiledResults_w7_K3_p0_ap9_t8_f2D.mat'};
outdir3 = 'S:\Gabriella\Dropbox\ProcessedEnrichmentData\GMPlots20210304\hbBAC-MS2-AllTemps\V3\';
GenerateStandardPlots(ResultsPaths3)
GenerateTwoProjectComps(ResultsPaths3, outdir3)

%%


ResultsPaths4 = { 'hbBAC-MS2-25C\cpHMM_results\compiledResults_w7_K3_p0_ap9_t8_f2D.mat',...
    'hbBAC-MS2-25C\cpHMM_results\compiledResults_w6_K3_p0_ap9_t8_f2D.mat'};
outdir4 = 'S:\Gabriella\Dropbox\ProcessedEnrichmentData\GMPlots20210304\hbBAC-MS2-AllTemps\V4\';

GenerateStandardPlots(ResultsPaths4)
GenerateTwoProjectComps(ResultsPaths4, outdir4)


%%


ResultsPaths5 = { 'hbBAC-MS2-20C\cpHMM_results\compiledResults_w7_K3_p0_ap9_t13_f2D.mat',...
    'hbBAC-MS2-20C\cpHMM_results\compiledResults_w9_K3_p0_ap9_t13_f2D.mat'};
outdir5 = 'S:\Gabriella\Dropbox\ProcessedEnrichmentData\GMPlots20210304\hbBAC-MS2-AllTemps\V5\';

GenerateStandardPlots(ResultsPaths5)
GenerateTwoProjectComps(ResultsPaths5, outdir5)

%%

ResultsPaths6 = { 'hbBAC-MS2-25C\cpHMM_results\compiledResults_w7_K3_p0_ap9_t8_f2D.mat',...
    'hbBAC-MS2-25C\cpHMM_results\compiledResults_w4_K3_p0_ap9_t8_f2D.mat'};
outdir6 = 'S:\Gabriella\Dropbox\ProcessedEnrichmentData\GMPlots20210304\hbBAC-MS2-AllTemps\V6\';

GenerateStandardPlots(ResultsPaths6)
GenerateTwoProjectComps(ResultsPaths6, outdir6)


%%

ResultsPaths7 = { 'hbBAC-MS2-25C\cpHMM_results\compiledResults_w7_K3_p0_ap9_t8_f2D.mat',...
    'hbBAC-MS2-25C\cpHMM_results\compiledResults_w6_K3_p0_ap9_t8_f2D.mat',...
    'hbBAC-MS2-25C\cpHMM_results\compiledResults_w5_K3_p0_ap9_t8_f2D.mat',...
    'hbBAC-MS2-25C\cpHMM_results\compiledResults_w4_K3_p0_ap9_t8_f2D.mat'};
outdir7 = 'S:\Gabriella\Dropbox\ProcessedEnrichmentData\GMPlots20210304\hbBAC-MS2-25C-allweights\V7\';

GenerateStandardPlots(ResultsPaths7, outdir7)
GenerateProjectSubplots(ResultsPaths7, outdir7)

%%

ResultsPaths7 = {'hbBAC-MS2-27_5C\NC13\cpHMM_results\compiledResults_w5_K3_p0_ap9_t1_f2D.mat',...
    'hbBAC-MS2-25C\NC13\cpHMM_results\compiledResults_w5_K3_p0_ap9_t1_f2D.mat',...
    'hbBAC-MS2-22_5C\NC13\cpHMM_results\compiledResults_w5_K3_p0_ap9_t1_f2D.mat',...
    'hbBAC-MS2-20C\NC13\cpHMM_results\compiledResults_w5_K3_p0_ap9_t1_f2D.mat',...
    'hbBAC-MS2-17_5C\NC13\cpHMM_results\compiledResults_w5_K3_p0_ap9_t1_f2D.mat'};
outdir7 = 'S:\Gabriella\Dropbox\ProcessedEnrichmentData\GMPlots20210304\hbBAC-MS2-AllTemps\NC13\';

GenerateStandardPlots(ResultsPaths7)
GenerateAPSubplots(ResultsPaths7, outdir7)
GenerateProjectSubplots(ResultsPaths7, outdir7)

%%

ResultsPaths8 = { 'hbBAC-MS2-20C\cpHMM_results\compiledResults_w7_K3_p0_ap9_t13_f2D.mat',...
    'hbBAC-MS2-20C\cpHMM_results\compiledResults_w9_K3_p0_ap9_t13_f2D.mat'};
outdir8 = 'S:\Gabriella\Dropbox\ProcessedEnrichmentData\GMPlots20210304\hbBAC-MS2-AllTemps\V8\';

GenerateStandardPlots(ResultsPaths8)
GenerateTwoProjectComps(ResultsPaths8, outdir8)

%%


ResultsPaths9 = {'hbBAC-MS2-27_5C\NC13\cpHMM_results\compiledResults_w7_K3_p0_ap9_t1_f2D.mat',...
    'hbBAC-MS2-25C\NC13\cpHMM_results\compiledResults_w7_K3_p0_ap9_t1_f2D.mat',...
    'hbBAC-MS2-22_5C\NC13\cpHMM_results\compiledResults_w7_K3_p0_ap9_t1_f2D.mat',...
    'hbBAC-MS2-20C\NC13\cpHMM_results\compiledResults_w7_K3_p0_ap9_t1_f2D.mat',...
    'hbBAC-MS2-17_5C\NC13\cpHMM_results\compiledResults_w7_K3_p0_ap9_t1_f2D.mat'};
outdir9 = 'S:\Gabriella\Dropbox\ProcessedEnrichmentData\GMPlots20210304\hbBAC-MS2-AllTemps\NC13_w7\';

GenerateStandardPlots(ResultsPaths9, outdir9)
GenerateAPSubplots(ResultsPaths9, outdir9)
GenerateProjectSubplots(ResultsPaths9, outdir9)


%%

ResultsPaths = {'hbBAC-MS2-27_5C\cpHMM_results\compiledResults_w7_K3_p0_ap9_t8_f2D.mat',...
    'hbBAC-MS2-27_5C\NC13\cpHMM_results\compiledResults_w7_K3_p0_ap9_t1_f2D.mat',...
    'hbBAC-MS2-27_5C\cpHMM_results\compiledResults_w4_K3_p0_ap9_t8_f2D_dt40.mat',...
    'hbBAC-MS2-27_5C\cpHMM_results\compiledResults_w5_K3_p0_ap9_t8_f2D_dt60.mat'};



