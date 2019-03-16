%% clear stuff

clc
clearvars

%% load the necessary data

config_file='config_template_rb2_oneHemi_analyzeGridRuns.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

%% load data

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_consensusRuns.mat' ] ;
load(loadName)

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_baseRes.mat' ] ;
load(loadName) ;

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_consensusCAs.mat' ] ;
load(loadName) ;

%% 

templateAdj = baseRes.rawData ;
templateAdj(~~isnan(templateAdj)) = 0 ;

%% and the modular model 

muMod = dummyvar(cons_ca.mod)' ;
[~,modularityModel] = wsbm(consensusCent.Data.Raw_Data, ...
    size(muMod,1), ...
    'W_Distr', consensusCent.W_Distr, ...
    'E_Distr', consensusCent.E_Distr, ...
    'alpha', consensusCent.Options.alpha, ...
    'mu_0', muMod , ...
    'verbosity', 0);

muWsbm = dummyvar(cons_ca.wsbm)' ;
[~,wsbmModel] = wsbm(consensusCent.Data.Raw_Data, ...
    size(muWsbm,1), ...
    'W_Distr', consensusCent.W_Distr, ...
    'E_Distr', consensusCent.E_Distr, ...
    'alpha', consensusCent.Options.alpha, ...
    'mu_0', muWsbm , ...
    'verbosity', 0);

%% eval gen call

numEval = 10000 ;

[eval_wsbm_B,eval_wsbm_E,eval_wsbm_K,eval_wsbm_EMD] = wsbm_eval_model_energy(wsbmModel,numEval,0,0);
[eval_mod_B,eval_mod_E,eval_mod_K,eval_mod_EMD] = wsbm_eval_model_energy(modularityModel,numEval,0,0);

[eval_wsbmRand_B,eval_wsbmRand_E,eval_wsbmRand_K,eval_wsbmRand_EMD] = wsbm_eval_model_energy(wsbmModel,numEval,1,0);
[eval_modRand_B,eval_modRand_E,eval_modRand_K,eval_modRand_EMD] = wsbm_eval_model_energy(modularityModel,numEval,1,0);

saveName = strcat(PROJECT_DIR,'/data/processed/',OUTPUT_STR, '_', GRID_RUN,'_evalGenReps.mat') ;
% load(saveName) ;
save(saveName,...
    'eval_wsbm_K','eval_wsbmRand_K','eval_wsbm_EMD','eval_wsbmRand_EMD',...
    'eval_mod_K','eval_modRand_K','eval_mod_EMD','eval_modRand_EMD' ) ;

