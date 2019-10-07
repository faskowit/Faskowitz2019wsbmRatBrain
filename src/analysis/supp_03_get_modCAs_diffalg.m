%% clear stuff

clc
clearvars

%% load the necessary data

config_file='config_template_rb2_oneHemi_analyzeGridRuns.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

%% load data

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_baseRes.mat' ] ;
load(loadName) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get the modular com struct

rng(123)

dat = baseRes.rawData ;
dat(isnan(dat)) = 0 ;

% initialize a variable
modCas = cell(length(REASONABLE_COM_RANGE_IND),1) ;
modQs = cell(length(REASONABLE_COM_RANGE_IND),1) ;
gammaRange = 0.1:0.005:6 ;

for idx = 1:length(REASONABLE_COM_RANGE_IND)

    disp([ newline newline num2str(idx) newline newline ])
    
    Klevel = REASONABLE_COM_RANGE_IND(idx) ;
    numRuns = size(baseRes.wsbm.ca_K{idx},2) ;
    
    [modCas{idx},modQs{idx}] = sweep_modularity_dir_atK(dat, Klevel,gammaRange,numRuns) ;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% take a look

% get central model
[centCa, viMat] = wsbm_cent_mod(modCas{baseRes.wsbm.bestKind}) ;

[x,y,indsort] = grid_communities(centCa) ;
imagesc(dat(indsort,indsort))
hold on
plot(x,y,'r','linewidth',2)
axis square
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% add to output

baseRes.mod.bestKcentmod_ca = centCa ;
baseRes.mod.ca_K = modCas ;
baseRes.mod.caQs_K = modQs ;

saveName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_altMod_baseRes.mat' ] ;
save(saveName,'baseRes') ;



