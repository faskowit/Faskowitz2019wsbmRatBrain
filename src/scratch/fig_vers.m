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
% loads a struct named 'baseRes'
load(loadName) ;

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_versatilityRes.mat' ] ;
% loads a struct named 'baseRes'
load(loadName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% lets plot

goodRunsRange = REASONABLE_COM_RANGE_IND ;

weidegrees = strengths_dir(dat) ;

% versatility across k
versMatWSBM = cell2mat(versWSBM') ;
sortIndWsbmVers = sortedInd(sum(versMatWSBM,2)) ;
imagesc(tiedrank(versMatWSBM(sortIndWsbmVers,goodRunsRange)))

versMatMOD = cell2mat(versMOD') ;
sortIndModVers = sortedInd(sum(versMatMOD,2)) ;
imagesc(tiedrank(versMatMOD(sortIndWsbmVers,goodRunsRange)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% agreement

agreeWSBM = agreement(baseRes.ca_K{14}) ./ 750 ;

agreeMOD = agreement(baseRes.bestKmodCAs{14}) ./ 750 ;

%%

cor_ver2deg_wsbm = zeros(1,length(goodRunsRange)) ;
cor_ver2deg_mod = zeros(1,length(goodRunsRange)) ;

for idx = goodRunsRange 

    cor_ver2deg_wsbm(idx) = corr(versMatWSBM(:,idx),weidegrees') ;
    cor_ver2deg_mod(idx) = corr(versMatMOD(:,idx),weidegrees') ;
end
