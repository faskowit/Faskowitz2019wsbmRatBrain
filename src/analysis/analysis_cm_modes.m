
clc
clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load

config_file='config_template_rb2_oneHemi_analyzeGridRuns.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_baseRes.mat' ] ;
% loads a struct named 'baseRes'
load(loadName) ;

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_motifAna.mat' ] ;
load(loadName) ;

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_consensusCAs.mat' ] ;
load(loadName) ;

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_netstatsRes.mat' ] ;
load(loadName) ;

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_motifEntropy.mat' ] ;
load(loadName) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get data

commNames = { 'wsbm' 'mod' } ;
%modes = struct() ;

nNodes = size(baseRes.rawData,1) ;
nanE = isnan(baseRes.rawData) ;

ind = baseRes.wsbm.bestKind ;
eMNames = { 'odMat' 'aMat' 'cMat' 'pMat' 'dMat' } ;

eMatStruct = struct() ;

for cn = 1:length(commNames)

    eMatStruct.(commNames{cn}).eMat = zeros(nNodes,nNodes,length(eMNames)) ;
    
    for idx = 1:length(eMNames)
        tmp = mean(motifAna.(commNames{cn}).motifEdgeMats{ind}.(eMNames{idx}) > 0,3) ;
        eMatStruct.(commNames{cn}).(eMNames{idx}) = tmp ;
        eMatStruct.(commNames{cn}).eMat(:,:,idx) = tmp ;
    end

    [~,eMatMax] = max(eMatStruct.(commNames{cn}).eMat,[],3) ;
    eMatMax = eMatMax .* ~nanE ;

    eMatStruct.(commNames{cn}).mode = get_block_mode(eMatMax,cons_ca.(commNames{cn}),0) ;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

imsc_grid_comm(iii,cons_ca.wsbm)

iii = iii .* ~nanE ;








