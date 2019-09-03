
clc
clearvars

%% load the necessary data

config_file='config_template_rb2_oneHemi_analyzeGridRuns.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

addpath(genpath(pwd))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_kResults.mat' ] ;
loadKResults = load(loadName) ;

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_baseRes.mat' ] ;
load(loadName) ;

% inline func useful
getLogEvid = @(x) x.Para.LogEvidence ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get the results at best k

% kResults = loadKResults.kResults ;
% 
% bestKRes = kResults{baseRes.wsbm.bestKind} ;
% bestKLogEvid = cellfun(getLogEvid,bestKRes) ;
% 
% bestKComs = cell2mat(cellfun(@(x)wsbm_community_assign(x),bestKRes,'UniformOutput',0));
% bestKComs = reshape(bestKComs,NUM_NODES, length(bestKRes) ) ;
% 
% bestKpdist = partition_distance(bestKComs) ;
% bestKdistances = sum(bestKpdist,2) ;

%% get other k

mk_ind = baseRes.wsbm.bestKind - 1;
pk_ind = baseRes.wsbm.bestKind + 1;

% minus k
mKRes = kResults{mk_ind} ;
mKLogEvid = cellfun(getLogEvid,mKRes) ;

mKComs = cell2mat(cellfun(@(x)wsbm_community_assign(x),mKRes,'UniformOutput',0));
mKComs = reshape(mKComs,NUM_NODES, length(mKRes) ) ;

mKpdist = partition_distance(mKComs) ;
mKdistances = sum(mKpdist,2) ;

% plus k
pKRes = kResults{pk_ind} ;
pKLogEvid = cellfun(getLogEvid,pKRes) ;

pKComs = cell2mat(cellfun(@(x)wsbm_community_assign(x),pKRes,'UniformOutput',0));
pKComs = reshape(pKComs,NUM_NODES, length(pKRes) ) ;

pKpdist = partition_distance(pKComs) ;
pKdistances = sum(pKpdist,2) ;

%% run consensus on models that are both central and high log evidences

mK_medianLogEvid = median(mKLogEvid) ;
mK_medianDist = median(mKdistances) ;

pK_medianLogEvid = median(pKLogEvid) ;
pK_medianDist = median(pKdistances) ;

% medianLogEvid = median(bestKLogEvid) ;
% medianDist = median(bestKdistances) ;

% get indicies for top quadrant
mK_optQuadInd = (mKLogEvid > mK_medianLogEvid) & (mKdistances < mK_medianDist) ;  
pK_optQuadInd = (pKLogEvid > pK_medianLogEvid) & (pKdistances < pK_medianDist) ;  

% get the top quadrant models
mk_optQuadModels = mKRes(mK_optQuadInd) ;
pk_optQuadModels = pKRes(pK_optQuadInd) ;

% put the data back in here
for idx = 1:length(mk_optQuadModels)
    mk_optQuadModels{idx}.Data.Raw_Data = baseRes.rawData  ;
    mk_optQuadModels{idx}.Data.n = size(baseRes.rawData ,1) ;
    
    pk_optQuadModels{idx}.Data.Raw_Data = baseRes.rawData  ;
    pk_optQuadModels{idx}.Data.n = size(baseRes.rawData ,1) ;
end

% % viz it
% figure
% scatter(bestKLogEvid,bestKdistances,[],optQuadInd)
% colormap(brewermap(2,'Spectral')) ;
% ylabel('VI Distance') 
% xlabel('Model Log Evidence')
% set(gca,'FontSize',14)

% % lets find centroid of these top models
% centralModel = wsbm_cent_mod(optQuadModels) ;

%% get the consensus models

% mk
mk_optQuadModelsStruct = cell2struct(mk_optQuadModels,{'Model'},2) ;

conReps = 5 ;
mk_conModels = cell(conReps,1) ;
mk_conModelsInfo = cell(conReps,1) ;

% time it
ttt = tic ;
for idx = 1:conReps

    [ mk_conModels{idx} , mk_conModelsInfo{idx} ] = wsbm_consensus_model(mk_optQuadModelsStruct) ;

end
elapsedTime = toc(ttt) ;

% pk
pk_optQuadModelsStruct = cell2struct(pk_optQuadModels,{'Model'},2) ;

conReps = 5 ;
pk_conModels = cell(conReps,1) ;
pk_conModelsInfo = cell(conReps,1) ;

% time it
ttt = tic ;
for idx = 1:conReps

    [ pk_conModels{idx} , pk_conModelsInfo{idx} ] = wsbm_consensus_model(pk_optQuadModelsStruct) ;

end
elapsedTime = toc(ttt) ;

%% how similar are all these conModels

% inputModels = mk_conModels ;
% 
% consensusComs = cell2mat(cellfun(@(x)wsbm_community_assign(x), inputModels ,'UniformOutput',0));
% consensusComs = reshape(consensusComs,NUM_NODES,length(inputModels)) ;
% 
% [ viDistMat , nmiDistMat ] = partition_distance(consensusComs) ;

mk_consensusCent = wsbm_cent_mod(mk_conModels);
pk_consensusCent = wsbm_cent_mod(pk_conModels);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save it

saveName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_mKpK_consensusRuns.mat' ] ;
save(saveName,'mk_consensusCent','pk_consensusCent',...
    '-v7.3') 




