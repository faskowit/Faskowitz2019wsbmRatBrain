
clc
clearvars

%% load the necessary data

config_file='config_template_rb2_analyzeGridRuns.m';
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

kResults = loadKResults.kResults ;

bestKRes = kResults{baseRes.wsbm.bestKind} ;
bestKLogEvid = cellfun(getLogEvid,bestKRes) ;

bestKComs = cell2mat(cellfun(@(x)wsbm_community_assign(x),bestKRes,'UniformOutput',0));
bestKComs = reshape(bestKComs,NUM_NODES, length(bestKRes) ) ;

bestKpdist = partition_distance(bestKComs) ;
bestKdistances = sum(bestKpdist,2) ;

%% run consensus on models that are both central and high log evidences

medianLogEvid = median(bestKLogEvid) ;
medianDist = median(bestKdistances) ;

% get indicies for top quadrant
optQuadInd = (bestKLogEvid > medianLogEvid) & (bestKdistances < medianDist) ;  

% get the top quadrant models
optQuadModels = bestKRes(optQuadInd) ;

% put the data back in here
for idx = 1:length(optQuadModels)
    optQuadModels{idx}.Data.Raw_Data = baseRes.rawData  ;
    optQuadModels{idx}.Data.n = size(baseRes.rawData ,1) ;
end

% viz it
figure
scatter(bestKLogEvid,bestKdistances,[],optQuadInd)
colormap(brewermap(2,'Spectral')) ;
ylabel('VI Distance') 
xlabel('Model Log Evidence')
set(gca,'FontSize',14)

% % lets find centroid of these top models
% centralModel = wsbm_cent_mod(optQuadModels) ;

%% get the consensus models

optQuadModelsStruct = cell2struct(optQuadModels,{'Model'},2) ;

conReps = 5 ;
conModels = cell(conReps,1) ;
conModelsInfo = cell(conReps,1) ;

% time it
ttt = tic ;
for idx = 1:conReps

    [ conModels{idx} , conModelsInfo{idx} ] = wsbm_consensus_model(optQuadModelsStruct) ;

end
elapsedTime = toc(ttt) ;

%% how similar are all these conModels

inputModels = conModels ;

consensusComs = cell2mat(cellfun(@(x)wsbm_community_assign(x), inputModels ,'UniformOutput',0));
consensusComs = reshape(consensusComs,NUM_NODES,length(inputModels)) ;

[ viDistMat , nmiDistMat ] = partition_distance(consensusComs) ;

consensusCent = wsbm_cent_mod(conModels);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save it

saveName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_consensusRuns.mat' ] ;
save(saveName,'conModels','conModelsInfo','consensusCent',...
    '-v7.3') 




