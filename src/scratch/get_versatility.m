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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get versatility per node, at each level

% initialize a variable
vers.wsbm = cell(length(COM_NUM_RANGE),1) ;
vers.mod = cell(length(COM_NUM_RANGE),1) ;

dat = baseRes.rawData ;
dat(isnan(dat)) = 0 ;

weiDegree = strengths_dir(dat) ;

for idx = 1:length(COM_NUM_RANGE)

    tmpWSBMCa = baseRes.wsbm.ca_K{idx} ;
    tmpMODCa = baseRes.mod.ca_K{idx} ;
    
    vers.wsbm{idx} = get_nodal_versatility(tmpWSBMCa)' ;
    vers.mod{idx} = get_nodal_versatility(tmpMODCa)' ;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% lets bootstrap some versatility readings
% turns out that versatility is an unbiased estimator (based on this boot)

% verResuktsWSBM = cell(length(baseRes.ca_K),1) ;
% verResuktsMOD = cell(length(baseRes.ca_K),1) ;
% 
% for idx = 1:length(baseRes.ca_K)
% 
%     numReps = 10000 ;
%     sampleSize = size(baseRes.ca_K{idx},2) ;
% 
%     tmpVersWSBM = zeros(NUM_NODES,numReps) ;
%     tmpVersMOD = zeros(NUM_NODES,numReps) ;
%     
%     for jdx = 1:numReps
% 
%         disp(jdx)
%         randInd = randsample(sampleSize,sampleSize,true) ;
% 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % wsbm
%         
%         tmpCA = baseRes.ca_K{idx}(:,randInd) ;
% 
%         tmpVersWSBM(:,jdx) = get_nodal_versatility(tmpCA) ;
% 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % now modular
%         
%         tmpCA = baseRes.bestKmodCAs{idx}(:,randInd) ;
% 
%         tmpVersMOD(:,jdx) = get_nodal_versatility(tmpCA) ;
%         
%     end
% 
%     verResuktsWSBM{idx} = prctile(tmpVersWSBM,[2.5 97.5],2) ;
%     verResuktsMOD{idx} = prctile(tmpVersMOD,[2.5 97.5],2) ;
%     
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save it

saveName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_versatilityRes.mat' ] ;
% loads a struct named 'baseRes'
save(saveName,...
    'vers')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% lets plot

goodRunsRange = REASONABLE_COM_RANGE_IND ;

weidegrees = strengths_dir(dat) ;

% versatility across k
versMatWSBM = cell2mat(versWSBM') ;
versMatWSBM(isnan(versMatWSBM)) = 0 ;
sortIndWsbmVers = sortedInd(sum(versMatWSBM,2)) ;
% sortIndWsbmWei = sortedInd(sum(weidegrees,1)) ;

imagesc(versMatWSBM(:,goodRunsRange))
imagesc( tiedrank(versMatWSBM(:,goodRunsRange)))

stdWsbmVers = std(versMatWSBM(:,goodRunsRange),[],2) ;

versMatMOD = cell2mat(versMOD') ;
sortIndModVers = sortedInd(sum(versMatMOD,2)) ;

imagesc(versMatMOD(:,goodRunsRange))
imagesc(tiedrank(versMatMOD(:,goodRunsRange)))

stdModVers = std(versMatMOD(sortIndModVers,goodRunsRange),[],2) ;

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






















