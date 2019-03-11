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
%% figure out optimal k

rng(123)

confInterval = [ 0.5 99.5 ] ;
numK = length(COM_NUM_RANGE) ;
robDiffMeansLogEvid = zeros(numK,2) ;

diffMeans = @(x,y) mean(x) - mean(y) ;

for idx = 1:(numK-1)
      
    logE1 = baseRes.wsbm.logEvid_K{idx} ;
    logE2 = baseRes.wsbm.logEvid_K{idx+1} ;
      
    if length(logE1) > length(logE2)
        % randomly select less models
        targSize = length(logE2) ;
        
        logE1 = datasample(logE1,targSize,'Replace',false) ;
    end
    
    robDiffMeansLogEvid(idx,:) = ...
        prctile(bootstrp(10000,diffMeans, ...
        logE1, logE2 ) ...
        ,confInterval) ;
    
end

% get the level at which difference in means does not credibly dip
tmpDiffs = find(robDiffMeansLogEvid(:,2) > 0) ;
ind = tmpDiffs(1) ;
if ~ind 
    error('no valid inds')
end
bestK = COM_NUM_RANGE(ind) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save baseRes again, with bestK

baseRes.wsbm.bestK = bestK ;
baseRes.wsbm.bestKind = ind ;

saveName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_baseRes.mat' ] ;
% loads a struct named 'baseRes'
save(saveName,'baseRes') ;






