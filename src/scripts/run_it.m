%% config
% first run the config file to get some variables in our environment

% clean enviroment
clc
clearvars

%% config

config_file='config_template.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

% addpath('~/JOSHSTUFF/projects/sbm3/src/tools/')
addpath('~/JOSHSTUFF/projects/blockmodeltools/src/tools/')

%% read data
% load the raw data

rawData = load('data/raw/CCTX_matrix.mat') ;

% make the input matrix
inputMat = rawData.CIJw ;
% inputMat = rawData.CIJ ;
inputMat = weight_conversion(inputMat,'autofix') ;
inputMat(inputMat==0) = NaN ;

inputMat = inputMat(1:77,1:77) ;

%% setup weight / edge pairs to iter over

weiEdgePair = cell(1,1);

% weiEdgePair{1} = [ setup_distr('exp') ...
%     setup_distr('poisson',[0 0.1]) ] ;
% weiEdgePair{2} = [ setup_distr('exp') ...
%     setup_distr('poisson',[0 0.01]) ] ;
% weiEdgePair{3} = [ setup_distr('exp') ...
%     setup_distr('poisson',[0 0.001]) ] ;
% weiEdgePair{4} = [ setup_distr('pareto',[0.001,0],1) ...
%     setup_distr('poisson',[0 0.1]) ] ;
% weiEdgePair{5} = [ setup_distr('pareto',[0.001,0],1) ...
%     setup_distr('poisson',[0 0.01]) ] ;
% weiEdgePair{6} = [ setup_distr('pareto',[0.001,0],1) ...
%     setup_distr('poisson',[0 0.001]) ] ;

weiEdgePair{1} = [ setup_distr(WEIGHT_DIST) ... 
    setup_distr(EDGE_DIST) ] ;

loopResults = cell(length(weiEdgePair),1) ;

%% setup looper
% lognormal, poisson

kIterOver = 7:1:14;    
kLoopIter = 50 ;

for paramIdx = 1:length(weiEdgePair)

    wsbmLooperInputs  = cell(numel(kIterOver),1);

    for idx = 1:numel(kIterOver)

        wsbmLooperInputs{idx} = { kIterOver(idx), ... 
            'W_Distr', weiEdgePair{paramIdx}(1), ...
            'E_Distr', weiEdgePair{paramIdx}(2), ...
            'numTrials',WSBM_NUM_TRIAL,...
            'alpha', 0.5, ...
            'mainMaxIter', WSBM_MAIN_ITER , ...
            'muMaxIter' , WSBM_MU_ITER,...
            'mainTol',0.01, ...
            'muTol', 0.01 ,...
            'verbosity', 0};

    end

    %% store results
    
    loopResults{paramIdx} = struct() ;
    
    [ loopResults{paramIdx}.kLoopRes , loopResults{paramIdx}.kLoopMod ] = wsbm_looper_wrapper(inputMat, ...
        wsbmLooperInputs, ...
        kLoopIter, []) ;

end

%% save the info

outName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR,'_looper_data.mat');
save(outName,'loopResults','-v7.3')

% consensus

consResults = cell(length(weiEdgePair),1) ;

for paramIdx = 1:length(weiEdgePair)

    % setup consensus result as struct
    consResults{paramIdx} = struct();
    
    % find the k at which logEvid is max
    kLoopLogEvidMean = mean(loopResults{paramIdx}.kLoopRes(:,2:end),2) ;
    [~,maxInd] = max(kLoopLogEvidMean) ;
    k = kIterOver(maxInd) ;

    % get all the models with best k
    kLoopBestModels = cell2struct(loopResults{paramIdx}.kLoopMod(maxInd,:), ...
        'Model', kLoopIter );

    % get the consensus
    [ consResults{paramIdx}.templateModel , consResults{paramIdx}.consensusStruct ] = ...
        wsbm_consensus_model(kLoopBestModels) ;

end

%% save the information

outName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR,'_basicData_v7p3.mat');
save(outName,...
    'rawData',...
    'inputMat',...
    'loopResults',...
    'consResults',...
    'config_file',...
    'PROJECT_DIR',...
    '-v7.3')

%% also save smaller file based on my results i want

keepRes = cell(2,1) ;
keepRes{1} = consResults{3} ;
keepRes{2} = consResults{6} ;

outName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR,'_basicDataSmall_v7p3.mat');
save(outName,...
    'rawData',...
    'inputMat',...
    'keepRes',...
    'config_file',...
    'PROJECT_DIR',...
    '-v7.3')







