
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

currCAs = baseRes.ca_K{9} ;

dat = baseRes.rawData ;
dat(isnan(dat)) = 0 ;

num=size(currCAs,2) ;
quickQ = zeros([num 1]) ;
for idx = 1:num
   quickQ(idx) = modularity_q(dat,currCAs(:,idx)) ;
end

viMat = partition_distance(currCAs) ;
h = plotSpace2(cca(viMat,2,1000),quickQ) ;

%%%%%%%%%%%%%%%%%%

currCAs = baseRes.caMod_K{9} ;

dat = baseRes.rawData ;
dat(isnan(dat)) = 0 ;

num=size(currCAs,2) ;
quickQ = zeros([num 1]) ;
for idx = 1:num
   quickQ(idx) = modularity_q(dat,currCAs(:,idx)) ;
end

viMat = partition_distance(currCAs) ;
h = plotSpace(cca(viMat,2,1000),quickQ) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% viz with the log evidence 

currBaseResInd = 9 ;
currCAs = baseRes.ca_K{currBaseResInd} ;
viMat = partition_distance(currCAs) ;
h = plotSpace2(cca(viMat,2,1000),baseRes.logEvid_K{currBaseResInd}) ;

%%%%%%%%%%%%%%%%%%

currCAs = baseRes.caMod_K{9} ;

dat = baseRes.rawData ;
dat(isnan(dat)) = 0 ;

num=size(currCAs,2) ;
quickQ = zeros([num 1]) ;
for idx = 1:num
   quickQ(idx) = modularity_q(dat,currCAs(:,idx)) ;
end

viMat = partition_distance(currCAs) ;
h = plotSpace(cca(viMat,2,1000),quickQ) ;
