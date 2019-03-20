
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

% loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_motifAna.mat' ] ;
% load(loadName) ;
% 
loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_consensusCAs.mat' ] ;
load(loadName) ;

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_netstatsRes.mat' ] ;
load(loadName) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% comm vs comm stats

rng(123)

% show that the two community structures different
[ emp_vi , emp_nmi] = partition_distance(cons_ca.wsbm,cons_ca.mod) ;

% but different than completely random
numPerms = 10000 ;
numNodes = length(cons_ca.mod) ;
viPermMod = zeros(numPerms,1) ;
viPermWsbm = zeros(numPerms,1) ;
nmiPermMod = zeros(numPerms,1) ;
nmiPermWsbm = zeros(numPerms,1) ;

for idx = 1:numPerms

    [~,tmpModSort] = sort(cons_ca.mod) ; 
    [~,tmpWsbmSort] = sort(cons_ca.wsbm) ;
    
%     tmpModInd = rand_cut(tmpModSort) ;
%     tmpWsbmInd = rand_cut(tmpWsbmSort) ; 
%     [ viPermMod(idx) , nmiPermMod(idx) ] = partition_distance(cons_ca.wsbm,cons_ca.mod(tmpModInd)) ;
%     [ viPermWsbm(idx), nmiPermWsbm(idx) ] = partition_distance(cons_ca.wsbm(tmpWsbmInd),cons_ca.mod) ;
    [ viPermMod(idx) , nmiPermMod(idx) ] = partition_distance(cons_ca.wsbm,cons_ca.mod(randperm(numNodes))) ;
    [ viPermWsbm(idx), nmiPermWsbm(idx) ] = partition_distance(cons_ca.wsbm(randperm(numNodes)),cons_ca.mod) ;
    
end

% histogram(mean([ viPermMod viPermWsbm],2))
% histogram(mean([ nmiPermMod nmiPermWsbm],2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% evaluate the 'modularity' of each cons comms

rng(123)

dat = baseRes.rawData ;
dat(isnan(dat)) = 0 ; 

[~,wsbmBl ] = get_block_mat(baseRes.rawData,cons_ca.wsbm) ;
[~,modBl ] = get_block_mat(baseRes.rawData,cons_ca.mod) ;

emp_wsbmOD = sum(diag(wsbmBl)) ;
emp_modOD = sum(diag(modBl)) ;
emp_wsbmODratio = emp_wsbmOD ./ ( sum(sum(wsbmBl)) - emp_wsbmOD ) ;
emp_modODratio = emp_modOD ./ ( sum(sum(modBl)) - emp_modOD ) ;

null_wsbmODratio = zeros(numPerms,1) ;
null_modODratio = zeros(numPerms,1) ;
null_diffODratio = zeros(numPerms,1) ;

numComms = size(wsbmBl,1) ;

for idx = 1:numPerms
    
    permvec1 = randperm(numNodes) ;
    permvec2 = randperm(numNodes) ;
%     permvec1 = rand_cut(1:numNodes) ;   
%     permvec2 = rand_cut(1:numNodes)  ;
    
    [~,tmpWsbmBl] = get_block_mat(baseRes.rawData,cons_ca.wsbm(permvec1)) ;
    [~,tmpModBl] = get_block_mat(baseRes.rawData,cons_ca.mod(permvec2)) ;
    
    tmpWsbmOD = sum(diag(tmpWsbmBl)) ;
    tmpModOD = sum(diag(tmpModBl)) ;
    
    null_wsbmODratio(idx) = tmpWsbmOD ./ ( sum(sum(tmpWsbmBl)) - tmpWsbmOD ) ;
    null_modODratio(idx) = tmpModOD ./ ( sum(sum(tmpModBl)) - tmpModOD ) ;

    null_diffODratio(idx) = null_modODratio(idx) - null_wsbmODratio(idx) ;
    
end





