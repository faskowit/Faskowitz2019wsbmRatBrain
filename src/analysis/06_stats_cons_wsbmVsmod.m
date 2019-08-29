
clc
clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load

config_file='config_template_rb2_analyzeGridRuns.m';
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
%% setup 

dat = baseRes.rawData ;
dat(isnan(dat)) = 0 ; 

numComms = size(max(cons_ca.mod),1) ;
numNodes = length(cons_ca.mod) ;

numPerms = 10000 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% comm vs comm stats

rng(123)

% show that the two community structures different
[ emp_vi , emp_nmi] = partition_distance(cons_ca.wsbm,cons_ca.mod) ;

% but different than completely random
viPermMod = zeros(numPerms,1) ;
viPermWsbm = zeros(numPerms,1) ;
nmiPermMod = zeros(numPerms,1) ;
nmiPermWsbm = zeros(numPerms,1) ;

for idx = 1:numPerms

    % the rand cut seems to be biased, based on initial node orderings
    % maybe something to explor in future... not here...
    
%     [~,tmpModSort] = sort(cons_ca.mod) ; 
%     [~,tmpWsbmSort] = sort(cons_ca.wsbm) ;
    
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

[~,wsbmBl ] = get_block_mat(baseRes.rawData,cons_ca.wsbm) ;
[~,modBl ] = get_block_mat(baseRes.rawData,cons_ca.mod) ;

emp_wsbmOD = sum(diag(wsbmBl)) ;
emp_modOD = sum(diag(modBl)) ;
emp_wsbmODratio = emp_wsbmOD ./ ( sum(sum(wsbmBl)) - emp_wsbmOD ) ;
emp_modODratio = emp_modOD ./ ( sum(sum(modBl)) - emp_modOD ) ;
emp_diffODratio = emp_modODratio - emp_wsbmODratio ;

null_wsbmODratio = zeros(numPerms,1) ;
null_modODratio = zeros(numPerms,1) ;

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
    
end

null_diffODratio = null_modODratio - null_wsbmODratio;

% histogram(null_diffODratio)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% a more conservative permutation?

rng(123)

numPerms = 10000 ;

[~,wsbmBl ] = get_block_mat(baseRes.rawData,cons_ca.wsbm) ;
[~,modBl ] = get_block_mat(baseRes.rawData,cons_ca.mod) ;

thr_vals = [ 0.25 0.2 0.15 0.1 0.05 0.01 ] ;

null_wsbmODratio_cell = cell(length(thr_vals),1) ;
null_modODratio_cell = cell(length(thr_vals),1) ;
null_diffODratio_cell = cell(length(thr_vals),1) ;

for pp = 1:length(thr_vals) 

null_wsbmODratio_cell{pp} = zeros(numPerms,1) ;
null_modODratio_cell{pp} = zeros(numPerms,1) ;

for idx = 1:numPerms
    
    if mod(idx,500) == 0 
        disp(idx)
    end
        
    randDat = randmio_dir_connected(dat,thr_vals(pp)) ;
    
    [~,tmpWsbmBl] = get_block_mat(randDat,cons_ca.wsbm) ;
    [~,tmpModBl] = get_block_mat(randDat,cons_ca.mod) ;
    
    tmpWsbmOD = sum(diag(tmpWsbmBl)) ;
    tmpModOD = sum(diag(tmpModBl)) ;
    
    null_wsbmODratio_cell{pp}(idx) = tmpWsbmOD ./ ( sum(sum(tmpWsbmBl)) - tmpWsbmOD ) ;
    null_modODratio_cell{pp}(idx) = tmpModOD ./ ( sum(sum(tmpModBl)) - tmpModOD ) ;
    
end

null_diffODratio_cell{pp} =  null_modODratio_cell{pp} - null_wsbmODratio_cell{pp} ;

end

%% vis it

for pp = 1:length(thr_vals) 

    histogram(null_diffODratio_cell{pp})
    hold on
    
end

% sum(null_diffODratio_cell{5} > emp_diffODratio)
% 
% ans =
% 
%      4
%
%
% 5 / 10001
% 
% ans =
% 
%    4.9995e-04
% 

% for both hemi--> all 0's except for last one
% sum(null_diffODratio_cell{6} > emp_diffODratio)
% 
% ans =
% 
%    423

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% confirm that wsbm still more random than chance

emp_wsbm_mod = modularity_q(dat,cons_ca.wsbm) ;

randq = zeros(numPerms,1) ;
randq2 = zeros(numPerms,1) ;

for idx = 1:numPerms
    
    permvec1 = randperm(numNodes) ;
    randq(idx) = modularity_q(dat,cons_ca.wsbm(permvec1)) ;
    randq2(idx) = modularity_q(randmio_dir_connected(dat,0.05),cons_ca.wsbm) ;
end

% yes, waaaaay passes bootstrap.


