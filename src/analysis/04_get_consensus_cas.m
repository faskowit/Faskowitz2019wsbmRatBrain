
clc
clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load

config_file='config_template_rb2_oneHemi_analyzeGridRuns.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_baseRes.mat' ] ;
load(loadName) ;

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_consensusRuns.mat' ] ;
load(loadName) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get the community assignments

nComm = size(consensusCent.Para.mu,1) ;
ca_wsbm_init = wsbm_community_assign(consensusCent) ;

figure
[newinds,breaks] = wsbm_plot_mat(consensusCent, [],@wsbm_reorder_avgbinblock_mod);

% the new wsbm node order
ca_wsbm_ra = lab_reassign(ca_wsbm_init,newinds) ;

ca_mod_init = baseRes.mod.bestKcentmod_ca ;

% the new mod order
ca_mod_ra = hungarianMatch(ca_wsbm_ra,ca_mod_init) ;

% make a struct
cons_ca =  struct() ;
cons_ca.wsbm = ca_wsbm_ra ;
cons_ca.mod = ca_mod_ra ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save it

saveName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_consensusCAs.mat' ] ;
save(saveName,'cons_ca') ;
