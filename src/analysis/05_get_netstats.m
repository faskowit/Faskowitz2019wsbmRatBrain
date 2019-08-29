%% clear stuff

clc
clearvars

%% load the necessary data

config_file='config_template_rb2_analyzeGridRuns.m';
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
netstats = struct() ;
commNames = { 'wsbm' 'mod' } ;
dat = baseRes.rawData ;
dat(isnan(dat)) = 0 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% community-wise ish

for cn = 1:length(commNames)

netstats.(commNames{cn}).pc.coef = cell(length(COM_NUM_RANGE),1) ;
netstats.(commNames{cn}).pc.rank = cell(length(COM_NUM_RANGE),1) ;

netstats.(commNames{cn}).modz.score = cell(length(COM_NUM_RANGE),1) ;
netstats.(commNames{cn}).modz.rank = cell(length(COM_NUM_RANGE),1) ;

netstats.(commNames{cn}).vers = cell(length(COM_NUM_RANGE),1) ;

for idx = 1:length(COM_NUM_RANGE)

    disp(idx)
    
    tmp_ca = baseRes.(commNames{cn}).ca_K{idx} ;
    
    numRuns = size(tmp_ca,2) ;
    tmp_parti = zeros(NUM_NODES,numRuns) ;
    tmp_modz = zeros(NUM_NODES,numRuns) ;
    
    for jdx = 1:numRuns
    
        % parti coeff
        tmp_parti(:,jdx) = (participation_coef(dat,tmp_ca(:,jdx),1) + ...
                            participation_coef(dat,tmp_ca(:,jdx),2)) ./ 2 ;

        % module zscore
        tmp_modz(:,jdx) = module_degree_zscore(dat,tmp_ca(:,jdx),3) ;
    end
    
    netstats.(commNames{cn}).pc.coef{idx} = tmp_parti ;
    netstats.(commNames{cn}).pc.rank{idx} = tiedrank(tmp_parti) ;
    
    netstats.(commNames{cn}).modz.score{idx} = tmp_modz ;
    netstats.(commNames{cn}).modz.rank{idx} = tiedrank(tmp_modz) ;
    
    % versatility
    netstats.(commNames{cn}).vers{idx} = get_nodal_versatility(tmp_ca)' ;

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% node-wise ish

[~,~,netstats.node.str] = strengths_dir(dat) ;
netstats.node.cc = clustering_coef_wd(dat) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save it

saveName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_netstatsRes.mat' ] ;
% loads a struct named 'baseRes'
save(saveName,...
    'netstats')
