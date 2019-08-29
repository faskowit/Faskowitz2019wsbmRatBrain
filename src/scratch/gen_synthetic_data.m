%% clear stuff

clc
clearvars

%% load the necessary data

config_file='config_template_rb2_analyzeGridRuns.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

%% load data

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_consensusRuns.mat' ] ;
lll = load(loadName) ;

templateModel = lll.consensusCent ; 

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_baseRes.mat' ] ;
lll = load(loadName) ;

ca_mod = lll.baseRes.mod.bestKcentmod_ca ;

%% actual data

ca_wsbm = wsbm_community_assign(templateModel) ;

templateAdj = baseRes.Raw_Data ; 
templateAdj(~~isnan(templateAdj)) = 0 ;

%% and the modular model 

muMod = dummyvar(ca_mod)' ;
[~,modularityModel] = wsbm(templateModel.Data.Raw_Data, ...
    size(muMod,1), ...
    'W_Distr', templateModel.W_Distr, ...
    'E_Distr', templateModel.E_Distr, ...
    'alpha', templateModel.Options.alpha, ...
    'mu_0', muMod , ...
    'verbosity', 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setup the generation

nNodes = size(templateAdj,1) ;
numPerms = 10000 ;

wsbm_synthmats = zeros([nNodes nNodes numPerms]) ;
mod_synthmats = zeros([nNodes nNodes numPerms]) ;

for idx = 1:numPerms
   
    disp(idx)
    wsbm_synthmats(:,:,idx) = wsbm_synth_adj_gen(templateModel,0) ;
    mod_synthmats(:,:,idx) = wsbm_synth_adj_gen(modularityModel,0) ;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save the data

saveName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_synmats.mat' ] ;
save(saveName,'wsbm_synthmats','mod_synthmats')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get strenght conf. int. 

wsbm_strengths = zeros([nNodes numPerms]) ;
mod_strengths = zeros([nNodes numPerms]) ;

for idx = 1:numPerms
   
    disp(idx)
    tmp = wsbm_synthmats(:,:,idx) ;
    tmp(isnan(tmp)) = 0 ;
    wsbm_strengths(:,idx) = strengths_dir(tmp) ;
    
    tmp = mod_synthmats(:,:,idx) ;
    tmp(isnan(tmp)) = 0 ;   
    mod_strengths(:,idx) = strengths_dir(tmp) ;
    
end

%% setup vars for perm tests
% select some metrics that can be reduce to scalar, to get distributions
% across many permutations

measures = { 'assortO' 'assortI' 'partiO' 'partiI' 'eff' 'trans' 'cc' 'dens' } ;
models = { 'wsbm' 'mod' 'rand' } ;

% defined above
% numPerms = 10000 ;
nNodes = templateModel.Data.n ;

permStruct = struct() ;
modelStruct = struct() ;

for idx = 1:3
    for jdx = 1:8
        permsStruct.(models{idx}).(measures{jdx}) = zeros([ numPerms 1 ]); 
    end
end

modelStruct.wsbm = templateModel ;
modelStruct.mod = modularityModel ;

genWeightCutoff = 10 ;

%% iterate through permuations

for idx=1:numPerms
   
    disp(idx)
    
    [~,tmpAdj] = wsbm_synth_adj_gen(modelStruct.wsbm,0);
    tmpAdj(1:nNodes+1:end)=0; %clear diagonal
    tmpAdj(tmpAdj > genWeightCutoff) = genWeightCutoff ;
    tmpAdj = double(tmpAdj);
    
    permsStruct.wsbm.(measures{1})(idx) = assortativity_wei(tmpAdj,3) ;
    permsStruct.wsbm.(measures{2})(idx) = assortativity_wei(tmpAdj,4) ;
    permsStruct.wsbm.(measures{3})(idx) = median(participation_coef(tmpAdj,ca_wsbm,1));
    permsStruct.wsbm.(measures{4})(idx) = median(participation_coef(tmpAdj,ca_wsbm,2));
    try
        permsStruct.wsbm.(measures{5})(idx) = efficiency_wei(tmpAdj);
    catch
        permsStruct.wsbm.(measures{5})(idx) = NaN ;
    end
    permsStruct.wsbm.(measures{6})(idx) = transitivity_wd(tmpAdj);
    permsStruct.wsbm.(measures{7})(idx) = median(clustering_coef_wd(tmpAdj)) ;
    permsStruct.wsbm.(measures{8})(idx) = density_dir(tmpAdj) ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [~,tmpAdj] = wsbm_synth_adj_gen(modelStruct.mod,0);
    tmpAdj(1:nNodes+1:end)=0; %clear diagonal
    tmpAdj(tmpAdj > genWeightCutoff) = genWeightCutoff ;
    tmpAdj = double(tmpAdj);
    
    permsStruct.mod.(measures{1})(idx) = assortativity_wei(tmpAdj,3) ;
    permsStruct.mod.(measures{2})(idx) = assortativity_wei(tmpAdj,4) ;
    permsStruct.mod.(measures{3})(idx) = median(participation_coef(tmpAdj,ca_mod,1));
    permsStruct.mod.(measures{4})(idx) = median(participation_coef(tmpAdj,ca_mod,2));
    try
        permsStruct.mod.(measures{5})(idx) = efficiency_wei(tmpAdj);
    catch
        permsStruct.mod.(measures{5})(idx) = NaN ;
    end
    permsStruct.mod.(measures{6})(idx) = transitivity_wd(tmpAdj);
    permsStruct.mod.(measures{7})(idx) = median(clustering_coef_wd(tmpAdj)) ;
    permsStruct.mod.(measures{8})(idx) = density_dir(tmpAdj) ;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % randomize template
    randModel = wsbm_randomize_model_params(modelStruct.wsbm,3);
    [~,tmpAdj] = wsbm_synth_adj_gen(randModel,0);
    tmpAdj(1:nNodes+1:end)=0; %clear diagonal
    tmpAdj(tmpAdj > genWeightCutoff) = genWeightCutoff ;
    tmpAdj = double(tmpAdj);

    permsStruct.rand.(measures{1})(idx) = assortativity_wei(tmpAdj,3) ;
    permsStruct.rand.(measures{2})(idx) = assortativity_wei(tmpAdj,4) ;
    permsStruct.rand.(measures{3})(idx) = median(participation_coef(tmpAdj,ca_wsbm,1));
    permsStruct.rand.(measures{4})(idx) = median(participation_coef(tmpAdj,ca_wsbm,2));
    try
        permsStruct.rand.(measures{5})(idx) = efficiency_wei(tmpAdj);
    catch
        permsStruct.rand.(measures{5})(idx) = NaN ;
    end
    permsStruct.rand.(measures{6})(idx) = transitivity_wd(tmpAdj);
    permsStruct.rand.(measures{7})(idx) = median(clustering_coef_wd(tmpAdj)) ;
    permsStruct.rand.(measures{8})(idx) = density_dir(tmpAdj) ;

end

%% empirical yo

permsStruct.empir.(measures{1}) = assortativity_wei(templateAdj,3) ;
permsStruct.empir.(measures{2}) = assortativity_wei(templateAdj,4) ;
permsStruct.empir.(measures{3}) = median(participation_coef(templateAdj,ca_wsbm,1));
permsStruct.empir.(measures{4}) = median(participation_coef(templateAdj,ca_wsbm,2));
try
    permsStruct.empir.(measures{5}) = efficiency_wei(templateAdj);
catch
    permsStruct.empir.(measures{5}) = NaN ;
end
permsStruct.empir.(measures{6}) = transitivity_wd(templateAdj);
permsStruct.empir.(measures{7}) = median(clustering_coef_wd(templateAdj)) ;
permsStruct.empir.(measures{8}) = density_dir(templateAdj) ;


%% save the perms structt

saveName = strcat(PROJECT_DIR,'/data/processed/',OUTPUT_STR, '_', GRID_RUN,'_permsStruct.mat') ;
save(saveName,...
        'permsStruct') ;

% load(saveName)
