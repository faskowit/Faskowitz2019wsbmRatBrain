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
load(loadName)

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_baseRes.mat' ] ;
load(loadName) ;

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_consensusCAs.mat' ] ;
load(loadName) ;

%% 

templateAdj = baseRes.rawData ;
templateAdj(~~isnan(templateAdj)) = 0 ;

%% and the modular model 

muWsbm = dummyvar(cons_ca.wsbm)' ;
[~,wsbmModel] = wsbm(consensusCent.Data.Raw_Data, ...
    size(muWsbm,1), ...
    'W_Distr', consensusCent.W_Distr, ...
    'E_Distr', consensusCent.E_Distr, ...
    'alpha', consensusCent.Options.alpha, ...
    'mu_0', muWsbm , ...
    'verbosity', 0);

muMod = dummyvar(cons_ca.mod)' ;
[~,modularityModel] = wsbm(consensusCent.Data.Raw_Data, ...
    size(muMod,1), ...
    'W_Distr', consensusCent.W_Distr, ...
    'E_Distr', consensusCent.E_Distr, ...
    'alpha', consensusCent.Options.alpha, ...
    'mu_0', muMod , ...
    'verbosity', 0);

%% eval gen call

numEval = 10000 ;

[~,~,eval_wsbm_K] = wsbm_eval_model_energy(wsbmModel,numEval,0,0);
[~,~,eval_mod_K] = wsbm_eval_model_energy(modularityModel,numEval,0,0);

[~,~,eval_wsbmRand_K] = wsbm_eval_model_energy(wsbmModel,numEval,1,0);
[~,~,eval_modRand_K] = wsbm_eval_model_energy(modularityModel,numEval,1,0);

saveName = strcat(PROJECT_DIR,'/data/processed/',OUTPUT_STR, '_', GRID_RUN,'_evalGenReps.mat') ;
% load(saveName) ;
save(saveName,...
    'eval_wsbm_K','eval_wsbmRand_K',...
    'eval_mod_K','eval_modRand_K') ;

% %%  eval with less...
% 
% rng(123)
% 
% % setup func handles for evaluations
% eFunc = cell(1,1) ;
% 
% % make sure all outputs are column vecs
% eFunc{1} = @(A) sum(A,2) + sum(A,1)';
% eFunc{2} = @(A) clustering_coef_wd(A);
% eFunc{3} = @(A) betweenness_wei(1 ./ A);
% 
% numEval = 10000 ;
% 
% [eval_wsbm_B2,...
%     eval_wsbm_E2,...
%     eval_wsbm_K2,...
%     eval_wsbm_EMD2] = wsbm_eval_model_energy(wsbmModel,numEval,0,0,eFunc);
% [eval_mod_B2,...
%     eval_mod_E2,...
%     eval_mod_K2,...
%     eval_mod_EMD2] = wsbm_eval_model_energy(modularityModel,numEval,0,0,eFunc);
% 
% [eval_wsbmRand_B2,...
%     eval_wsbmRand_E2,...
%     eval_wsbmRand_K2,...
%     eval_wsbmRand_EMD2] = wsbm_eval_model_energy(wsbmModel,numEval,1,0,eFunc);
% [eval_modRand_B2,...
%     eval_modRand_E2,...
%     eval_modRand_K2,...
%     eval_modRand_EMD2] = wsbm_eval_model_energy(modularityModel,numEval,1,0,eFunc);
% 
% saveName = strcat(PROJECT_DIR,'/data/processed/',OUTPUT_STR, '_', GRID_RUN,'_evalGenReps_lessFuncs.mat') ;
% % load(saveName) ;
% save(saveName,...
%     'eval_wsbm_K2','eval_wsbmRand_K2','eval_wsbm_EMD2','eval_wsbmRand_2',...
%     'eval_mod_K2','eval_modRand_K2','eval_mod_EMD2','eval_modRand_EMD2' ) ;

%% quick plot

figure

fontsize=12

numEFunc = length(1:5) ;

% calculate num bins
X = [ mean(eval_wsbm_K,2) ; mean(eval_mod_K,2) ] ;
binSize = 3.5*std(X(:))*numel(X)^(-1/3) ;

% if third needed, we can add
cmap = [0    0.4470    0.7410 ;
    0.8500    0.3250    0.0980 ;
    0.8500    0.3250    0.3920 ];

histogram(mean(eval_wsbm_K(:,1:numEFunc),2),...
    'normalization','probability',...
    'FaceColor',cmap(1,:),'EdgeAlpha',0.01,...
    'BinWidth',binSize) 
hold 
histogram(mean(eval_mod_K(:,1:numEFunc),2),...
    'normalization','probability',...
    'FaceColor',cmap(2,:),'EdgeAlpha',0.01,...
    'BinWidth',binSize)

lg = legend('WSBM','Modular') ;
lg.FontSize = fontsize ;

axis square
xl = xlabel('Mean energy') ;
yl = ylabel('Normalized frequency') ;
xl.FontSize = fontsize ;
yl.FontSize = fontsize ;

%%

figure

numEFunc = length(1:5) ;

curr_wsbm_dat = eval_wsbm_K ;
curr_wsbmRand_dat = eval_wsbmRand_K ;
curr_mod_dat = eval_mod_K ;
curr_modRand_dat = eval_modRand_K ;

set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.9, 0.50]);

% tight subplot
% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right])
sp = tight_subplot(2,numEFunc,[ .035 .04 ],[.12 .033],[.1 .1]);

measure_names = {'Strength' 'BDegree' 'Clustering' 'Betwness' 'BetwnessBin'} ;
measure_short_names = { 's' 'bd' 'c' 'b' 'bb' } ;

% measure_names = {'Strength' 'Clustering' 'Betwness'} ;
% measure_short_names = { 's' 'c' 'b' } ;

x_lim = cell(numEFunc,1) ;
binw = 0.03 ;

for idx = 1:numEFunc
   
    combo = [ curr_wsbmRand_dat(:,idx) ; curr_modRand_dat(:,idx) ; ...
        curr_wsbm_dat(:,idx) ; curr_mod_dat(:,idx) ] ;
    x_lim{idx} = [ 0 (max(combo) .* 1.02) ] ;
    
end

for idx = 1:numEFunc
    
    axes(sp(idx))
    
    histogram(curr_wsbm_dat(:,idx),'BinMethod','sturges',...
        'normalization','probability','EdgeAlpha',0.05,'FaceAlpha',0.6,...
        'FaceColor',cmap(1,:),'BinWidth',binw)
    hold
    histogram(curr_mod_dat(:,idx),'BinMethod','sturges',...
        'normalization','probability','EdgeAlpha',0.05,'FaceAlpha',0.6,...
        'FaceColor',cmap(2,:),'BinWidth',binw)
    axis square
    tl = title(measure_names{idx},'FontWeight','normal') ;
    tl.FontSize = fontsize ;
  
    set(gca,'XTickLabel',[])
    xlim(x_lim{idx}) 

end

for idx = 1:5
   
    axes(sp(idx+numEFunc))
    
    histogram(curr_wsbmRand_dat(:,idx),'normalization','probability',...
        'FaceColor',cmap(1,:),'FaceAlpha',0.25,'EdgeAlpha',0.01,'BinWidth',binw)
    hold
    histogram(curr_modRand_dat(:,idx),'normalization','probability',...
        'FaceColor',cmap(2,:),'FaceAlpha',0.25,'EdgeAlpha',0.01,'BinWidth',binw)
    axis square
    
    % and plot a line for the empirical
    ylimits = ylim ;
    
    tmp = median(mean(curr_wsbm_dat(:,idx)),2) ;
    plot([ tmp tmp ],...
        [ylimits(1) ylimits(2)*0.95],...
        'Color',[cmap(1,:) 0.9],'LineWidth',1.5)
    
    tmp = median(mean(curr_mod_dat(:,idx)),2) ;
    plot([ tmp tmp ],...
        [ylimits(1) ylimits(2)*0.95],...
        'Color',[cmap(2,:) 0.9],'LineWidth',1.5)
    
    ylim(ylimits)
    
    %x_lim{idx} = xlim() ;
    xlim(x_lim{idx}) ;
    
    xl = xlabel(strcat('KS({\it ',measure_short_names{idx},'})'))
    xl.FontSize = fontsize ;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% stats

rng(123)

nBoot = 10000 ;
nEnPerms = length(eval_wsbm_K) ;

wsbm_emp_energy = mean(eval_wsbm_K,2) ;
mod_emp_energy =mean(eval_mod_K,2) ;


bootdiff = zeros(nBoot,1);
% boot diff
for idx=1:nBoot

    tmpind = randi(nEnPerms,1,nEnPerms) ;
    bootdiff(idx) = mean(wsbm_emp_energy(tmpind)) - mean(mod_emp_energy(tmpind)) ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% boostrap difference of differences

rng(123)

wsbm_emp_energy = mean(eval_wsbm_K,2) ;
mod_emp_energy =mean(eval_mod_K,2) ;

wsbmRand_emp_energy = mean(eval_wsbmRand_K,2) ;
modRand_emp_energy =mean(eval_modRand_K,2) ;

bootdiffofdiff = zeros(nBoot,1);
% boot diff
for idx=1:nBoot

    tmpind = randi(nEnPerms,1,nEnPerms) ;
    diff1 = mean(wsbmRand_emp_energy(tmpind))- mean(wsbm_emp_energy(tmpind)) ;
    diff2 = mean(modRand_emp_energy(tmpind)) - mean(mod_emp_energy(tmpind)) ;
    
    bootdiffofdiff(idx) = diff2 - diff1 ;
end



