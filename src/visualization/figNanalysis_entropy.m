
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

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_consensusCAs.mat' ] ;
load(loadName) ;

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_netstatsRes.mat' ] ;
load(loadName) ;

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_motifEntropy.mat' ] ;
load(loadName) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fig stuff

FIGURE_NAME = 'figEntropy' ;
outputdir = strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/');
mkdir(outputdir) 

writeit = 0 ;

fontsize = 16 ;

% general graphics, this will apply to any figure you open
% (groot is the default figure object).
set(groot, ...
'DefaultFigureColor', 'w', ...
'DefaultAxesFontSize', 14, ...
'DefaultAxesFontName', 'Arial', ...
'DefaultLineLineWidth', 1, ...
'DefaultTextFontUnits', 'Points', ...
'DefaultTextFontSize', 16, ...
'DefaultTextFontName', 'Arial', ...
'DefaultAxesBox', 'off', ...
'DefaultAxesTickLength', [0.02 0.025]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setup some vars

nanE = isnan(baseRes.rawData) ;

dat = baseRes.rawData ;
dat(isnan(dat)) = 0 ;

binDat = dat > 0 ;

ind = baseRes.wsbm.bestKind ;

nodeWei = netstats.node.str ;
binNodeWei = strengths_dir(binDat) ; 

% ent_cm = brewermap(150,'RdPu') ;
% %ent_cm = brewermap(150,'YlOrRd') ;
% ent_cm = ent_cm(1:100,:) ;

ent_cm = brewermap(150,'Purples') ;
%ent_cm = brewermap(150,'YlOrRd') ;
ent_cm = ent_cm(1:100,:) ;

% comms_color = [ 0.5 0.5 0.5 ] ;
comms_color = [ .3 .3 .3 ] ;
comms_thick = 1.75 ;

nComm = baseRes.wsbm.bestK ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot entropy images

comTypes = { 'wsbm' 'mod' } ;

for cmTyp = 1:length(comTypes)

figure

[h,sss] = imsc_grid_comm(entr_K.(comTypes{cmTyp}).ent{ind}, ...
    cons_ca.(comTypes{cmTyp}),...
    comms_thick,comms_color) ;
set(h,'alphadata',nanE(sss,sss) == 0);

colormap(ent_cm)
colorbar

axis square

set(gca,'FontSize',fontsize)
set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.42, 0.5]);

xticks([])
yticks([]) 

%viz_comms_on_axes(cons_ca.(comTypes{cmTyp}), comm_cmap)
viz_labs_on_axes(cons_ca.(comTypes{cmTyp}))

%tightfig

if writeit
    fileName = strcat(comTypes{cmTyp}, '_edge_entropies.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% and the block

for cmTyp = 1:length(comTypes)

figure
    
[~,bl] = get_block_mat(entr_K.(comTypes{cmTyp}).ent{ind},cons_ca.(comTypes{cmTyp})) ;

%Plot the Matrix
h = imagesc(bl);
set(h,'alphadata',bl ~= 0);
colormap(ent_cm)
%h.CDataMapping = 'direct';

lo = min(entr_K.(comTypes{cmTyp}).ent{ind}(:)) ;
hi = max(entr_K.(comTypes{cmTyp}).ent{ind}(:)) ;

caxis([lo hi])
colorbar()

axis square

xticks([])
yticks([]) 

set(gca,'ytick',(1:10))
set(gca,'yticklabel',{1:nComm})
set(gca,'ticklength',[ 0 0]) 

set(gca,'xtick',(1:10))
set(gca,'xticklabel',{1:nComm})
set(gca,'ticklength',[ 0 0]) 

set(gca,'FontSize',fontsize)
set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.42, 0.5]);

if writeit
    fileName = strcat(comTypes{cmTyp}, '_blockedge_entropies.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% node weight vs node entropy

scatAlpha = 0.5 ;

rng(123)

% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right])
sp = tight_subplot(1,2,[ .12 .12 ],[.12 .12],[.12 .12]);

axes(sp(1))
s = scatter(nodeWei,entr_K.wsbm.sum{ind},'filled')
s.MarkerFaceAlpha = scatAlpha ;

axis square

%quick stats
% function [xvalR2, xvalsqErr, yhatLOOCV, coefStruct , lsFitStruct , permStruct ] = ... 
%    nc_FitAndEvaluateModels(y, x, model, crossvalidate, bootIter, params , permIter)
[wsbm_r2,~,~,wsbm_coefs_nodes] = nc_FitAndEvaluateModels( entr_K.wsbm.sum{ind}, nodeWei', ...
                    'linear', 1, 5000) ;
viz_FnE_Regression(wsbm_coefs_nodes)

yl = ylabel('WSBM node entropy')
yl.FontSize = fontsize ;

xl = xlabel('Node strength')
xl.FontSize = fontsize ;


%set(gca,'FontSize',fontsize)

%%%%%%%%%%%%%

axes(sp(2)) 
s= scatter(nodeWei,entr_K.mod.sum{ind},'filled')
s.MarkerFaceAlpha = scatAlpha ;

axis square

[mod_r2,~,~,mod_coefs_nodes] = nc_FitAndEvaluateModels( entr_K.mod.sum{ind}, nodeWei', ...
                    'linear', 1, 5000) ;
viz_FnE_Regression(mod_coefs_nodes)

yl = ylabel('Modular node entropy')
yl.FontSize = fontsize ;

xl = xlabel('Node strength')
xl.FontSize = fontsize ;

%set(gca,'FontSize',fontsize)

set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.42, 0.5]);

if writeit
    fileName = strcat('node_ent_scatter.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% lets bootstrap

rng(123)

nBoot = 5000 ;
boot1 = zeros(nBoot,1) ;
boot2 = zeros(nBoot,1) ;
nNodes = length(netstats.node.str) ;

for idx = 1:nBoot
    disp(idx)
    bootsamp = randi(nNodes,1,nNodes) ;
    boot1(idx) = corr(nodeWei(bootsamp)',entr_K.wsbm.sum{ind}(bootsamp),'Type','Spearman') ;
    boot2(idx) = corr(nodeWei(bootsamp)',entr_K.mod.sum{ind}(bootsamp),'Type','Spearman') ;
end

ci_1 = prctile(boot1,[ 2.5 97.5]) ;
ci_2 = prctile(boot2,[ 2.5 97.5]) ;

bootdiff = boot1 - boot2 ;

pvalue = mean(bootdiff<0);
pvalue = 2*min(pvalue,1-pvalue);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% edgewise 

scatAlpha = 0.05 ;

rng(123)

% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right])
sp = tight_subplot(1,2,[ .12 .12 ],[.12 .12],[.12 .12]);

axes(sp(1))
s = scatter(log(dat(~nanE)),entr_K.wsbm.ent{ind}(~nanE),50,'filled') ;
%colormap([ 0.4 0.4 0.4 ; 0.8 0.8 0.8] ) ;
s.MarkerFaceAlpha = scatAlpha ;
% ss = scatplot(dat(~nanE),entr_K.wsbm.ent{ind}(~nanE))

axis square

%quick stats
% function [xvalR2, xvalsqErr, yhatLOOCV, coefStruct , lsFitStruct , permStruct ] = ... 
%    nc_FitAndEvaluateModels(y, x, model, crossvalidate, bootIter, params , permIter)
[wsbm_r2,~,~,wsbm_coefs] = nc_FitAndEvaluateModels( entr_K.wsbm.ent{ind}(~nanE), log(dat(~nanE)), ...
                    'linear', 1, 5000) ;
viz_FnE_Regression(wsbm_coefs)

yl = ylabel('WSBM edge entropy')
yl.FontSize = fontsize ;

xl = xlabel('Log edge weight')
xl.FontSize = fontsize ;


%set(gca,'FontSize',fontsize)

%%%%%%%%%%%%%

axes(sp(2)) 
s = scatter(log(dat(~nanE)),entr_K.mod.ent{ind}(~nanE),50,'filled') ;
%colormap([ 0.4 0.4 0.4 ; 0.8 0.8 0.8] ) ;
s.MarkerFaceAlpha = scatAlpha ;
axis square

[mod_r2,~,~,mod_coefs] = nc_FitAndEvaluateModels( entr_K.mod.ent{ind}(~nanE), log(dat(~nanE)), ...
                    'linear', 1, 5000) ;
viz_FnE_Regression(mod_coefs)

yl = ylabel('Modular edge entropy')
yl.FontSize = fontsize ;

xl = xlabel('Log edge weight')
xl.FontSize = fontsize ;

%set(gca,'FontSize',fontsize)

set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.42, 0.5]);

if writeit
    fileName = strcat('edge_ent_scatter.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
end

%% WITH WEIGHTING

numBins = 10 ;

uniq_vals = unique(log(dat(~nanE))) ;
dat_nonan = log(dat(~nanE)) ;

den_cmap = brewermap(numBins,'RdYlBu') ;
datapointSz = 50 ;

nnn = { 'WSBM' 'Modular' } ;

% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right])
sp = tight_subplot(1,2,[ .12 .12 ],[.12 .12],[.12 .12]);

for cmTyp = 1:length(comTypes)

%figure
    
dattt = entr_K.(comTypes{cmTyp}).ent{ind}(~nanE) ;

axes(sp(cmTyp))

s1 = scatter(dat_nonan, dattt) ;
s1.MarkerEdgeAlpha = 0 ;

for idx = 1:length(uniq_vals)
    
    currval = uniq_vals(idx) ;
    currDat = dattt(dat_nonan==currval) ;
    % histogram it
    quantE = quantile(currDat,numBins-1) ;
%     [binCount,hc] = histcounts(currDat,numBins) ;
    [~,hc] = histcounts(currDat,[ min(currDat) quantE max(currDat) ]) ;
    discVals = discretize(currDat,hc) ;
    
    currValVec = ones(length(currDat),1) .* currval ;
    
    hold on
    s = scatter(currValVec,currDat,datapointSz,discVals,'filled') ;
    s.MarkerFaceAlpha = .1 ;
    
end

colormap(den_cmap)

yl = ylabel([ nnn{cmTyp} ' edge entropy'])
yl.FontSize = fontsize ;

xl = xlabel('Log edge weight')
xl.FontSize = fontsize ;

axis square

hold off

set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.42, 0.5]);

end

if writeit
    fileName = strcat(comTypes{cmTyp}, '_edge_ent_wColor.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
end

%% lets bootstrap again

rng(123)

dat_nonan = log(dat(~nanE)) ;

nBoot = 5000 ;
boot1 = zeros(nBoot,1) ;
boot2 = zeros(nBoot,1) ;
nEdges = sum(~nanE(:)) ;

for idx = 1:nBoot
    disp(idx)
    bootsamp = randi(nEdges,1,nEdges) ;
    tmp1 = entr_K.wsbm.ent{ind}(~nanE) ;
    tmp2 = entr_K.mod.ent{ind}(~nanE) ;
    boot1(idx) = corr(dat_nonan(bootsamp),tmp1(bootsamp),'Type','Spearman') ;
    boot2(idx) = corr(dat_nonan(bootsamp),tmp2(bootsamp),'Type','Spearman') ;
end

ci_1 = prctile(boot1,[ 2.5 97.5]) ;
ci_2 = prctile(boot2,[ 2.5 97.5]) ;

bootdiff = boot1 - boot2 ;

pvalue = mean(bootdiff<0);
pvalue = 2*min(pvalue,1-pvalue);

%%

numBins = 10
den_cmap = brewermap(numBins,'RdYlGn') ;

imagesc(1:numBins)
colormap(den_cmap)
colorbar()

% get percentages
ppp = ((1:numBins) ./ numBins) .* 100 ;
ppp = ppp(:) ;

labn = cell(numBins,1) ;
for idx = 1:numBins
    labn{idx} = sprintf('%i%%',round(ppp(idx))) ;
end

cmap_labs_discrete(labn)

set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.6, 0.5]);

if writeit
    fileName = strcat( 'colorbar.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% vertical figs

set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.42, 0.8]);


% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right])
sp = tight_subplot(2,1,[ .035 .12 ],[.12 .12],[.12 .12]);

axes(sp(1)) 

s = scatter(nodeWei,entr_K.wsbm.sum{ind},'filled')
s.MarkerFaceAlpha = scatAlpha ;

axis square

viz_FnE_Regression(wsbm_coefs_nodes)

xlim([0 32])

yl = ylabel('WSBM node entropy')
yl.FontSize = fontsize ;
ypos = yl.Position ;

% xl = xlabel('Node strength')
% xl.FontSize = fontsize ;

xticks([])

%mod

axes(sp(2))

s = scatter(nodeWei,entr_K.mod.sum{ind},'filled')
s.MarkerFaceAlpha = scatAlpha ;

axis square

viz_FnE_Regression(mod_coefs_nodes)

xlim([0 32])

yl = ylabel('Modular node entropy')
yl.FontSize = fontsize ;
ypos_tmp = yl.Position 
yl.Position = [ ypos(1) ypos_tmp(2:3) ] ;

xl = xlabel('Node strength')
xl.FontSize = fontsize ;

if writeit
    fileName = strcat('node_ent_scatter_vert.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
end

%%

numBins = 10 ;

uniq_vals = unique(log(dat(~nanE))) ;
dat_nonan = log(dat(~nanE)) ;

den_cmap = brewermap(numBins,'RdYlGn') ;
datapointSz = 50 ;

nnn = { 'WSBM' 'Modular' } ;

set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.42, 0.8]);

% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right])
sp = tight_subplot(2,1,[ .035 .12 ],[.12 .12],[.12 .12]);

for cmTyp = 1:length(comTypes)

%figure
    
dattt = entr_K.(comTypes{cmTyp}).ent{ind}(~nanE) ;

axes(sp(cmTyp))

s1 = scatter(dat_nonan, dattt) ;
s1.MarkerEdgeAlpha = 0 ;

for idx = 1:length(uniq_vals)
    
    currval = uniq_vals(idx) ;
    currDat = dattt(dat_nonan==currval) ;
    % histogram it
    quantE = quantile(currDat,numBins-1) ;
%     [binCount,hc] = histcounts(currDat,numBins) ;
    [~,hc] = histcounts(currDat,[ min(currDat) quantE max(currDat) ]) ;
    discVals = discretize(currDat,hc) ;
    
    currValVec = ones(length(currDat),1) .* currval ;
    
    hold on
    s = scatter(currValVec,currDat,datapointSz,discVals,'filled') ;
    s.MarkerFaceAlpha = .1 ;
        
end

colormap(den_cmap)

yl = ylabel([ nnn{cmTyp} ' edge entropy'])
yl.FontSize = fontsize ;

if cmTyp == 2
xl = xlabel('Log edge weight')
xl.FontSize = fontsize ;
else
    xticks([])
end

axis square

hold off

%set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.42, 0.5]);

end

if writeit
    fileName = strcat(comTypes{cmTyp}, '_edge_ent_wColor_vert.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
end

