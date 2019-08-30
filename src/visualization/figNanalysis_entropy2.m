
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

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_consensusCAs.mat' ] ;
load(loadName) ;

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_netstatsRes.mat' ] ;
load(loadName) ;

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_motifEntropy.mat' ] ;
load(loadName) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fig stuff

FIGURE_NAME = 'figEntropy_full' ;
outputdir = strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/');
mkdir(outputdir) 

writeit = 1 ;

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

wsbm_color = [         0    0.4470    0.7410] ;
mod_color = [0.8500    0.3250    0.0980] ;

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

set(gca,'ytick',(1:nComm))
set(gca,'yticklabel',{1:nComm})
set(gca,'ticklength',[ 0 0]) 

set(gca,'xtick',(1:nComm))
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
%% lets bootstrap

rng(123)

nBoot = 10000 ;
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

%% colorbar

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

rng(123)

scatAlpha = 0.5

set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.42, 0.8]);

[wsbm_r2,~,~,wsbm_coefs] = nc_FitAndEvaluateModels( entr_K.wsbm.sum{ind}, nodeWei', ...
                    'linear', 1, 5000) ;
[mod_r2,~,~,mod_coefs] = nc_FitAndEvaluateModels( entr_K.mod.sum{ind}, nodeWei', ...
                    'linear', 1, 5000) ;

% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right])
sp = tight_subplot(2,1,[ .035 .12 ],[.12 .12],[.12 .12]);

axes(sp(1)) 

s = scatter(nodeWei,entr_K.wsbm.sum{ind},'filled')
s.MarkerFaceAlpha = scatAlpha ;

axis square

viz_FnE_Regression(wsbm_coefs)

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

viz_FnE_Regression(mod_coefs)

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

%% edge weight vs edge entropy

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

%% lets bootstrap again

rng(123)

dat_nonan = log(dat(~nanE)) ;

nBoot = 10000 ;
boot1 = zeros(nBoot,1) ;
boot2 = zeros(nBoot,1) ;
nEdges = sum(~nanE(:)) ;

tmp1 = entr_K.wsbm.ent{ind}(~nanE) ;
tmp2 = entr_K.mod.ent{ind}(~nanE) ;

for idx = 1:nBoot
    disp(idx)
    bootsamp = randi(nEdges,1,nEdges) ;
    boot1(idx) = corr(dat_nonan(bootsamp),tmp1(bootsamp),'Type','Spearman') ;
    boot2(idx) = corr(dat_nonan(bootsamp),tmp2(bootsamp),'Type','Spearman') ;
end

ci_1 = prctile(boot1,[ 2.5 97.5]) ;
ci_2 = prctile(boot2,[ 2.5 97.5]) ;

bootdiff = boot1 - boot2 ;

pvalue = mean(bootdiff<0);
pvalue = 2*min(pvalue,1-pvalue);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% scratch
% 
% lll = load('data/raw/CCTX_matrix.mat') ;
% 
% [se,sortent] = sort(entr_K.wsbm.sum{ind},'descend') ;
% [ss,sortstr] = sort(nodeWei,'descend') ;
% 
% [s,sortentm] = sort(entr_K.mod.sum{ind},'descend') ;
% 
% lll.area_names(sortentm)'
% 
% lll.area_names(sortstr)'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% lets bootstrap between null and observed

[~,~,wwwd] = strengths_dir(dat);
[~,~,bbbd] = degrees_dir(dat) ;

rng(123)

nBoot = 10000 ;
boot1 = zeros(nBoot,1) ;
boot2 = zeros(nBoot,1) ;
nNodes = length(wwwd) ;

for idx = 1:nBoot
    %disp(idx)
    bootsamp = randi(nNodes,1,nNodes) ;
    tmp1 = entr_K.wsbm.sum{ind}(bootsamp) ;
    tmp2 = entr_K.randblW.sum{ind}(bootsamp) ;
    boot1(idx) = corr(bbbd(bootsamp)',tmp1,'Type','Spearman') ;
    boot2(idx) = corr(bbbd(bootsamp)',tmp2,'Type','Spearman') ;
end

ci_1 = prctile(boot1,[ 2.5 97.5]) ;
ci_2 = prctile(boot2,[ 2.5 97.5]) ;

bootdiff = boot2 - boot1;

pvalue = mean(bootdiff<0);
pvalue = 2*min(pvalue,1-pvalue);


%% trvial

rng(123)

nBoot = 10000 ;
boot1 = zeros(nBoot,1) ;
boot2 = zeros(nBoot,1) ;
nNodes = length(wwwd) ;

for idx = 1:nBoot
    %disp(idx)
    bootsamp = randi(nNodes,1,nNodes) ;
    tmp1 = entr_K.mod.sum{ind}(bootsamp) ;
    tmp2 = entr_K.randblM.sum{ind}(bootsamp) ;
    boot1(idx) = corr(bbbd(bootsamp)',tmp1,'Type','Spearman') ;
    boot2(idx) = corr(bbbd(bootsamp)',tmp2,'Type','Spearman') ;
end

ci_1 = prctile(boot1,[ 2.5 97.5]) ;
ci_2 = prctile(boot2,[ 2.5 97.5]) ;

bootdiff = boot2 - boot1;

pvalue = mean(bootdiff<0);
pvalue = 2*min(pvalue,1-pvalue);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% show that the percent changes are not due to degree

rng(123)

nBoot = 10000 ;
boot1 = zeros(nBoot,1) ;
boot2 = zeros(nBoot,1) ;
nNodes = length(wwwd) ;

prnctDiffWSBM = (entr_K.randblW.sum{ind} - entr_K.wsbm.sum{ind}) ./ entr_K.randblW.sum{ind} ;
prnctDiffMOD = (entr_K.randblM.sum{ind} - entr_K.mod.sum{ind}) ./ entr_K.randblM.sum{ind} ;

for idx = 1:nBoot
    %disp(idx)
    bootsamp = randi(nNodes,1,nNodes) ;
    boot1(idx) = corr(bbbd(bootsamp)',prnctDiffWSBM(bootsamp),'Type','Spearman') ;
    boot2(idx) = corr(bbbd(bootsamp)',prnctDiffMOD(bootsamp),'Type','Spearman') ;
end

ci_1 = prctile(boot1,[ 2.5 97.5]) ;
ci_2 = prctile(boot2,[ 2.5 97.5]) ;

% both include 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot 

rng(123)

scatAlpha = 0.5

[~,~,bdegree] = degrees_dir(dat) ;

set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.42, 0.5]);

tmp = (entr_K.wsbm.sum{ind} - entr_K.randblW.sum{ind}) ./ entr_K.randblW.sum{ind} ;

s = scatter(bdegree,tmp,'filled')
s.MarkerFaceAlpha = scatAlpha ;
s.MarkerFaceColor = wsbm_color ;

axis square

yl = ylabel('Change node entropy')
yl.FontSize = fontsize ;

%mod

tmp = (entr_K.mod.sum{ind} - entr_K.randblM.sum{ind}) ./ entr_K.randblM.sum{ind} ;

hold on
s = scatter(bdegree,tmp,'filled')
s.MarkerFaceAlpha = scatAlpha ;
s.MarkerFaceColor = mod_color ;

axis square

a=[cellstr(num2str(get(gca,'ytick')'*100))]; 
pct = char(ones(size(a,1),1)*'%'); 
new_yticks = [char(a),pct];
set(gca,'yticklabel',new_yticks)

xl = xlabel('Node degree')
xl.FontSize = fontsize ;

ll = legend({'WSBM' 'Modular'})
ll.FontSize = fontsize ;

if writeit
    fileName = strcat('node_ent_scatter_chgnprcnt.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot the null with the actual

rng(123)

scatAlpha1 = 0.5
scatAlpha2 = 0.2

[~,~,bdegree] = degrees_dir(dat) ;

set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.42, 0.8]);

% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right])
sp = tight_subplot(2,1,[ .035 .12 ],[.12 .12],[.12 .12]);

axes(sp(1)) 

s = scatter(bdegree,entr_K.wsbm.sum{ind},'filled')
s.MarkerFaceAlpha = scatAlpha1 ;
%s.MarkerEdgeColor = wsbm_color; 
s.MarkerFaceColor = wsbm_color ;
hold on
s = scatter(bdegree,entr_K.randblW.sum{ind},'filled')
s.MarkerFaceAlpha = scatAlpha2 ;
%s.MarkerEdgeColor = wsbm_color; 
s.MarkerFaceColor = wsbm_color ;

axis square

%xlim([0 32])
xl = xlim

yl = ylabel('WSBM node entropy')
yl.FontSize = fontsize ;
ypos = yl.Position ;
yl.Position = [ ypos(1)*1.1 ypos(2:3) ] ;
ypos = yl.Position ;

% xl = xlabel('Node strength')
% xl.FontSize = fontsize ;

xticks([])

%mod

axes(sp(2))

s = scatter(bdegree,entr_K.mod.sum{ind},'filled')
s.MarkerFaceAlpha = scatAlpha1 ;
%s.MarkerEdgeColor = mod_color; 
s.MarkerFaceColor = mod_color ;
hold on
s = scatter(bdegree,entr_K.randblM.sum{ind},'filled')
s.MarkerFaceAlpha = scatAlpha2 ;
%s.MarkerEdgeColor = mod_color; 
s.MarkerFaceColor = mod_color ;
axis square

xlim(xl)

yl = ylabel('Modular node entropy')
yl.FontSize = fontsize ;
ypos_tmp = yl.Position 
yl.Position = [ ypos(1) ypos_tmp(2:3) ] ;

xl = xlabel('Node degree')
xl.FontSize = fontsize ;

if writeit
    fileName = strcat('node_ent_scatter_vert.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

aaa = assortativity_bin(dat) ;




