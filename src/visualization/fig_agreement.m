
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

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_consensusCAs.mat' ] ;
load(loadName) ;

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_motifEntropy.mat' ] ;
load(loadName) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fig stuff

FIGURE_NAME = 'figAgreement_full' ;
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
%% agreement, for supplemental

ind = baseRes.wsbm.bestKind ;
nComm = baseRes.wsbm.bestK ;

% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right])
sp = tight_subplot(2,2,[ .07 .12 ],[.12 .12],[.12 .12]);

set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.42, 0.5]);

axes(sp(1))
imsc_grid_comm(agreement(baseRes.wsbm.ca_K{ind}) ./ 750,cons_ca.wsbm)
axis square
colorbar
caxis([0 1])
title('WSBM')

axes(sp(2))
imsc_grid_comm(agreement(baseRes.mod.ca_K{ind}) ./ 750,cons_ca.mod)
axis square
colorbar
caxis([0 1])
title('Modular')

axes(sp(3)) 
[~,tmp] = get_block_mat(agreement(baseRes.wsbm.ca_K{ind})./750,cons_ca.wsbm) ;
imagesc(tmp)
colorbar
caxis([0 1])
axis square
set(gca,'ytick',(1:10))
set(gca,'yticklabel',{1:nComm})
set(gca,'ticklength',[ 0 0]) 

axes(sp(4)) 
[~,tmp] = get_block_mat(agreement(baseRes.mod.ca_K{ind})./750,cons_ca.mod) ;
imagesc(tmp)
colorbar
caxis([0 1])
axis square
set(gca,'ytick',(1:10))
set(gca,'yticklabel',{1:nComm})
set(gca,'ticklength',[ 0 0]) 

if writeit
    fileName = strcat('agreement_mats.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
end

%%

% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right])
sp = tight_subplot(1,2,[ .07 .12 ],[.12 .12],[.12 .12]);

set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.42, 0.5]);

[~,~,wsbm_nodeAgree] = strengths_dir(agreement(baseRes.wsbm.ca_K{ind})./750) ;
[~,~,mod_nodeAgree] = strengths_dir(agreement(baseRes.mod.ca_K{ind})./750) ;

axes(sp(1)) 
scatter(wsbm_nodeAgree, entr_K.wsbm.sum{ind})
axis square
[a,b] = corr(wsbm_nodeAgree', entr_K.wsbm.sum{ind},'Type','Spearman')
title('WSBM')
yl = ylabel("Node entropy")
yl.FontSize = fontsize ;
xl = xlabel("Node agreement") 
xl.FontSize = fontsize ;

axes(sp(2)) 
scatter(mod_nodeAgree,entr_K.mod.sum{ind})
axis square
[a,b] = corr(mod_nodeAgree',entr_K.mod.sum{ind},'Type','Spearman')
title('Modular')
yl = ylabel("Node entropy")
yl.FontSize = fontsize ;
xl = xlabel("Node agreement") 
xl.FontSize = fontsize ;

if writeit
    fileName = strcat('agreement_corrs.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% agreement between wsbm n mod

cmap = brewermap(3,'PuBu') ;
cmap(1,:) = [ 1 1 1 ] ;

wm_agreement = agreement([cons_ca.wsbm cons_ca.mod]) ;

% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right])
sp = tight_subplot(1,2,[ .07 .12 ],[.12 .12],[.12 .12]);

set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.6, 0.5]);

axes(sp(1))
ims = imsc_grid_comm(wm_agreement+1,cons_ca.wsbm)
ims.CDataMapping = 'Direct'
%set(ims,'alphadata',wm_agreement(sss,sss)~=0);

axis square
cb = colorbar ;
colormap(cmap)
cb.Ticks = [ 1.5 2.5 3.5 ] ;
cb.TickLabels = {'n/a' 'one' 'both' }
%title('WSBM organzied')
xticks([])
yticks([]) 

% [a,b] = unique(sort(cons_ca.wsbm));
% hold on;
% for i=2:length(b)
%     ll = line([0 78],[b(i)-0.5 b(i)-0.5]);
%     uistack(ll,'bottom')
%     ll = line([b(i)-0.5 b(i)-0.5],[0 78]);
%     uistack(ll,'bottom')
% end

viz_labs_on_axes(cons_ca.wsbm)

axes(sp(2))
ims = imsc_grid_comm(wm_agreement+1,cons_ca.mod)
ims.CDataMapping = 'Direct'
axis square
cb = colorbar ;
colormap(cmap)
cb.Ticks = [ 1.5 2.5 3.5 ] ;
cb.TickLabels = {'n/a' 'one' 'both' }
%title('Modular organzied')
xticks([])
yticks([]) 

viz_labs_on_axes(cons_ca.mod)

if writeit
    fileName = strcat('agreement_btwn_models.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

dat = baseRes.rawData ;
dat(isnan(dat)) = 0 ;

w_parti = participation_coef(dat,cons_ca.wsbm) ;
m_parti = participation_coef(dat,cons_ca.mod) ;


