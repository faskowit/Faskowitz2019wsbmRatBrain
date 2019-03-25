
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

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_motifEntropy.mat' ] ;
load(loadName)

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_consensusCAs.mat' ] ;
load(loadName) ;

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_motifAna.mat' ] ;
load(loadName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fig stuff

FIGURE_NAME = 'figEdgeMotifProbs' ;
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
%% viz at particular scale

nanEdges = isnan(baseRes.rawData) ;

[~,wsbm_sorti] = sort(cons_ca.wsbm) ;
[~,mod_sorti] = sort(cons_ca.mod) ;

curScale = baseRes.wsbm.bestKind ;

edgeMats = { 'odMat' 'aMat' 'cMat' 'pMat' 'dMat' } ;
egdeMatNames = {'on-diag' 'assort' 'core' 'periph' 'disassort'} ;

mats_cm =  brewermap(5,'Paired') ;
matscm_width = 9 ;

probs_cm = brewermap(100,'OrRd') ;
probs_cm(1,:) = [ 1 1 1 ] ;

comms_color = [ 0.5 0.5 0.5 ] ;
comms_thick = 0.75 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right])
sp = tight_subplot(2,5,[ .02 .02 ],[.05 .1],[.02 .02]);

for idx = 1:length(edgeMats)

    tmp = mean(motifAna.wsbm.motifEdgeMats{curScale}.(edgeMats{idx}) > 0,3) ;
    %subplot(2,5,idx) ; 
    
    axes(sp(idx)) 
   
%     h = imagesc(tmp(wsbm_sorti,wsbm_sorti)) ;
    imsc_grid_comm(tmp,cons_ca.wsbm,comms_thick,comms_color)
    caxis([0 1])
    colormap(probs_cm)
    xticks([])
    yticks([]) 
    %set(h,'alphadata',nanEdges==0);
    viz_imsc_border(mats_cm(idx,:),matscm_width)
    axis square
    
    tt = title(egdeMatNames{idx}) ;
    tt.FontSize = fontsize ;
    tt.FontWeight = 'normal' ;
    pos = tt.Position ;
    tt.Position = [ pos(1) pos(2)*1.5 pos(3) ] ;
    
end

for idx = 1:length(edgeMats)

    tmp = mean(motifAna.mod.motifEdgeMats{curScale}.(edgeMats{idx}) > 0,3) ;
    %subplot(2,5,idx+5) ; 
    
    axes(sp(idx+5)) 
    
%     h = imagesc(tmp(mod_sorti,mod_sorti)) ;
    imsc_grid_comm(tmp,cons_ca.mod,comms_thick,comms_color)
    caxis([0 1])
    colormap(probs_cm)
    xticks([])
    yticks([]) 
    %set(h,'alphadata',nanEdges==0);
    viz_imsc_border(mats_cm(idx,:),matscm_width)
    axis square
end

set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.6, 0.5]);

if writeit
    fileName = strcat( 'edgeprobs_at_bestK.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get colormat

figure

imagesc(0:0.1:1)
colormap(probs_cm)
colorbar()

set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.6, 0.5]);

if writeit
    fileName = strcat( 'colorbar.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
end
