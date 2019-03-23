
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

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_commMotifModes.mat' ] ;
load(loadName)

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_consensusCAs.mat' ] ;
load(loadName) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fig stuff

FIGURE_NAME = 'figModes' ;
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
%% setup

nComm = baseRes.wsbm.bestK ;
nanE = isnan(baseRes.rawData ) ;
modes_cm =  brewermap(5,'Paired') ;

comm_cmap = brewermap(nComm,'Spectral') ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make some plots

comTypes = { 'wsbm' 'mod' } ;

for cmTyp = 1:length(comTypes)

figure

[h,sss] = imsc_grid_comm(eMatStruct.(comTypes{cmTyp}).max,cons_ca.(comTypes{cmTyp})) ;
set(h,'alphadata',nanE(sss,sss) == 0);
colormap(modes_cm)
h.CDataMapping = 'direct';
axis square

set(gca,'FontSize',fontsize)
set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.42, 0.5]);

xticks([])
yticks([]) 

viz_comms_on_axes(cons_ca.(comTypes{cmTyp}), comm_cmap)
viz_labs_on_axes(cons_ca.(comTypes{cmTyp}))

%tightfig

if writeit
    fileName = strcat(comTypes{cmTyp}, '_edge_modes.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
end

end

%%

figure
imagesc(1:5)
h = cmap_labs_discrete({'on-diag' 'assort' 'core' 'periph' 'disassort'}) ; 
set( h, 'YDir', 'reverse' );
colormap(modes_cm)

set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.5, 0.5]);

if writeit
    fileName = strcat( 'colorbar.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
end

%% modes per block

for cmTyp = 1:length(comTypes)

figure
    
bl = eMatStruct.(comTypes{cmTyp}).mode ;

%Plot the Matrix
h = imagesc(bl);
set(h,'alphadata',~isnan(bl));
colormap(modes_cm)
h.CDataMapping = 'direct';
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

tightfig

if writeit
    fileName = strcat(comTypes{cmTyp} , '_block_modes.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
end

end


