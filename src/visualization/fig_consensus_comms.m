
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

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_consensusCAs.mat' ] ;
load(loadName) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fig stuff

FIGURE_NAME = 'figComms' ;
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
%% image with the comms

comTypes = { 'wsbm' 'mod' } ;

for cmTyp = 1:length(comTypes)
    
figure

thisComm = cons_ca.(comTypes{cmTyp}) ;

nComm = length(unique(thisComm));

comm_cmap = brewermap(nComm,'Spectral') ;

[~,breaks] = wsbm_plot_mat(baseRes.rawData,dummyvar(thisComm)','') ;
hold on
viz_comms_on_axes(cons_ca.(comTypes{cmTyp}),comm_cmap)
       
xticks([])
yticks([]) 

% compute some yticks
breaks2 = [ 0 breaks ] ;
midlabelpoint = zeros([nComm 1]);
for idx = 1:length(breaks)
midlabelpoint(idx) = floor( (breaks2(idx+1) - breaks2(idx)) / 2) + breaks2(idx);  
end

set(gca,'ytick',midlabelpoint)
set(gca,'yticklabel',{1:nComm})
set(gca,'ticklength',[ 0 0]) 

hold off

set(gca,'FontSize',fontsize)
set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.42, 0.5]);
pbaspect([1 1 1])

%tightfig

if writeit
    fileName = strcat(comTypes{cmTyp}, '_comms.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% block matrix

ca_lim = [ 0 0.5 ] ;

for cmTyp = 1:length(comTypes)

figure
    
thisComm = cons_ca.(comTypes{cmTyp}) ;
    
[~,avgBl] = get_block_mat(baseRes.rawData, thisComm) ;

%Plot the Matrix
h = imagesc(avgBl,[min(avgBl(:))-.00001,max(avgBl(:))]);
set(h,'alphadata',avgBl~=0);
colorbar();
axis square

hold on
viz_comms_on_axes((1:10)',comm_cmap)

caxis(ca_lim)

xticks([])
yticks([]) 

set(gca,'ytick',(1:10))
set(gca,'yticklabel',{1:nComm})
set(gca,'ticklength',[ 0 0]) 

set(gca,'FontSize',fontsize)
set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.5, 0.5]);
pbaspect([1 1 1])

%tf = tightfig

if writeit
    fileName = strcat(comTypes{cmTyp} , '_blockMat_comms.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% and just visualize raw data

%Plot the Matrix
h = imagesc(baseRes.rawData);
set(h,'alphadata',~isnan(baseRes.rawData));
axis square

xticks([])
yticks([]) 

set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.5, 0.5]);

tightfig

if writeit
    fileName = strcat('rawData.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
end

