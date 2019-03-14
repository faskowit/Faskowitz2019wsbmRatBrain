
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

writeit = 1 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% image with the comms

figure

nComm = length(unique(cons_ca.wsbm));

comm_cmap = brewermap(nComm,'Spectral') ;

[~,breaks] = wsbm_plot_mat(baseRes.rawData,dummyvar(cons_ca.wsbm)','') ;
hold on
viz_comms_on_axes(cons_ca.wsbm,comm_cmap)
       
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

set(gca,'FontSize',15)
set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.5, 0.5]);
pbaspect([1 1 1])

if writeit
    fileName = strcat('wsbm_comms.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% wsbm block matrix

[~,avgBl] = get_block_mat(baseRes.rawData, cons_ca.wsbm) ;

%Plot the Matrix
h = imagesc(avgBl,[min(avgBl(:))-.00001,max(avgBl(:))]);
set(h,'alphadata',avgBl~=0);
colorbar();
axis square

hold on
viz_comms_on_axes((1:10)',comm_cmap)

xticks([])
yticks([]) 

set(gca,'ytick',(1:10))
set(gca,'yticklabel',{1:nComm})
set(gca,'ticklength',[ 0 0]) 

set(gca,'FontSize',15)
set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.5, 0.5]);
pbaspect([1 1 1])

if writeit
    fileName = strcat('wsbm_blockMat_comms.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% and mod

figure

[newinds,breaks] = wsbm_plot_mat(baseRes.rawData,dummyvar(cons_ca.mod)','') ;
hold on

viz_comms_on_axes(cons_ca.mod,comm_cmap)
       
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

set(gca,'FontSize',15)
set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.5, 0.5]);
pbaspect([1 1 1])

if writeit 
    fileName = strcat('mod_comms.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% mod block matrix

[~,avgBl] = get_block_mat(baseRes.rawData, cons_ca.mod) ;

%Plot the Matrix
h = imagesc(avgBl,[min(avgBl(:))-.00001,max(avgBl(:))]);
set(h,'alphadata',avgBl~=0);
colorbar();
axis square

hold on
viz_comms_on_axes((1:10)',comm_cmap)

xticks([])
yticks([]) 

set(gca,'ytick',(1:10))
set(gca,'yticklabel',{1:nComm})
set(gca,'ticklength',[ 0 0]) 

set(gca,'FontSize',15)
set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.5, 0.5]);
pbaspect([1 1 1])

if writeit
    fileName = strcat('mod_blockMat_comms.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
end
