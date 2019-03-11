
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

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_consensusRuns.mat' ] ;
load(loadName) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fig stuff

FIGURE_NAME = 'figComms' ;
outputdir = strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/');
mkdir(outputdir) 

writeit = 0 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% image with the comms

conModel = consensusCent ;

nComm = size(conModel.Para.mu,1) ;
ca_wsbm = wsbm_community_assign(conModel) ;

comm_cmap = brewermap(nComm,'Spectral') ;

[newinds,breaks] = wsbm_plot_mat(conModel,[],'mod') ;
hold on
ca_wsbm_ra = lab_reassign(ca_wsbm,newinds) ;

viz_comms_on_axes(ca_wsbm_ra,comm_cmap)
       
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
%% and mod

figure

ca_mod = baseRes.bestKcentmod ;
ca_mod = hungarianMatch(ca_wsbm_ra,ca_mod) ;

[newinds,breaks] = wsbm_plot_mat(conModel.Data.Raw_Data,dummyvar(ca_mod)','') ;
hold on

viz_comms_on_axes(ca_mod,comm_cmap)
       
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

