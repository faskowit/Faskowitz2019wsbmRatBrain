%% clear stuff

clc
clearvars

%% load the necessary data

config_file='config_template_rb2_oneHemi_analyzeGridRuns.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

%% load data

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_baseRes.mat' ] ;
% loads a struct named 'baseRes'
load(loadName) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fig stuff

FIGURE_NAME = 'figModelfitting' ;
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


%% plot what we got

rng(123)

figure
set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.3, 0.5]);

for idx = 1:length(COM_NUM_RANGE)

   tmpData = baseRes.wsbm.logEvid_K{idx} ;%cellfun(getLogEvid,kResults{idx}) ;
   
   tmpDists = mean(partition_distance(baseRes.wsbm.ca_K{idx}),2); 
   
   scatter(COM_NUM_RANGE(idx) .* ones(length(tmpData),1) + (0.4).*rand(length(tmpData),1) -0.2 ,...
       tmpData,[],tmpDists,'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.7)
   hold on
   
end

cb = colorbar ;
cb.Label.String = 'VI distance' ;
cb.Label.FontSize = 16 ;

xl = xlabel('Number of communities ({\itk})') ;
xl.FontSize = fontsize ;
yl = ylabel('Model log evidence') ;
yl.FontSize = fontsize ;

%set(gca,'FontSize',16)

ylim = get(gca,'Ylim') ;

plot([baseRes.wsbm.bestK baseRes.wsbm.bestK],ylim,'r','LineStyle',':','LineWidth',2)

if writeit 
    fileName = strcat('logevidence_K.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

ca_at_Kbest = baseRes.wsbm.ca_K{baseRes.wsbm.bestKind} ;
caDistances = sum(partition_distance(ca_at_Kbest),2) ;

logEvid_at_Kbest = baseRes.wsbm.logEvid_K{baseRes.wsbm.bestKind} ;

medianLogEvid = median(logEvid_at_Kbest) ;
medianDist = median(caDistances) ;

% get indicies for top quadrant
optQuadInd = (logEvid_at_Kbest > medianLogEvid) & (caDistances < medianDist) ;  

% viz it
figure
scatter(logEvid_at_Kbest,caDistances,[],optQuadInd)
colormap(brewermap(2,'Paired')) ;
yl = ylabel('VI Distance') ;
xl = xlabel('Model Log Evidence') ;
yl.FontSize = fontsize ;
xl.FontSize = fontsize ;

axis square

set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.38, 0.5]);
tightfig

if writeit 
    fileName = strcat('topModels.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
end

