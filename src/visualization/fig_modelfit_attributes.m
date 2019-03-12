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

%% plot what we got

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

xlabel('Number of communities ({\itk})')
ylabel('Model log evidence')

set(gca,'FontSize',16)

ylim = get(gca,'Ylim') ;

plot([baseRes.wsbm.bestK baseRes.wsbm.bestK],ylim,'r','LineStyle',':','LineWidth',2)

%%

if writeit 
    fileName = strcat('logevidence_K.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
%     fileName = strcat('logevidence_K_2.png');
%     ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
%     export_fig(ff,'-png')
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
ylabel('VI Distance') 
xlabel('Model Log Evidence')

set(gca,'FontSize',16)
set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.5, 0.5]);
pbaspect([1 1 1])

if writeit 
    fileName = strcat('topModels.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
end

