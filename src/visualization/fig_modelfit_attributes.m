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


