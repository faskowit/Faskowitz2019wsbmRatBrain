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

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_consensusCAs.mat' ] ;
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
caDistances = mean(partition_distance(ca_at_Kbest),2) ;

logEvid_at_Kbest = baseRes.wsbm.logEvid_K{baseRes.wsbm.bestKind} ;

medianLogEvid = median(logEvid_at_Kbest) ;
medianDist = median(caDistances) ;

% get indicies for top quadrant
optQuadInd = (logEvid_at_Kbest > medianLogEvid) & (caDistances < medianDist) ;  

% viz it
figure
s = scatter(logEvid_at_Kbest(~optQuadInd),caDistances(~optQuadInd),...
    [],caDistances(~optQuadInd),'filled') ;
%colormap([ 0.4 0.4 0.4 ; 0.8 0.8 0.8] ) ;
s.MarkerFaceAlpha = .25 ;
hold on
s = scatter(logEvid_at_Kbest(optQuadInd),caDistances(optQuadInd),...
    [],caDistances(optQuadInd),'filled') ;
%colormap([ 0.4 0.4 0.4 ; 0.8 0.8 0.8] ) ;
s.MarkerFaceAlpha = .95 ;

yl = ylabel('VI Distance') ;
xl = xlabel('Model Log Evidence') ;
yl.FontSize = fontsize ;
xl.FontSize = fontsize ;

axis square

set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.4, 0.5]);
tightfig

%qucik stats
% function [xvalR2, xvalsqErr, yhatLOOCV, coefStruct , lsFitStruct , permStruct ] = ... 
%    nc_FitAndEvaluateModels(y, x, model, crossvalidate, bootIter, params , permIter)

[r2,~,~,coefs] = nc_FitAndEvaluateModels(caDistances, logEvid_at_Kbest, ...
                    'linear', 1, 5000) ;

% get the confidence intervals
xVec = min(coefs.x):0.05:...
    max(coefs.x);
bootCI = zeros(size(coefs.boot,1),length(xVec)) ;
for idx = 1:size(coefs.boot,1)
    bootCI(idx,:) = polyval(coefs.boot(idx,:),xVec);  
end
p95 = prctile(bootCI,[2.5,97.5]);    
hold on;
fill([xVec fliplr(xVec)],[p95(1,:) fliplr(p95(2,:))],[ 0.75 0.75 0.75 ],'facealpha',.7,'edgealpha',0);

% and the line
yhat = polyval(coefs.full,xVec); 
plot(xVec,yhat,'Color',[0.5 0.5 0.5],'linewidth',2.5);


if writeit 
    fileName = strcat('topModels.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
end

%% plot what we got for modularity

dat = baseRes.rawData ;
dat(isnan(dat)) = 0 ;

rng(123)

figure
set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.3, 0.5]);

% !!! instead, use the Q that comes out of gen louvain !!! 
% tmpMods = cell(length(COM_NUM_RANGE),1) ;
% for idx = 1:length(COM_NUM_RANGE)
% 
%     tmpCAs = baseRes.mod.ca_K{idx} ;
%     tmpMods{idx} = zeros(size(tmpCAs,2),1) ;
%     for jdx = 1:size(tmpCAs,2)
%        tmpMods{idx}(jdx) = modularity_q(dat,tmpCAs(:,jdx)) ;
%     end
% end

for idx = 1:length(COM_NUM_RANGE)

   tmpData = baseRes.mod.caQs_K{idx} ;
   
   tmpDists = mean(partition_distance(baseRes.mod.ca_K{idx}),2); 
   
   scatter(COM_NUM_RANGE(idx) .* ones(length(tmpData),1) + (0.4).*rand(length(tmpData),1) -0.2 ,...
       tmpData,[],tmpDists,'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.7)
   hold on
   
end

cb = colorbar ;
cb.Label.String = 'VI distance' ;
cb.Label.FontSize = 16 ;

xl = xlabel('Number of communities ({\itk})') ;
xl.FontSize = fontsize ;
yl = ylabel('Modularity') ;
yl.FontSize = fontsize ;

%set(gca,'FontSize',16)
%ylim = get(gca,'Ylim') ;
%plot([baseRes.wsbm.bestK baseRes.wsbm.bestK],ylim,'r','LineStyle',':','LineWidth',2)

if writeit 
    fileName = strcat('modularity_K.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot the communities

figure 

nComm = baseRes.wsbm.bestK ;
comm_cmap = brewermap(nComm,'Spectral') ;

wsbmCaBest = zeros(size(baseRes.wsbm.ca_K{baseRes.wsbm.bestKind})) ;
wsbmCaBestLogE = baseRes.wsbm.logEvid_K{baseRes.wsbm.bestKind} ;
modCaBest = zeros(size(baseRes.mod.ca_K{baseRes.wsbm.bestKind})) ;
modCaBestLogE = baseRes.mod.caQs_K{baseRes.wsbm.bestKind} ;

nRep = size(wsbmCaBest,2) ;

for idx = 1:nRep 
    wsbmCaBest(:,idx) = hungarianMatch(cons_ca.wsbm,...
        baseRes.wsbm.ca_K{baseRes.wsbm.bestKind}(:,idx)) ;
end
for idx = 1:nRep 
    modCaBest(:,idx) = hungarianMatch(cons_ca.mod,...
        baseRes.mod.ca_K{baseRes.wsbm.bestKind}(:,idx)) ;
end

% tight_subplot(Nh, Nw, [gap_h gap_w], [lower upper], [left right])
sp = tight_subplot(2,1,[ .035 .04 ],[.12 .033],[.1 .1]);

axes(sp(1))
imagesc(wsbmCaBest(:,sortedInd(wsbmCaBestLogE))) ;
set(gca,'XTickLabel',[])
cmap_labs_discrete(1:nComm)

axes(sp(2))
imagesc(modCaBest(:,sortedInd(modCaBestLogE))) ;

colormap(comm_cmap) ;
cmap_labs_discrete(1:nComm)

if writeit 
    fileName = strcat('aligned_comms.png');
    ff = fullfile(strcat(outputdir,'/',OUTPUT_STR,'_',fileName)); 
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
end





