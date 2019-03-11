%% clear stuff

clc
clearvars

%% load the necessary data

config_file='config_template.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

%% load data

loadName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR,'_basicDataSmall_v7p3.mat');
load(loadName)

rawData = load('data/raw/CCTX_matrix.mat') ;

analyzeThis = 1 ;
templateModel = keepRes{analyzeThis}.templateModel ;

loadMod = load(strcat(PROJECT_DIR,'/data/raw/for_Josh.mat')) ;

ca_mod2 = loadMod.ci ;
ca_wsbm = wsbm_community_assign(templateModel);

%% also make mod data with K
imat = templateModel.Data.Raw_Data ;
imat(isnan(imat)) = 0 ;

[ca_mod,ca_modQ] = gsweep_mod_atK_dir(imat,templateModel.R_Struct.k) ;

%% align mu to the wsbm

ca_mod = hungarianMatch(ca_wsbm,ca_mod) ;
% ca_mod2 = hungarianMatch(ca_wsbm,ca_mod2) ;

%% and the modular geneative models

muMod = dummyvar(ca_mod)' ;
[~,modularityModel] = wsbm(templateModel.Data.Raw_Data, ...
    size(muMod,1), ...
    'W_Distr', templateModel.W_Distr, ...
    'E_Distr', templateModel.E_Distr, ...
    'alpha', templateModel.Options.alpha, ...
    'mu_0', muMod , ...
    'verbosity', 0);

muMod = dummyvar(ca_mod2)' ;
[~,modularityModel2] = wsbm(templateModel.Data.Raw_Data, ...
    size(muMod,1), ...
    'W_Distr', templateModel.W_Distr, ...
    'E_Distr', templateModel.E_Distr, ...
    'alpha', templateModel.Options.alpha, ...
    'mu_0', muMod , ...
    'verbosity', 0);

%%

commsCA = cell(3,1) ;
commsCA{1} = wsbm_community_assign(templateModel) ;
commsCA{2} = ca_mod ;
commsCA{3} = ca_mod2 ;

simMat = zeros(3,3) ;

for idx = 1:3
   for jdx = 1:3
       simMat(idx,jdx) = partition_distance(commsCA{idx},commsCA{jdx}) ;
   end
end

%% viz just matrix

% hemis
hemiVec = ones(154,1) ;
hemiVec(78:end) = 2 ;
hemiCmap = brewermap(2,'PrGn') ;

sp = tight_subplot(1,2,.05,[.01 .01],[.05 .05]);
set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.7, 0.7]);

% orginal weights
axes(sp(1))

plotMat = rawData.CIJ + 1;
plotMat(plotMat == 0) = NaN ;

hh = imagesc(plotMat) ;
% set(hh,'alphadata',~isnan(templateModel.Data.Raw_Data));
axis square
set(gca,'ytick',[])
set(gca,'xtick',[])
% set(gca,'visible','off')

% represent colormap as
hh.CDataMapping = 'direct' ;
colorbar('Ticks',1.5:1:9,'TickLabels',[ NaN rawData.wei ])
cmap = [ 0.5 0.5 0.5 ; flipud(brewermap(7,'spectral')) ] ;

colormap(sp(1),cmap)

set(gca,'FontSize',16)

% hemispheres
hold
viz_comms_on_axes(hemiVec,hemiCmap,16)


% scaled weights
axes(sp(2))

plotMat = rawData.CIJw ;
plotMat(plotMat == 0) = NaN ;

h = imagesc(plotMat) ;
% set(h,'alphadata',~isnan(templateModel.Data.Raw_Data));
axis square
set(gca,'ytick',[])
set(gca,'xtick',[])
% set(gca,'visible','off')

% represent colormap as
h.CDataMapping = 'scaled' ;
colorbar()
colormap(sp(2),[ 1 1 1 ; parula(1000)])

set(gca,'FontSize',16)

% hemispheres
hold
viz_comms_on_axes(hemiVec,hemiCmap,16)


%% viz models

comms = cell(3,1) ;
comms{1} = templateModel ;
comms{2} = modularityModel ;
comms{3} = modularityModel2 ;

for idx = 1:3

    tmpModel = comms{idx} ;
    
    figure
    set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.7, 0.7]);

    k = tmpModel.R_Struct.k ;

    plotWSBM(tmpModel)
    axis square
    xlabel('')
    ylabel('')

    viz_comm_labs_on_axes(wsbm_community_assign(tmpModel),1:8)

    set(gca,'FontSize',18)

end

%% viz info from model

for idx = 1:3

    tmpModel = comms{idx} ;
    k = tmpModel.R_Struct.k ;
    
    figure
    sp = tight_subplot(2,2,.05,[.1 .01],[.1 .1]);
    set(gcf,'Position', [0, 0, 800, 700]);

    axes(sp(1)) ;
    imagesc(reshape(tmpModel.Para.predict_e,k,k)')
    axis square
    caxis([0 1])
    xlabel('predicted edge-existence')
    colorbar
    set(gca,'FontSize',14)

    axes(sp(2)) ;
    imagesc(reshape(tmpModel.Para.predict_w,k,k)')
    axis square
    caxis([0 1])
    xlabel('predicted weight')
    colorbar
    set(gca,'FontSize',14)

    axes(sp(3)) ;
    
    tmpE = tmpModel.Para.predict_e ;
    tmpW = tmpModel.Para.predict_w ;
    tmpE(isnan(tmpE)) = 0 ;
    tmpW(isnan(tmpW)) = 0 ;
    
    scatter(tmpE,tmpW)
    axis square
    cb = colorbar();
    cb.Visible = 'off' ;
    xlabel('predicted edge-existence')
    ylabel('predicted weight')
    xlim([0 1])
    ylim([0 1])
    set(gca,'FontSize',14)
    
    numBoot = 500 ;
    bootVals = bootstrp(numBoot, @(x)corr(x(:,1),x(:,2),'type','spearman') , ...
            [ tmpE tmpW ]) ;
    bootCI = prctile(bootVals,[ 2.5 97.5]) ;
        
    tmpCorr = corr(tmpE,tmpW,'type','spearman') ;
    annotText = { strcat('r:',32,...
        num2str(round(tmpCorr,2,'significant')),32, ...
        '(' ,...
        num2str(round(bootCI(1),2,'significant')) ,'-', ...
        num2str(round(bootCI(2),2,'significant')), ...
        ')') } ;
    text(0.1,0.9,annotText,'FontSize',16,'VerticalAlignment','cap')
    
    axes(sp(4)) ;
    %function [weiBM,avgWeiBM,binBM,avgBinBM,stdWeiBM,sizesMat] = get_block_mat(CIJ,ca,excludeNaN)
    [~,tmp] = get_block_mat_fast(imat,wsbm_community_assign(tmpModel));
    imagesc(tmp)
    colorbar
    caxis([0 1])
    axis square
    xlabel('average weight btwn blocks')
    set(gca,'FontSize',14)

end

%% analyze hemisphere of blockmodel

hemiVec = ones(154,1) ;
hemiVec(78:end) = 2 ;

[~, sortCAidx] = sort(ca_wsbm) ;

k = templateModel.R_Struct.k ;

LIvec = zeros(k,1) ;
for idx = 1:k

    % measure the bi-laterality at each community
    tmpVec = hemiVec(ca_mod == idx) ;
    tmpL = sum(tmpVec == 1) ;
    tmpR = sum(tmpVec == 2) ;
    LIvec(idx) = (tmpL - tmpR) / (tmpL + tmpR) ;

end
%%
