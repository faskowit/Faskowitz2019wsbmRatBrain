
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
%%

dat = baseRes.rawData ;
dat(isnan(dat)) = 0 ;
ind = baseRes.wsbm.bestKind ;

%%

wsbm_ca = baseRes.wsbm.ca_K{ind} ;
mod_ca = baseRes.mod.ca_K{ind} ;

%wsbm_diffs = zeros(baseRes.wsbm.bestK,baseRes.wsbm.bestK,size(wsbm_ca,2)) ;
%mod_diffs = zeros(baseRes.wsbm.bestK,baseRes.wsbm.bestK,size(wsbm_ca,2)) ;
wsbm_diffs = zeros(45,size(wsbm_ca,2)) ;
mod_diffs = zeros(45,size(wsbm_ca,2)) ;

wsbm_sumdiff = zeros(size(wsbm_ca,2),1) ;
mod_sumdiff = zeros(size(wsbm_ca,2),1) ;

upp_mask = logical(triu(ones(baseRes.wsbm.bestK),1)) ;
low_mask = logical(tril(ones(baseRes.wsbm.bestK),-1)) ;

for idx = 1:size(wsbm_ca,2)
    
    tmpca_1 = wsbm_ca(:,idx) ;
    tmpca_2 = mod_ca(:,idx) ; 

    [~,tmpMat_1] = get_block_mat(dat,tmpca_1) ;
    [~,tmpMat_2] = get_block_mat(dat,tmpca_2) ;    
    
    wsbm_diffs(:,idx) = tmpMat_1(upp_mask) - tmpMat_1(low_mask) ;
    mod_diffs(:,idx) = tmpMat_2(upp_mask) - tmpMat_2(low_mask) ;
    
    wsbm_sumdiff(idx) = sum(tmpMat_1(upp_mask)) - sum(tmpMat_1(low_mask)) ;
    mod_sumdiff(idx) = sum(tmpMat_2(upp_mask)) - sum(tmpMat_2(low_mask)) ;
    
end

%%

wsbm_meanDiffs = mean(wsbm_diffs,2) ;
mod_meanDiffs = mean(mod_diffs,2) ;

%% boot the difference

rng(42)

nPerm = 10000 ; 
w_permRes = zeros(45,nPerm) ;
m_permRes = zeros(45,nPerm) ;
sz = size(wsbm_ca,2) ;
sumdiff_diff = zeros(nPerm,1) ;

for idx = 1:nPerm 
    
    disp(idx)

    randinds = datasample(1:sz,sz) ;

    w_permRes(:,idx) = mean(wsbm_diffs(:,randinds),2) ;
    m_permRes(:,idx) = mean(mod_diffs(:,randinds),2) ;

    sumdiff_diff(idx) = std(wsbm_sumdiff(randinds)) - std(mod_sumdiff(randinds)) ;
    
end

%% pvalus on the difference in ratios

pval = (1+ sum(sumdiff_diff<0)) / (1+nPerm) ;


%% get the percentiles that don't cross 0

w_p = prctile(w_permRes',[1 99]) ;
m_p = prctile(m_permRes',[1 99]) ;

w_nonz = (w_p(1,:) .* w_p(2,:)) > 0 ;
m_nonz = (m_p(1,:) .* m_p(2,:)) > 0 ;

w_asym = nan(10) ; m_asym = nan(10) ;

w_asym(upp_mask) = w_nonz' .* wsbm_meanDiffs ;
m_asym(upp_mask) = m_nonz' .* mod_meanDiffs ;

w_asym(w_asym==0) = nan ;
m_asym(m_asym==0) = nan ;

% w_nonz(~w_nonz) = nan ;
% m_nonz(~m_nonz) = nan ;

%%

cm = brewermap(50,'PuOr') ;

set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.42, 0.5]);
ts = tight_subplot(1,2,0.05) ;

cax = [ -0.2 0.2 ] ;

axes(ts(1)) 
h = imagesc(w_asym)
set(h,'alphadata',~isnan(w_asym));
axis square ; colorbar; colormap(cm) ; caxis(cax)
title('WSBM')

axes(ts(2)) 
h = imagesc(m_asym)
set(h,'alphadata',~isnan(m_asym));
axis square ; colorbar ; colormap(cm) ; caxis(cax)
title('Modular')




