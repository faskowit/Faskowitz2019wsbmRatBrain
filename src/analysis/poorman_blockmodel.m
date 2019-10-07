clc
clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load

config_file='config_template_rb2_oneHemi_analyzeGridRuns.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_baseRes.mat' ] ;
% loads a struct named 'baseRes'
load(loadName) ;

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_consensusCAs.mat' ] ;
load(loadName) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setup some vars

dat = baseRes.rawData ;
dat(isnan(dat)) = 0 ;

nComms = baseRes.wsbm.bestK ;
ind = baseRes.wsbm.bestKind ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get a couple of different similarities

[~,~,matchmat] = matching_ind(dat) ;
matchmat = matchmat + matchmat' ;

[~,~,jacmat] = jaccard_similarity(dat) ; 
jacmat = jacmat + jacmat' ;
jacmat(~~eye(size(jacmat,1))) = 0 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get modules

% function [ comsAcrossG2, qAcrossG2 ] = sweep_gamma_louvain_atK(CIJ,k,gRangeInit,numComsAtK,quiet)
[match_comms,match_qs] = sweep_gamma_louvain_atK(matchmat,nComms) ;
[jac_comms,jac_qs] = sweep_gamma_louvain_atK(jacmat,nComms) ;

% measure some distances
match_vi = partition_distance(match_comms) ;
jac_vi = partition_distance(jac_comms); 

% get some central communities
[~,sortidx] = sort(sum(match_vi),'ascend') ; 
match_centcomm = match_comms(:,sortidx(1)) ;
[~,sortidx] = sort(sum(jac_vi),'ascend') ; 
jac_centcomm = jac_comms(:,sortidx(1)) ;

% get a co-assignment mat
match_agreemat = agreement(match_comms) ./ size(match_comms,2) ;
jac_agreemat = agreement(jac_comms) ./ size(jac_comms,2) ;

wsbm_agreemat = agreement(baseRes.wsbm.ca_K{ind}) ./ size(baseRes.wsbm.ca_K{ind},2) ;
mod_agreemat = agreement(baseRes.mod.ca_K{ind}) ./ size(baseRes.mod.ca_K{ind},2) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% view it

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.4 0.6]);      

mask = logical(triu(ones(size(dat,1)),1)) ;

subplot(221)
imsc_grid_comm(dat,match_centcomm) 
colorbar
axis square
title({'Match community w/' 'rat connectivity'})

subplot(222)
imsc_grid_comm(match_agreemat,match_centcomm)
colorbar
axis square
title({'Match community w/' 'match agreement'})

subplot(223)
imsc_grid_comm(dat,jac_centcomm)
colorbar
axis square
title({'Jac community w/' 'rat connectivity'})

subplot(224)
imsc_grid_comm(jac_agreemat,jac_centcomm)
colorbar
axis square
title({'Jac community w/' 'jac agreement'})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% look at agreement patterns

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.4 0.55]);      

mask = logical(triu(ones(size(dat,1)),1)) ;

subplot(221)
scatter(wsbm_agreemat(mask),match_agreemat(mask))
cc = corr(wsbm_agreemat(mask),match_agreemat(mask),'type','spearman')
xlabel('WSBM agreement') ; ylabel('Match agreement')
annotText = [ '$\rho = ' num2str(round(cc,2)) '$' ] ;
xr = xlim ; yr = ylim ;
text(xr(2)*0.6,.2,annotText,'fontsize',14,'Interpreter','latex')

subplot(222)
scatter(wsbm_agreemat(mask),jac_agreemat(mask))
cc= corr(wsbm_agreemat(mask),jac_agreemat(mask),'type','spearman')
xlabel('WSBM agreement') ; ylabel('Jac agreement')
annotText = [ '$\rho = ' num2str(round(cc,2)) '$' ] ;
xr = xlim ; yr = ylim ;
text(xr(2)*0.7,.2,annotText,'fontsize',14,'Interpreter','latex')

subplot(223)
scatter(mod_agreemat(mask),match_agreemat(mask))
cc = corr(mod_agreemat(mask),match_agreemat(mask),'type','spearman')
xlabel('MOD agreement') ; ylabel('Match agreement')
annotText = [ '$\rho = ' num2str(round(cc,2)) '$' ] ;
xr = xlim ; yr = ylim ;
text(xr(2)*0.6,.2,annotText,'fontsize',14,'Interpreter','latex')

subplot(224)
scatter(mod_agreemat(mask),jac_agreemat(mask))
cc= corr(mod_agreemat(mask),jac_agreemat(mask),'type','spearman')
xlabel('MOD agreement') ; ylabel('Jac agreement')
annotText = [ '$\rho = ' num2str(round(cc,2)) '$' ] ;
xr = xlim ; yr = ylim ;
text(xr(2)*0.6,.2,annotText,'fontsize',14,'Interpreter','latex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% lets bootstrap some coors
rng(123)

nPerm = 10000 ;
agree_coors = zeros(nPerm,4) ; 

maskinds = find(mask) ;
for idx = 1:nPerm
   disp(idx)
   bootmask = randsample(maskinds,floor(length(maskinds)/4),'true') ;
   
   agree_coors(idx,1) = corr(wsbm_agreemat(bootmask),match_agreemat(bootmask),'type','spearman');
   agree_coors(idx,2) = corr(wsbm_agreemat(bootmask),jac_agreemat(bootmask),'type','spearman');
   agree_coors(idx,3) = corr(mod_agreemat(bootmask),match_agreemat(bootmask),'type','spearman');
   agree_coors(idx,4) = corr(mod_agreemat(bootmask),jac_agreemat(bootmask),'type','spearman');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

aa = abs(wsbm_agreemat(mask) - match_agreemat(mask)) ; 
bb = abs(wsbm_agreemat(mask) - jac_agreemat(mask)) ; 
cc = abs(mod_agreemat(mask) - match_agreemat(mask)) ; 
dd = abs(mod_agreemat(mask) - jac_agreemat(mask)) ; 

histogram(aa) ;
hold on
histogram(bb) ; 
histogram(cc) ;
histogram(dd) ;


%%
m_nonNanBl = dummyvar(cons_ca.mod)' * ~isnan(baseRes.rawData) * dummyvar(cons_ca.mod) ;
w_nonNanBl = dummyvar(cons_ca.wsbm)' * ~isnan(baseRes.rawData) * dummyvar(cons_ca.wsbm) ;

nanE = isnan(baseRes.rawData) ;
numNonNan = length(baseRes.rawData(~nanE)) ;

[~,wSort] = sort(cons_ca.wsbm) ;
[~,mSort] = sort(cons_ca.mod) ;

wsbm_blmat = get_block_mat(entr_K.wsbm.ent{ind}, cons_ca.wsbm) ;
mod_blmat = get_block_mat(entr_K.mod.ent{ind}, cons_ca.mod) ;

% mod edge probs
m_od = mean(motifAna.mod.motifEdgeMats{ind}.odMat>0,3) ;
m_assort = mean(motifAna.mod.motifEdgeMats{ind}.aMat>0,3) ;

[~,m_od_bl] = get_block_mat(m_od,cons_ca.mod) ;

%%

avg_nonz_val = zeros(nComms,1) ;
bbb2 = zeros(nComms,1) ;

for idx = 1:nComms

    tmpInd = cons_ca.mod == idx ;
    
    tmp1 = m_od(tmpInd,tmpInd) ;
    tmp2 = nanE(tmpInd,tmpInd) ;
    
    avg_nonz_val(idx) = mean(tmp1(~tmp2)) ;
    
end

%%

w_od = mean(motifAna.wsbm.motifEdgeMats{ind}.odMat>0,3) ;
w_assort = mean(motifAna.wsbm.motifEdgeMats{ind}.aMat>0,3) ;
w_core = mean(motifAna.wsbm.motifEdgeMats{ind}.cMat>0,3) ;
w_per = mean(motifAna.wsbm.motifEdgeMats{ind}.pMat>0,3) ;
w_dis = mean(motifAna.wsbm.motifEdgeMats{ind}.dMat>0,3) ;

%%

[~,m_ent_bl] = get_block_mat(entr_K.mod.ent{ind}, cons_ca.mod)
[~,w_ent_bl] = get_block_mat(entr_K.wsbm.ent{ind}, cons_ca.wsbm)

[~,w_avg_bl] = get_block_mat(dat,cons_ca.wsbm) ;
[~,m_avg_bl] = get_block_mat(dat,cons_ca.mod) ;


