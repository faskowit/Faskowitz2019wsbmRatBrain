
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

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_motifAna.mat' ] ;
load(loadName) ;

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_consensusCAs.mat' ] ;
load(loadName) ;

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_netstatsRes.mat' ] ;
load(loadName) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% entropy across K

commNames = { 'wsbm' 'mod' } ;

entr_K = struct() ;
entr_K.wsbm.ent = cell(size(baseRes.wsbm.ca_K,1),1) ;
entr_K.mod.ent = cell(size(baseRes.wsbm.ca_K,1),1) ;

entr_K.wsbm.sum = cell(size(baseRes.wsbm.ca_K,1),1) ;
entr_K.mod.sum = cell(size(baseRes.wsbm.ca_K,1),1) ;

for idx = 1:length(entr_K.wsbm.ent) 

    for cn = 1:length(commNames)
    
        tmpMot = motifAna.(commNames{cn}).motifEdgeMats{idx} ;

        entr_K.(commNames{cn}).ent{idx} = get_comm_motif_entropy(mean(tmpMot.aMat>0,3),...
                                mean(tmpMot.cMat>0,3),...
                                mean(tmpMot.pMat>0,3),...
                                mean(tmpMot.dMat>0,3),...
                                mean(tmpMot.odMat>0,3) ) ;    
        entr_K.(commNames{cn}).sum{idx} = ( sum(entr_K.(commNames{cn}).ent{idx})' + ...
                                sum(entr_K.(commNames{cn}).ent{idx},2) ) ./ 2 ;
       
    end                                         
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% lets look at the edge-wise vals

nanEdges = isnan(baseRes.rawData) ;

sss = sortedInd(cons_ca.wsbm) ;
sss2 = sortedInd(cons_ca.mod) ;
% sss = 1:77 ;
% sss2 = 1:77 ;

curScale = baseRes.wsbm.bestKind ;
curScale = 9 ;

edgeMats = { 'odMat' 'aMat' 'cMat' 'pMat' 'dMat' } ;

for idx = 1:length(edgeMats)

    tmp = mean(motifAna.wsbm.motifEdgeMats{curScale}.(edgeMats{idx}) > 0,3) ;
    subplot(2,6,idx) ; 
    h = imagesc(tmp(sss,sss)); caxis([0 1])
    %set(h,'alphadata',nanEdges==0);
end
subplot(266) ; 
%h = imagesc( entr_K.wsbm.ent{curScale}(sss,sss) ) ; colorbar
h = imsc_grid_comm( entr_K.wsbm.ent{curScale}, cons_ca.wsbm ) ; colorbar
%set(h,'alphadata',nanEdges(sss,sss)==0);

for idx = 1:length(edgeMats)

    tmp = mean(motifAna.mod.motifEdgeMats{curScale}.(edgeMats{idx}) > 0,3) ;
    subplot(2,6,idx+6) ; 
    h= imagesc(tmp(sss2,sss2)); caxis([0 1])
    %set(h,'alphadata',nanEdges==0);
end
subplot(2,6,12) ;  
%h = imagesc( entr_K.mod.ent{curScale}(sss2,sss2) ) ; colorbar
h = imsc_grid_comm( entr_K.mod.ent{curScale}, cons_ca.mod ) ; colorbar

%set(h,'alphadata',nanEdges(sss2,sss2)==0);

%%

dat = baseRes.rawData ;
dat(isnan(dat)) = 0 ;

weis = strengths_dir(dat) ;
clustc = clustering_coef_wd(dat) ;

wsbmSum = ( sum(entr_K.wsbm{curScale})' + sum(entr_K.wsbm{curScale},2) ) ./ 2 ;
modSum = ( sum(entr_K.mod{curScale})' + sum(entr_K.mod{curScale},2) ) ./ 2  ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

oo = cellfun(@(x)sum(sum(x)),entr_K.mod,'UniformOutput',true) ;

iii = [ 1 REASONABLE_COM_RANGE_IND] ;

scatter(iii,oo(iii))

%%

mmm = zeros(length(REASONABLE_COM_RANGE_IND),length(REASONABLE_COM_RANGE_IND)) ;
mmm2 = zeros(length(REASONABLE_COM_RANGE_IND),length(REASONABLE_COM_RANGE_IND)) ;


for idx = 1:length(REASONABLE_COM_RANGE_IND)
    
    for jdx = 1:length(REASONABLE_COM_RANGE_IND)
    
        ind = REASONABLE_COM_RANGE_IND(idx) ;

        w_ent1 = entr_K.wsbm.sum{idx} ;
        m_ent1 = entr_K.mod.sum{idx} ;
        
        w_ent2 = entr_K.wsbm.sum{jdx} ;
        m_ent2 = entr_K.mod.sum{jdx} ;
 
        mmm(idx,jdx) = corr(w_ent1,w_ent2) ;
        mmm2(idx,jdx) = corr(m_ent1,m_ent2) ;
       
    end
end

%%

curScale = baseRes.wsbm.bestKind ;

ag_curScale = partition_distance(baseRes.wsbm.ca_K{curScale}) ; 

cc = consensus_und(vi_curScale,0.1,100)

%%

wsbm_od_probs = mean(motifAna.wsbm.motifEdgeMats{10}.odMat > 0,3) ;
mod_od_probs = mean(motifAna.mod.motifEdgeMats{10}.odMat > 0,3) ;

wsbm_a_probs = mean(motifAna.wsbm.motifEdgeMats{10}.aMat > 0,3) ;
mod_a_probs = mean(motifAna.mod.motifEdgeMats{10}.aMat > 0,3) ;

dat = CIJ .* 1; 
dat(isnan(dat)) = 0 ;

dd = strengths_dir(dat) ;

subplot(221)
imagesc(wsbm_od_probs)
subplot(222)
imagesc(mod_od_probs)
subplot(223)
imagesc(wsbm_a_probs)
subplot(224)
imagesc(mod_a_probs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

curScale = motifAna.wsbm.motifEdgeMats{baseRes.wsbm.bestKind} ;

ttt = get_comm_motif_entropy(mean(curScale.aMat>0,3),mean(curScale.cMat>0,3),...
                             mean(curScale.pMat>0,3),mean(curScale.dMat>0,3)) ;

ttt2 = get_comm_motif_entropy(mean(curScale.aMat>0,3),mean(curScale.cMat>0,3),...
                             mean(curScale.pMat>0,3),mean(curScale.dMat>0,3),...
                             mean(curScale.odMat>0,3)) ;                         
      
ttt2(isnan(ttt2)) = 0 ;

[~,~,mmm] = matching_ind(dat) ;
mmm = mmm + mmm' ;                         

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

pp = mean(curScale.aMat>0,3) + mean(curScale.cMat>0,3) + mean(curScale.pMat>0,3) + mean(curScale.dMat>0,3) + mean(curScale.odMat>0,3)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now lets compare each profile 
% this analysis does not really show anything... 

% blues_cm = brewermap(19,'blues') ;
% reds_cm = brewermap(19,'reds') ;
% purps_cm = brewermap(19,'purples');
% greens_cm = brewermap(19,'greens') ;
% 
% dat = CIJ ;
% dat(isnan(dat)) = 0 ;
% rawDataStr = strengths_dir(dat) ;
% [weiSort,strSortInd] = sort(rawDataStr,'descend') ;
% 
% figure 
% hold
% for idx = 1:19
%     
%     meanAssort = mean(motifAna.motifStructs{idx}.assort,2) ;
%     plot(meanAssort(strSortInd),'Color',blues_cm(idx,:)) ;
% 
%     meanAssort = mean(motifAna.motifStructs{idx}.core,2) ;
%     plot(meanAssort(strSortInd),'Color',reds_cm(idx,:)) ;
% 
%     meanAssort = mean(motifAna.motifStructs{idx}.peri,2) ;
%     plot(meanAssort(strSortInd),'Color',purps_cm(idx,:)) ;
% 
%     meanAssort = mean(motifAna.motifStructs{idx}.disort,2) ;
%     plot(meanAssort(strSortInd),'Color',greens_cm(idx,:)) ;
% 
%     
% end
% 
% meanCore = mean(motifAna.motifStructs{12}.core,2) ;
% meanPeri = mean(motifAna.motifStructs{12}.peri,2) ;
% meanDisort = mean(motifAna.motifStructs{12}.disort,2) ;

%%

a = zeros(154,19) ;
b = zeros(154,19) ;
c = zeros(154,19) ;
d = zeros(154,19) ;
e = zeros(154,19) ;

for idx = 1:19

    a(:,idx) = mean(motifAna.motifStructs{idx}.assort,2) ;
    b(:,idx) = mean(motifAna.motifStructs{idx}.core,2) ;
    c(:,idx) = mean(motifAna.motifStructs{idx}.peri,2) ;
    d(:,idx) = mean(motifAna.motifStructs{idx}.disort,2) ;
    e(:,idx) = -1 .* ( a(:,idx).*log2(a(:,idx)) + ...
        b(:,idx).*log2(b(:,idx)) + ...
        c(:,idx).*log2(c(:,idx)) + ...
        d(:,idx).*log2(d(:,idx)) ) ; 
end
e(isnan(e)) = 0 ;

[~,entSort] = sort(mean(e,2)) ;

spect_cm = brewermap(19,'spectral') ;

% plot it
figure
set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.6, 0.6]);

hold

for idx = 1:19

    plot(e(entSort,idx),'Color',spect_cm(idx,:));
    
end

ylabel('Motif Entropy')

cb = colorbar('Ticks',0.5:1:18.5,'TickLabels', 2:20 ) ;
colormap(spect_cm)
caxis([0 19])
cb.Label.String = 'Numeber of communities (k)' ;
set(gca,'FontSize',16)

% ddd = load('data/raw/CCTX_matrix.mat') ;
% ddd.area_names
% 
% ppp = gca ;
% ppp.XTick = 1:2:154 ;
% xticklabels(ddd.area_names)
% 
% xtickangle(90)

figure
hold
tmp = mean(a,2) ;
plot(tmp(strSortInd),'Color',blues_cm(10,:)) ;

tmp = mean([ b c ],2) ;
plot(tmp(strSortInd),'Color',reds_cm(10,:)) ;

% tmp = mean(c,2) ;
% plot(tmp(strSortInd),'Color',purps_cm(10,:)) ;

tmp = mean(d,2) ;
plot(tmp(strSortInd),'Color',greens_cm(10,:)) ;

legend({ 'assort.' 'core-periphery' 'dissassort.'})

ylabel('Probability')
xlabel('Nodes, by increasing degree')

%%

config_file='config_template_analyzeGridRuns.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_consensusRuns.mat' ] ;
bothHemiConsensus = load(loadName) ;

BH_templateModel = bothHemiConsensus.consensusCent ;

%% wsbm

ca_wsbm = wsbm_community_assign(BH_templateModel) ;

[X,Y,INDSORT] = grid_communities(ca_wsbm) ;

entropyMat = squareform(pdist(e,'cityblock')) ;

figure
set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.4, 0.8]);

imagesc(entropyMat(INDSORT,INDSORT)) ;
hold
plot(X,Y,'r','linewidth',2)
axis square
colorbar

set(gca,'ytick',[])
set(gca,'xtick',[])
set(gca,'visible','off')
set(gca,'FontSize',16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% modular

[ca_mod_BH,ca_modQ] = gsweep_mod_atK_dir(datStruct.both ,size(BH_templateModel.Options.mu_0,1)) ;

ca_mod_BH = hungarianMatch(ca_wsbm,ca_mod_BH) ;
[X,Y,INDSORT] = grid_communities(ca_mod_BH) ;

entropyMat = squareform(pdist(e,'cityblock')) ;

figure
set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.4, 0.8]);

imagesc(entropyMat(INDSORT,INDSORT)) ;
hold
plot(X,Y,'r','linewidth',2)
axis square

colorbar

set(gca,'ytick',[])
set(gca,'xtick',[])
set(gca,'visible','off')
set(gca,'FontSize',16)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save the data at the very least!!

% todo