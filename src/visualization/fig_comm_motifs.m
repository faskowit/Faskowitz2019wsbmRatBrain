
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

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_consensusRuns.mat' ] ;
load(loadName) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% lets look at the edge-wise vals

ca_wsbm = wsbm_community_assign(consensusCent) ;
ca_mod = baseRes.bestKcentmod ;

sort_wsbm = sortedInd(ca_wsbm) ;
sort_mod = sortedInd(ca_mod) ;

wsbm_od_probs = mean(motifAna.motifEdgeMats{10}.odMat > 0,3) ;
mod_od_probs = mean(modMotifAna.motifEdgeMats{10}.odMat > 0,3) ;

wsbm_a_probs = mean(motifAna.motifEdgeMats{10}.aMat > 0,3) ;
mod_a_probs = mean(modMotifAna.motifEdgeMats{10}.aMat > 0,3) ;

dat = CIJ .* 1; 
dat(isnan(dat)) = 0 ;

dd = strengths_dir(dat) ;

subplot(221)
imagesc(wsbm_od_probs(sort_wsbm,sort_wsbm))
subplot(222)
imagesc(mod_od_probs(sort_mod,sort_mod))
subplot(223)
imagesc(wsbm_a_probs(sort_wsbm,sort_wsbm))
subplot(224)
imagesc(mod_a_probs(sort_mod,sort_mod))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

curScale = motifAna.motifEdgeMats{10} ;

ttt = get_comm_motif_entropy(mean(curScale.aMat>0,3),mean(curScale.cMat>0,3),...
                             mean(curScale.pMat>0,3),mean(curScale.dMat>0,3)) ;

ttt2 = get_comm_motif_entropy(mean(curScale.aMat>0,3),mean(curScale.cMat>0,3),...
                             mean(curScale.pMat>0,3),mean(curScale.dMat>0,3),...
                             mean(curScale.odMat>0,3)) ;                         
      
ttt2(isnan(ttt2)) = 0 ;

[~,~,mmm] = matching_ind(dat) ;
mmm = mmm + mmm' ;                         

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