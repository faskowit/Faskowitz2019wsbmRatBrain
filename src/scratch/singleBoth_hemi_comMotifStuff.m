
clc
clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this is where i load data

% addpath(genpath('/home/jfaskowi/JOSHSTUFF/projects/ratbrain'))

config_file='config_template_rb2_analyzeGridRuns.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_baseRes.mat' ] ;
bothHemi = load(loadName) ;
ls(loadName)

config_file='config_template_rb2_oneHemi_analyzeGridRuns.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_baseRes.mat' ] ;
oneHemi = load(loadName) ;
ls(loadName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this is where i analyze data

hemis = { 'both' 'sing' } ;
numNodes = { 154 77 } ;

inData = struct() ;
inData.both = bothHemi.baseRes ;
inData.sing = oneHemi.baseRes ;

% %% get in data that is top 50
% 
% % getLogEvid = @(x) x.Para.LogEvidence ;
% 
% inDataTop = struct() ; 
% inDataTop.sing.logEvid = cell(19,1) ;
% inDataTop.both.logEvid = cell(19,1) ;
% inDataTop.sing.evidSort = cell(19,1) ;
% inDataTop.both.evidSort = cell(19,1) ;
% 
% for idx = 1:19
% 
%     inDataTop.sing.logEvid{idx} = inData.sing.logEvid_K{idx} ;
%     inDataTop.both.logEvid{idx} = inData.both.logEvid_K{idx} ;
%     
%     [~,inDataTop.sing.evidSort{idx}] = sort(inDataTop.sing.logEvid{idx},'descend') ;
%     [~,inDataTop.both.evidSort{idx}] = sort(inDataTop.both.logEvid{idx},'descend') ;   
% 
% %     tmp = inData.sing.kResults{idx}(inDataTop.sing.evidSort{idx}) ;
% %     inDataTop.sing.kResults{idx} = tmp(1:40) ;
% %     tmp = inData.both.kResults{idx}(inDataTop.both.evidSort{idx}) ;
% %     inDataTop.both.kResults{idx} = tmp(1:40) ;
% 
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% gather the data into a struct
 
% getLogEvid = @(x) x.Para.LogEvidence ;

idatstruct = inData ;

gather = struct() ;
gather.both = struct() ;
gather.sing = struct() ;

for hemi = 1:length(hemis)

    numLevels = size(idatstruct.(hemis{hemi}).logEvid_K,1) ;
    
    % initialize variables
    gather.(hemis{hemi}).wsbmCon = zeros(numNodes{hemi},numNodes{hemi},numLevels) ;
    gather.(hemis{hemi}).versAcrossK = zeros(numNodes{hemi},numLevels) ;

    % loop across all the kResults for this hemi
    for idx = 1:numLevels
 
        % already gathered
        gather.(hemis{hemi}).wsbmComs{idx} = idatstruct.(hemis{hemi}).ca_K{idx} ;
        gather.(hemis{hemi}).logEvid{idx} = idatstruct.(hemis{hemi}).logEvid_K{idx} ;
        
        numRuns = size(idatstruct.(hemis{hemi}).logEvid_K{idx},1) ;
        
        % recover the agreement matrix (make into probability by dividing
        % by number of algo. repeitions
        gather.(hemis{hemi}).wsbmCon(:,:,idx) = ...
            agreement(gather.(hemis{hemi}).wsbmComs{idx}) ./ ...
            size(gather.(hemis{hemi}).wsbmComs{idx},2) ;

        % versatility across the community assignments 
        gather.(hemis{hemi}).versAcrossK(:,idx) = ...
            get_nodal_versatility(gather.(hemis{hemi}).wsbmComs{idx});
           
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do the motif analysis

motifVarNAmes = { 'assort' 'core' 'peri' 'disort' } ; 

ddd = load('data/raw/CCTX_matrix.mat') ;

datStruct = struct() ;
datStruct.sing = ddd.CIJw(1:77,1:77) ;
datStruct.both = ddd.CIJw ;

idatstruct = inData ;

motifAna = struct() ; 
motifAna.both = struct() ;
motifAna.sing = struct() ;

for hemi = 1:length(hemis)
    
    numLevels = size(idatstruct.(hemis{hemi}).ca_K,1) ;
    
    % allocate cell array
    motifAna.(hemis{hemi}).blMatStacks = cell([numLevels 1]) ;
    motifAna.(hemis{hemi}).comMotifMat = cell([numLevels 1]) ;
    motifAna.(hemis{hemi}).motifStructs = cell([numLevels 1]) ;
    
    % loop across all the kResults for this hemi
    for idx = 1:length(gather.(hemis{hemi}).wsbmComs)

        disp([ 'idx: ' int2str(idx) ])

        numReps = size(idatstruct.(hemis{hemi}).logEvid_K{idx},1) ;
        
        % allocate kxk matrix at this idx (k = idx+1)
        motifAna.(hemis{hemi}).blMatStacks{idx} =  ...
            zeros(idx+1,idx+1,numReps) ;
        
        % allocate kxk com motif matrix
        motifAna.(hemis{hemi}).comMotifMat{idx} =  ...
            zeros(idx+1,idx+1,numReps) ;
        
        % allocate com motif struct
        motifAna.(hemis{hemi}).motifStructs{idx} =  ...
             struct('assort',zeros([numNodes{hemi} numReps]), ...
                    'core',zeros([numNodes{hemi} numReps]),...
                    'peri',zeros([numNodes{hemi} numReps]),...
                    'disort',zeros([numNodes{hemi} numReps])...
                    ) ;
        tmp = fieldnames(motifAna.(hemis{hemi}).motifStructs{idx}) ;
        motifAna.(hemis{hemi}).motifStructs{idx}.varNames = tmp ;
        
        % get the stack
        for jdx = 1:size(gather.(hemis{hemi}).wsbmComs{idx},2)
        
            [~,tmp] = get_block_mat(datStruct.(hemis{hemi}),...
                gather.(hemis{hemi}).wsbmComs{idx}(:,jdx)) ;
            
            % get the avg block matrix
            motifAna.(hemis{hemi}).blMatStacks{idx}(:,:,jdx) = tmp ;
            motifAna.(hemis{hemi}).comMotifMat{idx}(:,:,jdx) = ...
                wsbm_offDiag_motif(tmp) ;
            
        end 
    
        % get the node-wise rate of com motif participation
        for jdx = 1:size(gather.(hemis{hemi}).wsbmComs{idx},2)
            
            tmpMat = motifAna.(hemis{hemi}).comMotifMat{idx}(:,:,jdx) ;
            
            % for 1:4, assortative, cor, periph, dis
            for commMotifInd = 1:length(motifVarNAmes)
                 % map the motif participation to the nodes (k = idx+1)
                 k = (idx+1) ;
                 for kdx = 1:k
                 
                    % select the nodes in the k community
                    tmpCommSelect = (kdx == gather.(hemis{hemi}).wsbmComs{idx}(:,jdx)) ;
                    
                    % compute proportion that community participates in this
                    % specific community motif
                    tmpMotifVal = ( sum(tmpMat(kdx,:) == commMotifInd) + ...
                        sum(tmpMat(:,kdx) == commMotifInd) ) / ...
                        ((k-1)*2) ;
                    
                    % get proportion by dividing by total off-diag blocked consider
                    tmpName = motifVarNAmes{commMotifInd} ;
                    motifAna.(hemis{hemi}).motifStructs{idx}.(tmpName)(tmpCommSelect,jdx) = ...
                        tmpMotifVal ;
                
                 end % kdx
            end % commMotifInd
            
        end
        
    end % 1:length(gather.(hemis{hemi}).wsbmComs)
end %  hemi = 1:length(hemis)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now lets compare each profile 

blues_cm = brewermap(19,'blues') ;
reds_cm = brewermap(19,'reds') ;
purps_cm = brewermap(19,'purples');
greens_cm = brewermap(19,'greens') ;

rawDataStr = strengths_dir(ddd.CIJw) ;
[weiSort,strSortInd] = sort(rawDataStr,'descend') ;

figure 
hold
for idx = 1:19
    
    meanAssort = mean(motifAna.both.motifStructs{idx}.assort,2) ;
    plot(meanAssort(strSortInd),'Color',blues_cm(idx,:)) ;

    meanAssort = mean(motifAna.both.motifStructs{idx}.core,2) ;
    plot(meanAssort(strSortInd),'Color',reds_cm(idx,:)) ;

    meanAssort = mean(motifAna.both.motifStructs{idx}.peri,2) ;
    plot(meanAssort(strSortInd),'Color',purps_cm(idx,:)) ;

    meanAssort = mean(motifAna.both.motifStructs{idx}.disort,2) ;
    plot(meanAssort(strSortInd),'Color',greens_cm(idx,:)) ;

    
end

meanCore = mean(motifAna.both.motifStructs{12}.core,2) ;
meanPeri = mean(motifAna.both.motifStructs{12}.peri,2) ;
meanDisort = mean(motifAna.both.motifStructs{12}.disort,2) ;

%%

a = zeros(154,19) ;
b = zeros(154,19) ;
c = zeros(154,19) ;
d = zeros(154,19) ;
e = zeros(154,19) ;

for idx = 1:19

    a(:,idx) = mean(motifAna.both.motifStructs{idx}.assort,2) ;
    b(:,idx) = mean(motifAna.both.motifStructs{idx}.core,2) ;
    c(:,idx) = mean(motifAna.both.motifStructs{idx}.peri,2) ;
    d(:,idx) = mean(motifAna.both.motifStructs{idx}.disort,2) ;
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




