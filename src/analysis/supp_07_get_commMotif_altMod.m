
clc
clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this is where i load data

config_file='config_template_rb2_oneHemi_analyzeGridRuns.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_altMod_baseRes.mat' ] ;
% loads a struct named 'baseRes'
load(loadName) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% gather the data into a struct
%  
% gather = struct() ;
% 
% numLevels = size(baseRes.wsbm.logEvid_K,1) ;
% 
% % initialize variables
% gather.wsbmCon = zeros(NUM_NODES,NUM_NODES,numLevels) ;
% gather.versAcrossK = zeros(NUM_NODES,numLevels) ;
% 
% % loop across all the kResults for this hemi
% for idx = 1:numLevels
% 
%     % already gathered
%     gather.wsbmComs{idx} = baseRes.wsbm.ca_K{idx} ;
%     gather.logEvid{idx} = baseRes.wsbm.logEvid_K{idx} ;
% 
%     numRuns = size(baseRes.wsbm.logEvid_K{idx},1) ;
% 
%     % recover the agreement matrix (make into probability by dividing
%     % by number of algo. repeitions
%     gather.wsbmCon(:,:,idx) = ...
%         agreement(gather.wsbmComs{idx}) ./ ...
%         size(gather.wsbmComs{idx},2) ;
% 
%     % versatility across the community assignments 
%     gather.versAcrossK(:,idx) = ...
%         get_nodal_versatility(gather.wsbmComs{idx});
% 
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make some rand communities

rng(123)

comTypes = { 'mod' 'wsbm' } ;
nComTypes = { 'randM' 'randW' };

nNodes = size(baseRes.rawData,1);

for cdx = 1:length(comTypes)

    baseRes.(nComTypes{cdx}).ca_K = cell(length(comTypes),1) ;
    
    % randomize each community strucutre
    for idx = 1:length(REASONABLE_COM_RANGE_IND)

        % get the ca_K struct
        tmpCAs = baseRes.(comTypes{cdx}).ca_K{idx} ;
        randCAs = zeros(size(tmpCAs));

        for jdx = 1:size(tmpCAs,2)
            randp = randperm(nNodes) ;
            randCAs(:,jdx) = tmpCAs(randp,idx) ;
        end

        baseRes.(nComTypes{cdx}).ca_K{idx} = randCAs ;
        
    end
end

% % also randomize based on shuffling labels
% 
% comTypes = { 'wsbm' 'mod' } ;
% nComTypes = { 'randblW' 'randblM' };
% 
% nNodes = size(baseRes.rawData,1);
% 
% for cdx = 1:length(comTypes)
% 
%     baseRes.(nComTypes{cdx}).ca_K = cell(length(comTypes),1) ;
%     
%     % randomize each community strucutre
%     for idx = 1:size(baseRes.wsbm.ca_K,1)
%         
%         % get the ca_K struct
%         tmpCAs = baseRes.(comTypes{cdx}).ca_K{idx} ;
%         randCAs = zeros(size(tmpCAs));
% 
%         nComms = length(unique(tmpCAs(:,1))) ;
%         
%         for jdx = 1:size(tmpCAs,2)
%             randp = randperm(nComms) ;
%             tmp = tmpCAs(:,idx) ;
%             for kdx = 1:nComms
%                tmp(tmp==kdx) = randp(kdx) ; 
%             end
%             
%             randCAs(:,jdx) = tmpCAs(:,idx) ;
%         end
% 
%         baseRes.(nComTypes{cdx}).ca_K{idx} = randCAs ;
%         
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do the motif analysis

comTypes = { 'mod' 'randM' 'wsbm' 'randW' } ;
motifVarNAmes = { 'assort' 'core' 'peri' 'disort' } ; 
CIJ = baseRes.rawData ;

% allocate struct
motifAna = struct() ; 

% numLevels = size(baseRes.wsbm.ca_K,1) ;
numLevels = length(REASONABLE_COM_RANGE_IND) ;

for cdx = 1:length(comTypes)

    % allocate cell array
    motifAna.(comTypes{cdx}).blMatStacks = cell([numLevels 1]) ;
    motifAna.(comTypes{cdx}).comMotifMat = cell([numLevels 1]) ;
    motifAna.(comTypes{cdx}).motifStructs = cell([numLevels 1]) ;
    motifAna.(comTypes{cdx}).motifEdgeMats = cell([numLevels 1]) ;

    % loop across all the kResults for this hemi
    for idx = 1:numLevels

        disp([ 'idx: ' int2str(idx) ])

        numReps = size(baseRes.(comTypes{cdx}).ca_K{idx},1) ;

        % allocate kxk matrix at this idx (k = idx+1)
        motifAna.(comTypes{cdx}).blMatStacks{idx} =  ...
            zeros(idx+1,idx+1,numReps) ;

        % allocate kxk com motif matrix
        motifAna.(comTypes{cdx}).comMotifMat{idx} =  ...
            zeros(idx+1,idx+1,numReps) ;

        % allocate com motif struct
        motifAna.(comTypes{cdx}).motifStructs{idx} =  ...
             struct('assort',zeros([NUM_NODES numReps]), ...
                    'core',zeros([NUM_NODES numReps]),...
                    'peri',zeros([NUM_NODES numReps]),...
                    'disort',zeros([NUM_NODES numReps])...
                    ) ;
        tmpName = fieldnames(motifAna.(comTypes{cdx}).motifStructs{idx}) ;
        motifAna.(comTypes{cdx}).motifStructs{idx}.varNames = tmpName ;

        % allocate the com motif edge-wise mats
        motifAna.(comTypes{cdx}).motifEdgeMats{idx}.aMat = zeros([NUM_NODES NUM_NODES numReps]) ;
        motifAna.(comTypes{cdx}).motifEdgeMats{idx}.cMat = zeros([NUM_NODES NUM_NODES numReps]) ;
        motifAna.(comTypes{cdx}).motifEdgeMats{idx}.pMat = zeros([NUM_NODES NUM_NODES numReps]) ;
        motifAna.(comTypes{cdx}).motifEdgeMats{idx}.dMat = zeros([NUM_NODES NUM_NODES numReps]) ;
        motifAna.(comTypes{cdx}).motifEdgeMats{idx}.odMat = zeros([NUM_NODES NUM_NODES numReps]) ;

        % get the stack
        for jdx = 1:size(baseRes.(comTypes{cdx}).ca_K{idx},2)
            
            comInds = baseRes.(comTypes{cdx}).ca_K{idx}(:,jdx) ;
            [~,tmpBl] = get_block_mat(CIJ,comInds) ;
            % get the avg block matrix
            motifAna.(comTypes{cdx}).blMatStacks{idx}(:,:,jdx) = tmpBl ;

            % the motif analysys
            % function [ motifM, aMat, cMat, pMat, dMat, odMat ] = wsbm_comm_motif(cij,ca)
            [motifAna.(comTypes{cdx}).comMotifMat{idx}(:,:,jdx),...
             motifAna.(comTypes{cdx}).motifEdgeMats{idx}.aMat(:,:,jdx),...
             motifAna.(comTypes{cdx}).motifEdgeMats{idx}.cMat(:,:,jdx),...
             motifAna.(comTypes{cdx}).motifEdgeMats{idx}.pMat(:,:,jdx),...
             motifAna.(comTypes{cdx}).motifEdgeMats{idx}.dMat(:,:,jdx),...
             motifAna.(comTypes{cdx}).motifEdgeMats{idx}.odMat(:,:,jdx),...
            ] = wsbm_comm_motif(CIJ,comInds) ;

        end 

        % get the node-wise rate of com motif participation
        for jdx = 1:size(baseRes.(comTypes{cdx}).ca_K{idx},2)

            tmpMat = motifAna.(comTypes{cdx}).comMotifMat{idx}(:,:,jdx) ;

            % for 1:4, assortative, cor, periph, dis
            for commMotifInd = 1:length(motifVarNAmes)
                 % map the motif participation to the nodes (k = idx+1)
                 k = (idx+1) ;
                 for kdx = 1:k

                    % select the nodes in the k community
                    tmpCommSelect = (kdx == baseRes.(comTypes{cdx}).ca_K{idx}(:,jdx)) ;

                    % compute proportion that community participates in this
                    % specific community motif
                    tmpMotifVal = ( sum(tmpMat(kdx,:) == commMotifInd) + ...
                        sum(tmpMat(:,kdx) == commMotifInd) ) / ...
                        ((k-1)*2) ;

                    % get proportion by dividing by total off-diag blocked consider
                    tmpName = motifVarNAmes{commMotifInd} ;
                    motifAna.(comTypes{cdx}).motifStructs{idx}.(tmpName)(tmpCommSelect,jdx) = ...
                        tmpMotifVal ;

                 end % kdx
            end % commMotifInd
        end % loop number of reps
    end % 1:numLevels

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do a randomized motif analysis

rng(123)

iComm = { 'mod' 'wsbm' } ;
comTypes = { 'randblM' 'randblW' } ;
motifVarNAmes = { 'assort' 'core' 'peri' 'disort' } ; 
CIJ = baseRes.rawData ;

% allocate struct
%motifAna = struct() ; 

% numLevels = size(baseRes.wsbm.ca_K,1) ;

for cdx = 1:length(comTypes)

    % allocate cell array
    motifAna.(comTypes{cdx}).blMatStacks = cell([numLevels 1]) ;
    motifAna.(comTypes{cdx}).comMotifMat = cell([numLevels 1]) ;
    motifAna.(comTypes{cdx}).motifStructs = cell([numLevels 1]) ;
    motifAna.(comTypes{cdx}).motifEdgeMats = cell([numLevels 1]) ;

    % loop across all the kResults for this hemi
    for idx = 1:numLevels

        disp([ 'idx: ' int2str(idx) ])

        numReps = size(baseRes.(iComm{cdx}).ca_K{idx},1) ;

        % allocate kxk matrix at this idx (k = idx+1)
        motifAna.(comTypes{cdx}).blMatStacks{idx} =  ...
            zeros(idx+1,idx+1,numReps) ;

        % allocate kxk com motif matrix
        motifAna.(comTypes{cdx}).comMotifMat{idx} =  ...
            zeros(idx+1,idx+1,numReps) ;

        % allocate com motif struct
        motifAna.(comTypes{cdx}).motifStructs{idx} =  ...
             struct('assort',zeros([NUM_NODES numReps]), ...
                    'core',zeros([NUM_NODES numReps]),...
                    'peri',zeros([NUM_NODES numReps]),...
                    'disort',zeros([NUM_NODES numReps])...
                    ) ;
        tmpName = fieldnames(motifAna.(comTypes{cdx}).motifStructs{idx}) ;
        motifAna.(comTypes{cdx}).motifStructs{idx}.varNames = tmpName ;

        % allocate the com motif edge-wise mats
        motifAna.(comTypes{cdx}).motifEdgeMats{idx}.aMat = zeros([NUM_NODES NUM_NODES numReps]) ;
        motifAna.(comTypes{cdx}).motifEdgeMats{idx}.cMat = zeros([NUM_NODES NUM_NODES numReps]) ;
        motifAna.(comTypes{cdx}).motifEdgeMats{idx}.pMat = zeros([NUM_NODES NUM_NODES numReps]) ;
        motifAna.(comTypes{cdx}).motifEdgeMats{idx}.dMat = zeros([NUM_NODES NUM_NODES numReps]) ;
        motifAna.(comTypes{cdx}).motifEdgeMats{idx}.odMat = zeros([NUM_NODES NUM_NODES numReps]) ;

        % get the stack
        for jdx = 1:size(baseRes.(iComm{cdx}).ca_K{idx},2)
                        
            comInds = baseRes.(iComm{cdx}).ca_K{idx}(:,jdx) ;
            [~,tmpBl] = get_block_mat(CIJ,comInds) ;
            % get the avg block matrix
            motifAna.(comTypes{cdx}).blMatStacks{idx}(:,:,jdx) = tmpBl ;

            % the motif analysys
            % function [ motifM, aMat, cMat, pMat, dMat, odMat ] = wsbm_comm_motif(cij,ca)
            [motifAna.(comTypes{cdx}).comMotifMat{idx}(:,:,jdx),...
             motifAna.(comTypes{cdx}).motifEdgeMats{idx}.aMat(:,:,jdx),...
             motifAna.(comTypes{cdx}).motifEdgeMats{idx}.cMat(:,:,jdx),...
             motifAna.(comTypes{cdx}).motifEdgeMats{idx}.pMat(:,:,jdx),...
             motifAna.(comTypes{cdx}).motifEdgeMats{idx}.dMat(:,:,jdx),...
             motifAna.(comTypes{cdx}).motifEdgeMats{idx}.odMat(:,:,jdx),...
            ] = wsbm_comm_motif(CIJ,comInds,1) ;

        end 

        % get the node-wise rate of com motif participation
        for jdx = 1:size(baseRes.(iComm{cdx}).ca_K{idx},2)

            tmpMat = motifAna.(comTypes{cdx}).comMotifMat{idx}(:,:,jdx) ;

            % for 1:4, assortative, cor, periph, dis
            for commMotifInd = 1:length(motifVarNAmes)
                 % map the motif participation to the nodes (k = idx+1)
                 k = (idx+1) ;
                 for kdx = 1:k

                    % select the nodes in the k community
                    tmpCommSelect = (kdx == baseRes.(iComm{cdx}).ca_K{idx}(:,jdx)) ;

                    % compute proportion that community participates in this
                    % specific community motif
                    tmpMotifVal = ( sum(tmpMat(kdx,:) == commMotifInd) + ...
                        sum(tmpMat(:,kdx) == commMotifInd) ) / ...
                        ((k-1)*2) ;

                    % get proportion by dividing by total off-diag blocked consider
                    tmpName = motifVarNAmes{commMotifInd} ;
                    motifAna.(comTypes{cdx}).motifStructs{idx}.(tmpName)(tmpCommSelect,jdx) = ...
                        tmpMotifVal ;

                 end % kdx
            end % commMotifInd
        end % loop number of reps
    end % 1:numLevels

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save it

saveName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_altMod_motifAna.mat' ] ;
save(saveName,'CIJ','motifAna','-v7.3') ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% entropy across K

commNames = { 'mod' 'randM' 'randblM' 'wsbm' 'randW' 'randblW' } ;

entr_K = struct() ;

for cn = 1:length(commNames)

    entr_K.(commNames{cn}).ent = cell(numLevels,1) ;
    entr_K.(commNames{cn}).sum = cell(numLevels,1) ;
    
    for idx = 1:numLevels
        
            disp(idx)
        
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
%%

dat = baseRes.rawData ;
dat(isnan(dat)) = 0 ;

nanE = isnan(baseRes.rawData) ;

[~,~,bd] = degrees_dir(dat) ;
[~,~,wd] = strengths_dir(dat) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save it

saveName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_altMod_motifEntropy.mat' ] ;
save(saveName,'entr_K','-v7.3') ;


