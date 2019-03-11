
clc
clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this is where i load data

config_file='config_template_rb2_oneHemi_analyzeGridRuns.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_baseRes.mat' ] ;
% loads a struct named 'baseRes'
load(loadName) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% gather the data into a struct
 
gather = struct() ;

numLevels = size(baseRes.wsbm.logEvid_K,1) ;

% initialize variables
gather.wsbmCon = zeros(NUM_NODES,NUM_NODES,numLevels) ;
gather.versAcrossK = zeros(NUM_NODES,numLevels) ;

% loop across all the kResults for this hemi
for idx = 1:numLevels

    % already gathered
    gather.wsbmComs{idx} = baseRes.wsbm.ca_K{idx} ;
    gather.logEvid{idx} = baseRes.wsbm.logEvid_K{idx} ;

    numRuns = size(baseRes.wsbm.logEvid_K{idx},1) ;

    % recover the agreement matrix (make into probability by dividing
    % by number of algo. repeitions
    gather.wsbmCon(:,:,idx) = ...
        agreement(gather.wsbmComs{idx}) ./ ...
        size(gather.wsbmComs{idx},2) ;

    % versatility across the community assignments 
    gather.versAcrossK(:,idx) = ...
        get_nodal_versatility(gather.wsbmComs{idx});

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do the motif analysis

comTypes = { 'wsbm' 'mod' } ;
motifVarNAmes = { 'assort' 'core' 'peri' 'disort' } ; 
CIJ = baseRes.rawData ;

% allocate struct
motifAna = struct() ; 

numLevels = size(baseRes.wsbm.ca_K,1) ;

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
%% save it

saveName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_motifAna.mat' ] ;
save(saveName,'CIJ','gather','motifAna','-v7.3') ;
