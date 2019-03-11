
clc
clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this is where i load data

addpath(genpath('/home/jfaskowi/JOSHSTUFF/projects/ratbrain'))

config_file='config_template_analyzeGridRuns.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_kResults.mat' ] ;
bothHemi = load(loadName) ;
ls(loadName)

config_file='config_template_analyzeGridRuns_oneHemi.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_kResults.mat' ] ;
oneHemi = load(loadName) ;
ls(loadName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this is where i analyze data

hemis = { 'both' 'sing' } ;
numNodes = { 154 77 } ;

inData = struct() ;
inData.both = bothHemi ;
inData.sing = oneHemi ;

gather = struct() ;
gather.both = struct() ;
gather.sing = struct() ;

for hemi = 1:length(hemis)

    gather.(hemis{hemi}).wsbmComs = cell(size(inData.(hemis{hemi}).kResults)) ;
    gather.(hemis{hemi}).wsbmCon = zeros(numNodes{hemi},numNodes{hemi},length(inData.(hemis{hemi}).kResults)) ;
    
    gather.(hemis{hemi}).versAcrossK = zeros(numNodes{hemi},length(inData.(hemis{hemi}).kResults)) ;

    for idx = 1:length(gather.(hemis{hemi}).wsbmComs)
        tmp = cell2mat(cellfun(@(x)wsbm_community_assign(x), ...
            inData.(hemis{hemi}).kResults{idx},'UniformOutput',0));
        gather.(hemis{hemi}).wsbmComs{idx} = reshape(tmp,numNodes{hemi}, length(inData.(hemis{hemi}).kResults{idx}) ) ;
        gather.(hemis{hemi}).wsbmCon(:,:,idx) = agreement(gather.(hemis{hemi}).wsbmComs{idx}) ./ length(inData.(hemis{hemi}).kResults{idx}) ;

        gather.(hemis{hemi}).versAcrossK(:,idx) = get_nodal_versatility(gather.(hemis{hemi}).wsbmComs{idx});
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% more analysis

res = struct(); 

for hemi = 1:length(hemis)
    
    res.(hemis{hemi}).consensus = nanmean(gather.(hemis{hemi}).wsbmCon,3) ;
    res.(hemis{hemi}).versatility = nanmean(gather.(hemis{hemi}).versAcrossK,2) ;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load the template model from oneHemi

config_file='config_template_analyzeGridRuns_oneHemi.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_consensusRuns.mat' ] ;
oneHemiConsensus = load(loadName) ;

OH_templateModel = oneHemiConsensus.consensusCent ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% modularity

templateAdj = OH_templateModel.Data.Raw_Data ; 
templateAdj(~~isnan(templateAdj)) = 0 ;

[ca_mod,ca_modQ] = gsweep_mod_atK_dir(templateAdj,OH_templateModel.R_Struct.k) ;

%% compare 77 node ca with the runs from the both hemi

ca_tempModel = wsbm_community_assign(OH_templateModel) ;

% make a ca for both hemisphere
caDoub_tempModel = [ ca_tempModel ; ca_tempModel ] ;

% now compare this to the ca for each model of kResults
viToBothHemi = cell(size(gather.both.wsbmComs)) ;
nmiToBothHemi = cell(size(gather.both.wsbmComs)) ;

for idx = 1:length(gather.both.wsbmComs)
   
    tmp = gather.both.wsbmComs{idx} ;
    
    [viToBothHemi{idx},nmiToBothHemi{idx}] = partition_distance(caDoub_tempModel,tmp) ;
        
end

% to look at it
figure 
violinplot(cell2mat(viToBothHemi(1:15))')
ppp = gca ;
xticklabels(2:16)
xlabel('{\itk}')
ylabel('VI Distance')
set(gca,'FontSize',16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load the template model from bothHemi

config_file='config_template_analyzeGridRuns.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

loadName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_consensusRuns.mat' ] ;
bothHemiConsensus = load(loadName) ;

BH_templateModel = bothHemiConsensus.consensusCent ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

ca_tempModel = wsbm_community_assign(BH_templateModel) ;

% make a ca for both hemisphere
caHalf_tempModel = ca_tempModel(1:77) ;

% now compare this to the ca for each model of kResults
viToSingHemi = cell(size(gather.both.wsbmComs)) ;
nmiToSingHemi = cell(size(gather.both.wsbmComs)) ;

for idx = 1:length(gather.both.wsbmComs)
   
    tmp = gather.sing.wsbmComs{idx} ;
    
    [viToSingHemi{idx},nmiToSingHemi{idx}] = partition_distance(caHalf_tempModel,tmp) ;
        
end

% to viz
% figure; violinplot(cell2mat(viToSingHemi(1:12))')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make a figure of the consensus

figure
set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.45, 0.8]);

sp = tight_subplot(4,4) ;

for idx = 2:17
   
    disp(idx)
    axes(sp(idx-1)) ;
    imagesc(gather.sing.wsbmCon(:,:,idx))
    title(sp(idx-1),['{\itk}=' int2str(idx)])
    axis square
    xticks('')
    yticks('')
    caxis([0 1])
end

figure;
imagesc(gather.sing.wsbmCon(:,:,idx))
colorbar





