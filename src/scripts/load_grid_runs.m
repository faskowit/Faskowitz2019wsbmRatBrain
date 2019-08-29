
clc
clearvars

%% load the necessary data

config_file='config_template_rb2_analyzeGridRuns.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

addpath(genpath(pwd))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load

origKResults = read_WSBMcompile_results(GRID_OUTPUT_DIR) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% remove communities where incorrect number of coms was rendered

coms = cell(size(origKResults)) ;
kResults1 = cell(size(origKResults)) ;

% get the coms
for idx = 1:length(COM_NUM_RANGE)
    tmp = cell2mat(cellfun(@(x)wsbm_community_assign(x), ...
        origKResults{idx},'UniformOutput',0));
    tmp2 = reshape(tmp,NUM_NODES, length(origKResults{idx}) ) ;
    
    incVec = ones(length(origKResults{idx}),1) ;
    
    for jdx = 1:length(origKResults{idx})
       
        nnn = length(unique(tmp2(:,jdx))) ;
        
        if nnn ~= COM_NUM_RANGE(idx) 
           disp('wrong number communities')
           incVec(jdx) = 0 ;
        end
    end
    
    % possibly exclude some
    kResults1{idx} = origKResults{idx}(~~incVec) ;
    
end

%% for fairness, randomly select same number of sample models for

rng(123)

kResults2 = cell(size(kResults1)) ;

for idx = 1:length(RESAMP_DIVS)
    for jdx = 1:length(RESAMP_DIVS{idx}.range)
    
        kdx = RESAMP_DIVS{idx}.range(jdx);
        
        tmp = datasample(1:length(kResults1{kdx}),RESAMP_DIVS{idx}.size,'Replace',RESAMP_DIVS{idx}.resamp) ;
        kResults2{kdx} = kResults1{kdx}(tmp) ;
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% lets try to plot the log evidences

% inline func
getLogEvid = @(x) x.Para.LogEvidence ;

numK = size(kResults2,1) ;
actualK = 2:20 ;

figure

for idx = 1:numK

   tmpData = cellfun(getLogEvid,kResults2{idx}) ;
   scatter(actualK(idx) .* ones(length(tmpData),1),tmpData)
   hold on
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% before we save, get rid of data to make the file smaller

rawData = origKResults{1}{1}.Data.Raw_Data ;

for idx = 1:length(kResults2)
   for jdx = 1:length(kResults2{idx})
      
       kResults2{idx}{jdx}.Data = [] ; 
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save the kResults for later usage

kResults = kResults2 ;

saveName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_kResults.mat' ] ;
save(saveName,'kResults','rawData','-v7.3') 
ls(saveName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% finally, make a data with just some basic info
% take less time loading it in future!

baseRes = struct() ;
baseRes.wsbm.ca_K = cell(size(kResults)) ;
baseRes.wsbm.logEvid_K = cell(size(kResults)) ;
baseRes.wsbm.predE_K = cell(size(kResults)) ;
baseRes.wsbm.predW_K = cell(size(kResults)) ;
baseRes.wsbm.exmplModel_K = cell(size(kResults)) ;

getPredE = @(x) x.Para.predict_e ;
getPredW = @(x) x.Para.predict_w ;

baseRes.rawData = rawData ;

for idx = 1:length(baseRes.wsbm.ca_K)
    
   baseRes.wsbm.ca_K{idx} = align_com_labeling(...
       cell2mat(cellfun(@(x)wsbm_community_assign(x), kResults{idx}' ,'UniformOutput',0)) ...
       ) ;
   
   baseRes.wsbm.logEvid_K{idx} = cellfun(getLogEvid,kResults{idx}) ;
   baseRes.wsbm.predE_K{idx} = cell2mat(cellfun(getPredE,kResults{idx},'UniformOutput',0)') ;
   baseRes.wsbm.predW_K{idx} = cell2mat(cellfun(getPredW,kResults{idx},'UniformOutput',0)') ;
   baseRes.wsbm.exmplModel_K{idx} = kResults{idx}{1} ;
   
end

%%

saveName = [ PROJECT_DIR '/data/processed/' OUTPUT_STR '_' GRID_RUN '_baseRes.mat' ] ;
save(saveName,'baseRes','-v7.3') 
ls(saveName)









