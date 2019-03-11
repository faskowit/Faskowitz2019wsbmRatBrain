function [ kResults ] = read_WSBMcompile_results(inputPath)
% function to read the results of the 

if nargin < 1
   error('need one arg') 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inputPath = [ inputPath '/' ] ;

% first, read the folders at the input path
kDirs = dir(inputPath) ;
kDirs(1:2) = [] ;

% get only the directories
kDirsInd = [ kDirs.isdir ] ;
kDirs = kDirs(kDirsInd) ;

% how many kDirs are there?
numDirs = length(kDirs) ;

% store the results in a cell array
kResults = cell(numDirs,1) ;

for idx = 1:numDirs

    currDir = [ kDirs(idx).folder '/' kDirs(idx).name '/' ] ;
    
    % mats for this dir
    currMatsStruct = dir(currDir) ;
    currMatsStruct(1:2) = [] ;
    
    currMatsNames = arrayfun(@(x) [ x.folder '/' x.name ] , ...
        currMatsStruct,'UniformOutput',0) ;
    
    numMats = length(currMatsStruct) ;
    
    if numMats > 1000
       numMats = 1000 ; 
    end
    
    res = cell(numMats,1) ;
    
    for jdx = 1:numMats
    
        disp(jdx)
        tmp = load(currMatsNames{jdx}) ;
        res{jdx} = tmp.Model ;
        
    end
    
    kResults{idx} = res ;
    
end


