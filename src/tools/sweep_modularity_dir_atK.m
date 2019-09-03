function [ comsAcrossG2, qAcrossG2 ] = sweep_modularity_dir_atK(CIJ,k,gRangeInit,numComsAtK,quiet)
% simple function to return community partitions from louvain, using the 
% functions of BCT, at a specific k

if nargin < 2
   error('need at least two args: CIJ, k') 
end

if ~exist('gRangeInit','var') || isempty(gRangeInit)
    gRangeInit = 0.1:0.005:3 ;
end

if ~exist('numComsAtK','var') || isempty(numComsAtK)
    numComsAtK = 1000 ;
end

if ~exist('quiet','var') || isempty(quiet)
    quiet = false ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initial gamma loop

nNodes = size(CIJ,1);

% first lets sweep gamma
comsAcrossG = zeros([nNodes length(gRangeInit)]);
for idx = 1:length(gRangeInit)
    
    if ~quiet
        disp(idx)
    end
    comsAcrossG(:,idx) = modularity_dir(CIJ,gRangeInit(idx));
end

% get the num coms that had the largests plateau
% get num coms at each sweep iteration 
numComsAcrossG = max(comsAcrossG)' ;
kInd = find(numComsAcrossG == k) ;

if isempty(kInd) || (kInd(1) == kInd(end))
   error('gamma range not big enough') 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% gamma loop for k

% number of finaly com assignments we want
numGReps = numComsAtK ;

startG = gRangeInit(kInd(1)) ;
endG = gRangeInit(kInd(end)) ;

% get gammas uniformally fomr 
gRange2 = rand_in_range(startG,endG,numGReps) ;

% preallocate the output vars
comsAcrossG2 = zeros([nNodes numGReps]);
qAcrossG2 = zeros([numGReps 1]) ;

% runit

idx = 1 ;
lattempt = 1 ;
maxlattempt = 10 ; 
while(idx <= numGReps)

    if lattempt > maxlattempt
        % get a new gamma here
        gRange2(idx) = rand_in_range(startG,endG,1) ;
        lattempt = 1 ;
    end

    [tmp,tmpq] = modularity_dir(CIJ,gRange2(idx));
    if max(tmp) ~= k
       lattempt = lattempt + 1 ;
       continue
    end
    if ~quiet
        disp([ 'generated com:' num2str(idx) ])
    end
    comsAcrossG2(:,idx) = tmp ;
    qAcrossG2(idx) = tmpq ;
    idx = idx + 1 ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% helper functions

% function [M,qmod] = com_louvain_Qmax(adj,gam)
% 
%     % use the found gamma and squeeze any extra modularity outta it
%     % Iterative community finetuning.
%     % W is the input connection matrix.
%     Q0 = -1; Q1 = 0;            % initialize modularity values
%     while Q1-Q0>1e-2            % while modularity increases
%         Q0 = Q1; 
%         [M,qmod] = community_louvain(adj, gam);
%     end
% 
% end % end com_louvain_Qmax

function retsamps = rand_in_range(minVal,maxVal,numSamps)
   
    if maxVal <= minVal
       error('max and min vals wacky') 
    end
   
    % https://www.mathworks.com/help/matlab/math/floating-point-numbers-within-specific-range.html
    retsamps = (maxVal-minVal).*rand(numSamps,1) + minVal;
end


end % end the main function
