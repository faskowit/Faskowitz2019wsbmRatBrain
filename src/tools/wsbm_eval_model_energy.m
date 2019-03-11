function [B,E,K,EMD] = wsbm_eval_model_energy(origModel,numSims,randModelParam,symmetric,userEvalFuncs)
% EVAL_GENWSBM_MODEL     Generation and evaluation of synthetic networks
%
%   [B,E,K] = EVALUATE_GENERATIVE_MODEL(A,Atgt,D,m,modeltype,modelvar,params) 
%
%   Generates synthetic networks and evaluates their energy function (see
%   below) using the models described in the study by Betzel et al (2016)
%   in Neuroimage.
%
%   Inputs:
%           wsbmModel,  WSBM model struct output from wsbm fitting
%           numSims,    Number of simulated networks to generate (1000)
%
%   Outputs:
%           B,          n x n x numSims matrix of synthetic networks
%           E,          energy for each synthetic network
%           K,          Kolmogorov-Smirnov statistics for each synthetic
%                       network.
%
%
%   Note: Energy is calculated in exactly the same way as in Betzel et
%   al (2016). There are four components to the energy are KS statistics
%   comparing degree, clustering coefficient, betweenness centrality, and 
%   edge length distributions. Energy is calculated as the maximum across
%   all four statistics.
%
%   Reference: Betzel et al (2016) Neuroimage 124:1054-64.
%
%   Richard Betzel, Indiana University/University of Pennsylvania, 2015
%   Josh Faskowitz edited

if nargin < 2
    numSims = 100 ;
end

if nargin < 3
    randModelParam = 0; 
    disp('will not randomize')
end

if nargin < 4
    symmetric = 1; 
    disp('will use symmetric synth mats')
end

if ~exist('userEvalFuncs','var') || isempty(userEvalFuncs)
    
    % setup func handles for evaluations
    evalFuncs = cell(1,1) ;

    % make sure all outputs are column vecs
    evalFuncs{1} = @(A) sum(A,2) + sum(A,1)';
    evalFuncs{2} = @(A) sum(A ~= 0,2) + sum(A ~= 0,1)';
    evalFuncs{3} = @(A) clustering_coef_wd(A);
    evalFuncs{4} = @(A) betweenness_wei(1 ./ A);
    evalFuncs{5} = @(A) betweenness_bin(A)' ;

else
    evalFuncs = userEvalFuncs ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Atgt = origModel.Data.Raw_Data ;
%replace NaN's with 0
Atgt(~~isnan(Atgt)) = 0 ;
%remove diagonal to be safe
nNodes = size(Atgt,1) ;
Atgt(1:nNodes+1:end)=0; 
% Atgt = single(Atgt~=0); % DONT BINARIZE

% emperical stats
X = cell(length(evalFuncs),1);
for idx = 1:length(X)
 
    X{idx} = evalFuncs{idx}(Atgt) ;
end

% record num stats to eval on
numStats = length(X);

% record K-S, EMD stats
K = zeros(numSims,numStats);
EMD = zeros(numSims,numStats);

% record simulated networks
B = zeros([ nNodes nNodes numSims]);

for idx = 1:numSims
        
    % rand the model?
    if randModelParam == 1
        wsbmModel = wsbm_randomize_model_params(origModel,3);
    else
        wsbmModel = origModel; 
    end
    
    % recover the no-NaN output
    [~,b] = wsbm_synth_adj_gen(wsbmModel,symmetric) ;
    b=single(b);         % DONT binarize!
    
    B(:,:,idx) = b ;
    
    Y = cell(numStats,1);
    for jdx = 1:length(Y)
        Y{jdx} = evalFuncs{jdx}(b) ;
    end
    
    for jdx = 1:numStats
        try
            [K(idx,jdx),EMD(idx,jdx)] = fcn_ks(X{jdx},Y{jdx});
        catch
            K(idx,jdx) = NaN ;
            EMD(idx,jdx) = NaN ;
        end
    end
    
    disp(idx)
end

E = max(K,[],2);

function [kstat,emd] = fcn_ks(x1,x2)
binEdges    =  [-inf ; sort([x1;x2]) ; inf];

binCounts1  =  histc (x1 , binEdges, 1);
binCounts2  =  histc (x2 , binEdges, 1);

sumCounts1  =  cumsum(binCounts1)./sum(binCounts1);
sumCounts2  =  cumsum(binCounts2)./sum(binCounts2);

sampleCDF1  =  sumCounts1(1:end-1);
sampleCDF2  =  sumCounts2(1:end-1);

deltaCDF  =  abs(sampleCDF1 - sampleCDF2);
kstat = max(deltaCDF);
emd = sum(deltaCDF);








