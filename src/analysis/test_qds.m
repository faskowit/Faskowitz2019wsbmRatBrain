%% prep the data

addpath('~/JOSHSTUFF/projects/ratbrain2/src/external/WSBM_v1.2/')
addpath('~/JOSHSTUFF/projects/ratbrain2/data/')

ll = '/home/jfaskowi/JOSHSTUFF/projects/ratbrain2/data/processed/gridRuns1_results_expBern_baseRes.mat' ;
load(ll)
dat = baseRes.rawData ;
dathemi = dat(1:77,1:77) ;

raw_list = Adj2Edg(dathemi) ;
% writeit
writematrix(raw_list,'./ratbrain.csv','Delimiter','\t')

%% code 

import java.util.*;

% Please change the following two variables according to the properties the
% networks
isUnweighted = false;
isUndirected = false;

% The path and filename of the networks
networkFile = 'ratbrain.csv';
% The file to store the detected communities
communityResultFile = 'example_data/result.txt';

[net, totalEdge, totalWeight] = network.getNetwork(networkFile, isUnweighted, isUndirected);

disp(['Number of nodes = ' num2str(net.size) ', number of edges = ' num2str(totalEdge)]);

if ~isUndirected
    revnet=reversedNetwork.getReversedNetwork(networkFile, isUnweighted, isUndirected);
else
    revnet=0;
end

%% loop it to see if its stochastic
% its not

numReps = 5 ; 
COMMS = zeros(net.size,numReps) ; 

% LOOP IT
for idx = 1:5

    % The communities
    % Get communities from file or treat all nodes as a whole community
    communities = ArrayList;
    community = HashSet;
    nodeIter = net.entrySet.iterator;
    while nodeIter.hasNext
        item = nodeIter.next;
        community.add(item.getKey);
    end
    communities.add(community);

    % Fine-tuned Qds
    communities = fineTuneQds_sort(net,revnet,totalEdge,totalWeight,isUndirected,communities);

    % Output community detection result
    communitySize=communities.size;
    for i=0:communitySize-1
       community=communities.get(i);
       disp(community);
       comNodeIter=community.iterator;
       while comNodeIter.hasNext
           nodeId=comNodeIter.next;
           COMMS(nodeId,idx) = i+1 ;
           disp(nodeId)
       end
    end

end

%% visualize it

%   >> load AIJ;                                % load adjacency matrix
%   >> [C,Q] = modularity_louvain_und(AIJ);     % get community assignments
%   >> [X,Y,INDSORT] = fcn_grid_communities(C); % call function
%   >> imagesc(AIJ(INDSORT,INDSORT));           % plot ordered adjacency matrix
%   >> hold on;                                 % hold on to overlay community visualization
%   >> plot(X,Y,'r','linewidth',2);             % plot community boundaries
[X,Y,INDSORT] = grid_communities(COMMS(:,1)) ;
imagesc(dathemi(INDSORT,INDSORT))
hold on
plot(X,Y,'r','linewidth',2)
axis square


