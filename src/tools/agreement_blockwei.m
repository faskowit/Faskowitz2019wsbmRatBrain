function D = agreement_blockwei(ci,buffsz,data)
%AGREEMENT      Agreement matrix from clusters
% but weighted differently

nReps = size(ci,2);
k = max(max(ci)) ;
numNodes = size(ci,1) ;

if nargin < 2
    buffsz = 1000;
end

if nReps <= buffsz
    
    ind = dummyvar(ci);
%     D = ind*ind';

    mods = zeros([k nReps]) ;
    betwn = zeros([k nReps]) ;
    for idx=1:nReps
       [tmp] = get_block_mat(data,ci(:,idx)) ;
       mods(:,idx) = tmp(~~eye(k)) ;
       betwn(:,idx) = sum(tmp,2) - mods(:,idx) ;
    end
    
    modRatio = mods ./ betwn ;
    
    ind2 = ind.*modRatio(:)' ;
    
    D = ind*ind2' ;
    
else
    
    a = 1:buffsz:nReps;
    b = buffsz:buffsz:nReps;
    
    if length(a) ~= length(b)
        b = [b, nReps];
    end
    
    x = [a' b'];
    nbuff = size(x,1);
    
    D = zeros(size(ci,1));
    for i = 1:nbuff
       y = ci(:,x(i,1):x(i,2));
       ind = dummyvar(y);
       D = D + ind*ind';
    end
    
end

D = D.*~eye(length(D));
