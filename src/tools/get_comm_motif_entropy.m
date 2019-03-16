function [entMat] = get_comm_motif_entropy(aMat,cMat,pMat,dMat,odMat)
% get the comm motif entropy for edge-wise probabilities
% eq. 10, https://doi.org/10.1038/s41467-017-02681-z
% 
% altered to include within-community probability

if nargin < 5
    
    warning(['since odMat not included, rescaling probabilities of aMat'... 
    ' cMat pMat and dMat to sum to 1']) ;
    
    % need to re-normalize, so the probs sum to 1
    tot = aMat + cMat + pMat + dMat ;
    aMat2 = aMat ./ tot ;
    cMat2 = cMat ./ tot ;
    pMat2 = pMat ./ tot ;
    dMat2 = dMat ./ tot ;
    
    entMat = -1 .* ( aMat2.*log2(aMat2) + ...
            cMat2.*log2(cMat2) + ...
            pMat2.*log2(pMat2) + ...
            dMat2.*log2(dMat2) ) ; 
else
    entMat = -1 .* ( aMat.*log2(aMat) + ...
        cMat.*log2(cMat) + ...
        pMat.*log2(pMat) + ...
        dMat.*log2(dMat) + ...
        odMat.*log2(odMat)) ; 
end

