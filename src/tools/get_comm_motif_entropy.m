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
    
    entMat = -1 .* ( aMat2.*log2or0(aMat2) + ...
            cMat2.*log2or0(cMat2) + ...
            pMat2.*log2or0(pMat2) + ...
            dMat2.*log2or0(dMat2) ) ; 
else
    % check to make sure the sums add to 1
    check_probs(cat(3,aMat,cMat,pMat,dMat,odMat))
    
    entMat = -1 .* ( aMat.*log2or0(aMat) + ...
        cMat.*log2or0(cMat) + ...
        pMat.*log2or0(pMat) + ...
        dMat.*log2or0(dMat) + ...
        odMat.*log2or0(odMat)) ; 
end

function check_probs(stackedMat)
    % check within a tol (because matlab has decimal issues sometimes)
    
    sums = sum(stackedMat,3) ;
    
    % get nonzero
    nz = nonzeros(sums) ;
    ismem = abs(ones(length(nz),1) - nz) < 10e-6 ;

    if sum(ismem) ~= length(ismem)
       error('probabilies not good') 
    end
end

function out = log2or0(in)  
    % function to calculate log2, but if 0, return 0
    out = log2(in) ;
    out(isinf(out)) = 0 ;
end

end % main function