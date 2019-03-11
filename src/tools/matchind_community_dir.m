function [ca] = matchind_community_dir(dat,numK,reps)

if nargin < 3
   reps = 1000 ;
end

% first convert the data to matching index
[~,~,matchind] = matching_ind(dat) ;

matchind = matchind + matchind' ;

% get communities
matchModCA = sweep_gamma_louvain_atK(matchind,numK,[],reps,1) ;

% return the central model
ca = wsbm_cent_mod(matchModCA) ;









