function [] = viz_labs_on_axes(ca,labs)

hold on

if ~iscolumn(ca)
   ca = ca' ; 
end
  
comms = unique(ca) ;
nc = length(comms) ;

if nargin < 2
    labs = comms ;
end

bb = diff(sort(ca)) ;
breaks = find(bb)' ;

% compute some yticks
breaks2 = [ 0 breaks length(ca) ] ;
midlabelpoint = zeros([nc 1]);
for idx = 1:length(midlabelpoint)
midlabelpoint(idx) = floor( (breaks2(idx+1) - breaks2(idx)) / 2) + breaks2(idx);  
end

set(gca,'ytick',midlabelpoint)
set(gca,'yticklabel',labs)
set(gca,'ticklength',[ 0 0 ]) 

hold off
