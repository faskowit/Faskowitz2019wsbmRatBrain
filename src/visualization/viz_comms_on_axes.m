function [] = viz_comms_on_axes(ca,cmap,lineWidth)
% visualize communities on the axis of an imagesc
%
% ca -> community affiliation vector
% cmap -> color map to use, which should be a (num colors)x3 matrix 
% lineWidth -> width of the line to draw

hold on

if nargin < 3
    lineWidth = 8 ;
end

if ~iscolumn(ca)
   ca = ca' ; 
end

ca = sort(ca,'ascend') ;

nc = length(unique(ca)) ;
numNodes = length(ca) ;

for i = 1:nc
    ind = find(ca == i);
    if ~isempty(ind)
        mn = min(ind) - 0.5;
        mx = max(ind) + 0.5;
        x = [ 0.5 0.5 ];
        y = [ mn mx ];
        p = plot(x,y,'Color',cmap(i,:),'linewidth',lineWidth);      
        x = [ numNodes+0.5 numNodes+0.5 ] ;
        pp = plot(y,x,'Color',cmap(i,:),'linewidth',lineWidth);
                
        uistack(pp,'bottom')
        uistack(p,'bottom')
    end
end

hold off


