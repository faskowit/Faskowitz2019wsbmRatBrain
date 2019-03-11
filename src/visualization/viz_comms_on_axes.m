function [] = viz_comms_on_axes(CA,cmap,lineWidth)

if nargin < 3
    lineWidth = 8 ;
end

if ~iscolumn(CA)
   CA = CA' ; 
end

CA = sort(CA,'ascend') ;

nc = length(unique(CA)) ;
numNodes = length(CA) ;

for i = 1:nc
    ind = find(CA == i);
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

