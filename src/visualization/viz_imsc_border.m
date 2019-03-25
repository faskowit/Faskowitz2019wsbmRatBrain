function [p] = viz_imsc_border(color,lineWidth)
% color borders on imagesc

hold on

if nargin < 2
    lineWidth = 8 ;
end

x_lims = get(gca,'XLim') ;
y_lims = get(gca,'Ylim') ;

x = [ x_lims(1) x_lims(2) x_lims(2) x_lims(1) x_lims(1) x_lims(2) ];
y = [ y_lims(1) y_lims(1) y_lims(2) y_lims(2) y_lims(1) y_lims(1) ] ;
p = plot(x,y,'Color',color,'linewidth',lineWidth);   

uistack(p,'bottom')

hold off


