function [ pl ] = viz_FnE_Regression(coefs)

hold on;

rr = range(coefs.x) ;

% get the confidence intervals
xVec = (min(coefs.x)*0.99 )...
    :(rr/200):...
    (max(coefs.x)*1.01);
bootCI = zeros(size(coefs.boot,1),length(xVec)) ;

% bootstrap
for idx = 1:size(coefs.boot,1)
    bootCI(idx,:) = polyval(coefs.boot(idx,:),xVec);  
end

% 95%
p95 = prctile(bootCI,[2.5,97.5]);    

% fill that 95%
fill([xVec fliplr(xVec)],[p95(1,:) fliplr(p95(2,:))],[ 0.75 0.75 0.75 ],'facealpha',.7,'edgealpha',0);

% and the line
yhat = polyval(coefs.full,xVec); 
plot(xVec,yhat,'Color',[0.5 0.5 0.5],'linewidth',2.5);











