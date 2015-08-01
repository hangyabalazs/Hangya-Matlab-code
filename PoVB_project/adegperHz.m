%%

x = get(findobj(gca,'type','line'),'XData');
y = get(findobj(gca,'type','line'),'YData');
x = x{2};
y = y{2};

% figure;plot(x,y,'.')

%%

[b,bint,r,rint,stats] = regress(y',[ones(length(x),1) x']);
icp = b(1);
gr = b(2)
xx = [min(x):0.01:max(x)];
yy = xx .* gr + icp;
hold on
plot(xx,yy)

%%

x = get(findobj(gca,'type','line'),'XData');
y = get(findobj(gca,'type','line'),'YData');
x = x{1};
y = y{1};

% figure;plot(x,y,'.')

%%

[b,bint,r,rint,stats] = regress(y',[ones(length(x),1) x']);
icp = b(1);
gr = b(2)
xx = [min(x):0.01:max(x)];
yy = xx .* gr + icp;
hold on
plot(xx,yy)