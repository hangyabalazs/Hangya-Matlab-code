%% transmission success rate - complementarity correlation

% f:\balazs\_analysis\Czurko\discriminated2\newcomp.mat should be loaded!

figure
plot(pos(:,2),pos(:,1),'r.','MarkerSize',30)
hold on
plot(neg(:,2),neg(:,1),'b.','MarkerSize',30)

x = [neg(:,2); pos(:,2)];
y = [neg(:,1); pos(:,1)];
X = [ones(size(x)) x];
[b,bint,r,rint,stats] = regress(y,X);
corrcoef(x,y)
R = sqrt(stats(1))         % correlation coefficient (R-value of the regression)
F = stats(2);           % F-test for H0: all coeff.-s are zero
p = stats(3)           % F-test significance

[gr icp err] = linefit(x,y);
xx = (min(x)-0.1):0.01:(max(x)+0.1);
yy = xx .* gr + icp;
plot(xx,yy,'LineWidth',2,'Color','black')
% set(gca,'LineWidth',2,'XTick',[0 0.2 0.4 0.6 0.8 1],'YTick',[0 0.14],...
%     'FontSize',16,'FontWeight','bold','YLim',[0 0.14])
box off

%% transmission success rate - spatial corr. correlation

figure
plot(pos(:,4),pos(:,1),'r.','MarkerSize',30)
hold on
plot(neg(:,4),neg(:,1),'b.','MarkerSize',30)

x = [neg(:,4); pos(:,4)];
y = [neg(:,1); pos(:,1)];
X = [ones(size(x)) x];
[b,bint,r,rint,stats] = regress(y,X);
corrcoef(x,y)
R = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
F = stats(2);           % F-test for H0: all coeff.-s are zero
p = stats(3)           % F-test significance

[gr icp err] = linefit(x,y);
xx = (min(x)-0.1):0.01:(max(x)+0.1);
yy = xx .* gr + icp;
plot(xx,yy,'LineWidth',2,'Color','black')
% set(gca,'LineWidth',2,'XTick',[-0.4 -0.2 0 0.2 0.4 0.6],'YTick',[0 0.14],...
%     'FontSize',16,'FontWeight','bold','YLim',[0 0.14],'XLim',[-0.5 0.75])
box off