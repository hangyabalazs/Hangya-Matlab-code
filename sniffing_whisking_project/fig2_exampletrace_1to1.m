%% open

uiopen('C:\Balazs\_analysis\SniffWhisk\Fig_examples\example2_1to1_sniff.fig',1)
ln1k = findobj(allchild(gca),'type','line','color','k');
ln1r = findobj(allchild(gca),'type','line','color','r');

%%

uiopen('C:\Balazs\_analysis\SniffWhisk\Fig_examples\example2_1to1_whisk.fig',1)
ln2 = findobj(allchild(gca),'type','line','color','r');

%%

figure
A = set_subplots(2,1,0,0);
A1 = A(1);
A2 = A(2);
linkaxes([A1 A2],'x')
axis([A1 A2],'off')
xlim([390 1300])

% A1 = subplot(211);
lnn1k = copyobj(ln1k,A1);
set(lnn1k,'Color',[255 204 0]/255,'LineWidth',4)
lnn1r = copyobj(ln1r,A1);
set(lnn1r,'Color','k','LineWidth',2)

% A2 = subplot(212);
lnn2 = copyobj(ln2,A2);
set(lnn2,'Color',[0 153 0]/255,'LineWidth',2)  % ,'Color',[153 0 255]/255)      % blue: [0 153 255]/255, green: [0 153 0]/255

x_lim = xlim;
y_lim = ylim;
line([x_lim(1) x_lim(1)+100],[y_lim(1) y_lim(1)],'Color','k','LineWidth',8)