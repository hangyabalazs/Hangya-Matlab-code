%% obsolete

uiopen('C:\Balazs\_analysis\SniffWhisk\Fig_examples\example1_1to1_sniff.fig',1)
seg = getappdata(gcf,'segment');
snwexampleplotter2(seg)

%%

uiopen('C:\Balazs\_analysis\SniffWhisk\Fig_methods\sniffing_recolored.fig',1)

xlim([470 775])
setmyplot(gca)
box off
axs = findobj(allchild(gcf),'type','axes');
axis(axs,'off')
set(gcf,'Color','w')


%% obsolete

uiopen('C:\Balazs\_analysis\SniffWhisk\Fig_examples\example1_1to1_whisk.fig',1)

xlim([470 775])

%%

uiopen('C:\Balazs\_analysis\SniffWhisk\Fig_methods\whisking_recolored.fig',1)

xlim([470 775])
setmyplot(gca)
box off
axis off
set(gcf,'Color','w')
xl = xlim;
line([xl(1) xl(1)+50],[-12 -12],'LineWidth',6,'Color','k')