%% open

uiopen('C:\Balazs\_analysis\SniffWhisk\Fig_meanphase\meanphase_final1.fig',1)

%%

ln = findobj(allchild(gcf),'type','line');
set(ln,'LineWidth',4)
L1 = xlabel('Phase (deg)');
setmyplot(gca,L1)
set(gcf,'renderer','painters')