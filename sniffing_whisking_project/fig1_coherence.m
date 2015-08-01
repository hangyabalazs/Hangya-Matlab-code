%% open

uiopen('C:\Balazs\_analysis\SniffWhisk\Fig_coherence\averaged_pellet_triggered_coherence4.fig',1)

%%

ax = findobj(allchild(gcf),'type','axes');

for k = 1:length(ax)
    colormap('jet')
    axes(ax(k))
    line([10 10],ylim,'Color','w','LineStyle',':','LineWidth',2)
    line([14 14],ylim,'Color','w','LineStyle',':','LineWidth',2)
    L(1) = xlabel('Time - tone onset (s)');
    L(2) = ylabel('Frequency (Hz)');
    setmyplot(ax(k),L)
    
end

set(gcf,'Color','w')