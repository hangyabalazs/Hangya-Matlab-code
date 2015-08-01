%% open

uiopen('C:\Balazs\_analysis\SniffWhisk\Fig_meancoherence\allcoh_sel_segments_new.fig',1)

%%

ln = findobj(allchild(gcf),'type','line');

cmp = bone(length(ln)+1);
cmp = cmp(1:end-1,:);

for k = 1:length(ln)
    set(ln(k),'Color',cmp(k,:),'LineWidth',3)
end

L(1) = xlabel('Frequency (Hz)');
L(2) = ylabel('Coherence');

setmyplot(gca,L)
xlim([0 25])
axis square