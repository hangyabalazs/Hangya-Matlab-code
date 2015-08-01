%%

Pairs = cell(0,2);
for k = 1:184
    Pairs{k,1} = PairOfCells{k}{1}{1};
    Pairs{k,2} = PairOfCells{k}{2}{1};
end

%%

ntinx = [];
for k = 1:152
    pinx1 = strcmp(poc(k,1),Pairs(:,1));
    pinx2 = strcmp(poc(k,2),Pairs(:,2));
    pinx = pinx1 & pinx2;
    ntinx(k) = find(pinx);
end

%%

for iC = ntinx
    nCCR(iC,:) = (CCR(iC,:) - MeanH0(iC)) / SDH0(iC); %        (nanmean(CCR(iC,:))));
end

%%
figure
% time_series = [-30:30];
time_series = [0:60];
xlim=[0 60];xtick=[0:10:60];xticklabel=[-30:10:30];ylim=[-1 4]; ytick=[-1:1:4];yticklabel=[-1:1:4];
SECCR=std(nCCR)/sqrt(size(nCCR,1));

errorshade(time_series,nanmean(nCCR),SECCR,'LineColor',[1 0 0]) %'ShadeColor',[1 0 0]
set(Sub_handle(3),'xlim',xlim,'xtick',xtick,'xticklabel',[-30:10:30],'ylim',ylim,'ytick',ytick,'yticklabel',yticklabel,'TickDir','out','fontsize',12,'fontname','Arial')
