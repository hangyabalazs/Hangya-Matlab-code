%% 

% Load CCGs
load ('C:\Dropbox\KepecsLab\_Sachin\mPFC\CCG analysis\CCG3\pv_pv_H02');

% Calculate normalized CCGs
NumPairs = size(CCR,1);
nCCR = nan(size(CCR));
for iC = 1:NumPairs
    nCCR(iC,:) = (CCR(iC,:) - MeanH0(iC)) / SDH0(iC);
end

% Plot effect size
figure
for k = 1:size(nCCR,1)
    if any(CCR(k,27:35) < LCCR(k,27:35))
        line([1 2],[mean(CCR(k,[1:26 36:61])) min(CCR(k,27:35))],'Color','b')
    end
end

for k = 1:size(nCCR,1)
    if any(CCR(k,27:35)>UCCR(k,27:35))
        line([1 2],[mean(CCR(k,[1:26 36:61])) max(CCR(k,27:35))],'Color','r')
    end
end

%%

clear CCR nCCR LCCR UCCR PairOfCells UniqueCellPairs MeanH0 SDH0

load ('C:\Dropbox\KepecsLab\_Sachin\mPFC\CCG analysis\CCG3\pv_nontagged_H02');

for iC = 1:size(CCR,1),
    nCCR(iC,:) = (CCR(iC,:) - MeanH0(iC)) / SDH0(iC); %        (nanmean(CCR(iC,:))));

end

 KeepFigure=zeros(1,size(CCR,1));

% pv to nontagged
figure(2)
% time_series = [-30:30];
time_series = [0:60];
xlim=[10 50];xtick=[10:10:50];xticklabel=[-20:10:20];ylim=[-1 2]; ytick=[-1:1:2];yticklabel=[-1:1:2];
SECCR=std(nCCR)/sqrt(size(nCCR,1));

errorshade(time_series,nanmean(nCCR),SECCR,'LineColor',[1 0 0]) %'ShadeColor',[1 0 0]
set(Sub_handle(3),'xlim',xlim,'xtick',xtick,'xticklabel',[-30:10:30],'ylim',ylim,'ytick',ytick,'yticklabel',yticklabel,'TickDir','out','fontsize',12,'fontname','Arial')
Ch1=allchild(gca);
copyobj(Ch1,Sub_handle(3))

clear nCCR CCR LCCR UCCR PairOfCells UniqueCellPairs MeanH0 SDH0


%som-som pairs
load('C:\Dropbox\KepecsLab\_Sachin\mPFC\CCG analysis\CCG3\som_som_H02.mat')
ifinvert = zeros(1,15);ifinvert(1,7)=1;
for iC =  1:size(CCR,1), 

        nCCR(iC,:) = (CCR(iC,:) - MeanH0(iC)) / SDH0(iC);%        (nanmean(CCR(iC,:))));

end
 nCCR=nCCR(any(nCCR,2),:);


figure(3)

time_series = [0:60];
xlim=[10 50];xtick=[10:10:50];xticklabel=[-20:10:20];ylim=[-6 9]; ytick=[-6:4:9];yticklabel=[-6:4:9];
SECCR=std(nCCR)/sqrt(size(nCCR,1));

errorshade(time_series,nanmean(nCCR),SECCR,'LineColor',[0 0 1]) %'ShadeColor',[1 0 0]
set(Sub_handle(2),'xlim',xlim,'xtick',xtick,'xticklabel',[],'ylim',ylim,'ytick',ytick,'yticklabel',[],'TickDir','out','fontsize',12,'fontname','Arial')
Ch1=allchild(gca);
copyobj(Ch1,Sub_handle(2))
% ',ylim,'ytick',ytick,'yticklabel',yticklabel,'fontsize',8)
clear CCR nCCR LCCR UCCR PairOfCells UniqueCellPairs MeanH0 SDH0

load('C:\Dropbox\KepecsLab\_Sachin\mPFC\CCG analysis\CCG3\som_nontagged_H02.mat')

for iC = 1:size(CCR,1),
    nCCR(iC,:) = (CCR(iC,:) - MeanH0(iC)) / SDH0(iC);%        (nanmean(CCR(iC,:))));

end
    

figure(4)
% time_series = [-30:30];
time_series = [0:60];
xlim=[10 50];xtick=[10:10:50];xticklabel=[-20:10:20];ylim=[-1 2]; ytick=[-1:1:2];yticklabel=[-1:1:2];
SECCR=std(nCCR)/sqrt(size(nCCR,1));
% errorbar(nanmean(nCCR),SECCR,'color','b')
errorshade(time_series,nanmean(nCCR),SECCR,'LineColor',[0 0 1]) %'ShadeColor',[1 0 0]
set(Sub_handle(4),'xlim',xlim,'xtick',xtick,'xticklabel',xticklabel,'ylim',ylim,'ytick',ytick,'yticklabel',[],'TickDir','out','fontsize',12,'fontname','Arial')
Ch1=allchild(gca);
copyobj(Ch1,Sub_handle(4))

% F10H=figure(10);
% set( F10H,'renderer','painter');
