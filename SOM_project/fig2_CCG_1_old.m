% from Duda 
% this scropt plots normalized and avareged CCG for pv som and pyramidal

figure(10)
intw=0.05;inth=0.05;
Sub_handle=set_subplots(2,2,intw,inth);


clear CCR nCCR LCCR UCCR PairOfCells UniqueCellPairs MeanH0 SDH0

load C:\Dropbox\KepecsLab\_Duda\Analyses\CCG2\pv_pv_H0;
ifinvert = [0 0 0 0 0 0 0 0 0];
for iC = 1:size(CCR,1),
    if ifinvert(iC) == 1,
        nCCR(iC,:) = (fliplr(CCR(iC,:)) - MeanH0(iC)) / SDH0(iC);%        (nanmean(CCR(iC,:))));

    else
        nCCR(iC,:) =  (CCR(iC,:) - MeanH0(iC)) / SDH0(iC);

    end
end

%pv-pv
figure(1)
% time_series = [-30:30];
time_series = [0:60];
SECCR=std(nCCR)/sqrt(size(nCCR,1));

errorshade(time_series,nanmean(nCCR),SECCR,'LineColor',[1 0 0]) %'ShadeColor',[1 0 0]
xlim=[0 60];xtick=[0:10:60];xticklabel=[-30:10:30];ylim=[-8 8]; ytick=[-8:4:8];yticklabel=[-8:4:8];
set(Sub_handle(2),'xlim',xlim,'xtick',xtick,'xticklabel',[],'ylim',ylim,'ytick',ytick,'yticklabel',yticklabel,'TickDir','out','fontsize',12,'fontname','Arial')
Ch1=allchild(gca);
copyobj(Ch1,Sub_handle(2))

for k = 1:size(nCCR,1)
    if any(CCR(k,:)<LCCR(k,:))
        figure(k);plot(CCR(k,:),'k')
        hold on
        plot(LCCR(k,:),'r')
        plot(UCCR(k,:),'r')
    end
end

clear CCR nCCR LCCR UCCR PairOfCells UniqueCellPairs MeanH0 SDH0

load C:\Dropbox\KepecsLab\_Duda\Analyses\CCG2\pv_nontagged_H0;

for iC = 1:size(CCR,1),
    nCCR(iC,:) = (CCR(iC,:) - MeanH0(iC)) / SDH0(iC); %        (nanmean(CCR(iC,:))));

end

 KeepFigure=zeros(1,size(CCR,1));

% pv to nontagged
figure(2)
% time_series = [-30:30];
time_series = [0:60];
xlim=[0 60];xtick=[0:10:60];xticklabel=[-30:10:30];ylim=[-1 4]; ytick=[-1:1:4];yticklabel=[-1:1:4];
SECCR=std(nCCR)/sqrt(size(nCCR,1));

errorshade(time_series,nanmean(nCCR),SECCR,'LineColor',[1 0 0]) %'ShadeColor',[1 0 0]
set(Sub_handle(3),'xlim',xlim,'xtick',xtick,'xticklabel',[-30:10:30],'ylim',ylim,'ytick',ytick,'yticklabel',yticklabel,'TickDir','out','fontsize',12,'fontname','Arial')
Ch1=allchild(gca);
copyobj(Ch1,Sub_handle(3))

clear nCCR CCR LCCR UCCR PairOfCells UniqueCellPairs MeanH0 SDH0


%som-som pairs
load('C:\Dropbox\KepecsLab\_Duda\Analyses\CCG2\som_som_H0.mat')
ifinvert = zeros(1,15);ifinvert(1,7)=1;
for iC =  1:size(CCR,1), 

        nCCR(iC,:) = (CCR(iC,:) - MeanH0(iC)) / SDH0(iC);%        (nanmean(CCR(iC,:))));

end
 nCCR=nCCR(any(nCCR,2),:);


figure(3)

time_series = [0:60];
xlim=[0 60];xtick=[0:10:60];xticklabel=[-30:10:30];ylim=[-8 8]; ytick=[-8:4:8];yticklabel=[-8:4:8];
SECCR=std(nCCR)/sqrt(size(nCCR,1));

errorshade(time_series,nanmean(nCCR),SECCR,'LineColor',[0 0 1]) %'ShadeColor',[1 0 0]
set(Sub_handle(2),'xlim',xlim,'xtick',xtick,'xticklabel',[],'ylim',ylim,'ytick',ytick,'yticklabel',[],'TickDir','out','fontsize',12,'fontname','Arial')
Ch1=allchild(gca);
copyobj(Ch1,Sub_handle(2))
',ylim,'ytick',ytick,'yticklabel',yticklabel,'fontsize',8)
clear CCR nCCR LCCR UCCR PairOfCells UniqueCellPairs MeanH0 SDH0

load('C:\Users\rattmanuser\Dropbox\KepecsLab\_Duda\Analyses\CCG2\som_nontagged_H0.mat')

for iC = 1:size(CCR,1),
    nCCR(iC,:) = (CCR(iC,:) - MeanH0(iC)) / SDH0(iC);%        (nanmean(CCR(iC,:))));

end
    

figure(4)
% time_series = [-30:30];
time_series = [0:60];
xlim=[0 60];xtick=[0:10:60];xticklabel=[-30:10:30];ylim=[-1 4]; ytick=[-1:1:4];yticklabel=[-1:1:4];
SECCR=std(nCCR)/sqrt(size(nCCR,1));
% errorbar(nanmean(nCCR),SECCR,'color','b')
errorshade(time_series,nanmean(nCCR),SECCR,'LineColor',[0 0 1]) %'ShadeColor',[1 0 0]
set(Sub_handle(4),'xlim',xlim,'xtick',xtick,'xticklabel',xticklabel,'ylim',ylim,'ytick',ytick,'yticklabel',[],'TickDir','out','fontsize',12,'fontname','Arial')
Ch1=allchild(gca);
copyobj(Ch1,Sub_handle(4))

% F10H=figure(10);
% set( F10H,'renderer','painter');
