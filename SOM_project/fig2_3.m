% from Duda

%  function [bar_h1 bar_h2 bar_h3 bar_h4 bar_h5 bar_h6]=fig2c(non_tagged_pv,non_tagged_som)
figure(3)
xlim=[-30 30];xtick=[-30:10:30];xticklabel=[-30:10:30]; ytick=[0.5:0.1:1.5];yticklabel=[0.5:0.1:1.5];
clear CCR nCCR LCCR UCCR PairOfCells UniqueCellPairs MeanH0 SDH0

load('C:\My Dropbox\KepecsLab\_Duda\Analyses\newCCG\pv_pv_H0.mat')
load('C:\My Dropbox\KepecsLab\_Duda\Analyses\newCCG\pv_pv2.mat')
ifinvert = [1 0 0 0 0 0 1 1 0];
for iC = 1:size(CCR,1),
    if ifinvert(iC) == 1,
        nCCR(iC,:) = (fliplr(CCR(iC,:)) - MeanH0(iC)) / SDH0(iC);%        (nanmean(CCR(iC,:))));
%         nCCR(iC,:) = fliplr(CCR(iC,:));
%         sCCR(iC,:)=filtfilt(b,a,nCCR(iC,:));
%         nCCR(iC,:) = sCCR(iC,:)./(nanmean(sCCR(iC,:)));
    else
        nCCR(iC,:) =  (CCR(iC,:) - MeanH0(iC)) / SDH0(iC);
%         nCCR(iC,:) = CCR(iC,:);
%         sCCR(iC,:)=filtfilt(b,a,nCCR(iC,:));
%         nCCR(iC,:) = sCCR(iC,:)./(nanmean(sCCR(iC,:)));
    end
end
% for iC = 1:size(CCR,1),
%    figure(iC)
%     cla
%     st_h = bar([-30:30],CCR(iC,:),1,'FaceColor','r','EdgeColor','none');
%     set(gca,'TickDir','out','xtick',[],'ytick',[], 'xlim',[-30 30],'box','off')
%     hold on
%     stairs([-30:30],LCCR(iC,:),'b')
%     stairs([-30:30],UCCR(iC,:),'b')
% %     if max(CCR(iC,:))>max(UCCR(iC,:)) | min(CCR(iC,:))<min(LCCR(iC,:))
% %        KeepFigure(1,iC)=1;
% %      end
% end


%pv-pv
subplot(221)
time_series = [-30:30];
SECCR=std(nCCR)/sqrt(size(nCCR,1));
errorbar(nanmean(nCCR),SECCR)
ylim([-7 7.5]);
% stairs(time_series,nanmean(nCCR([1 2 4:9],:)),'color','r');
% set(sub_h(1),'box','off','xlim',xlim,'xtick',xtick,'TickDir','out','xticklabel',[],'ylim',ylim,'ytick',ytick,'yticklabel',yticklabel,'fontsize',8)
clear CCR nCCR LCCR UCCR PairOfCells UniqueCellPairs MeanH0 SDH0

load('C:\My Dropbox\KepecsLab\_Duda\Analyses\newCCG\pv_nontagged_H0.mat')
load('C:\My Dropbox\KepecsLab\_Duda\Analyses\newCCG\pv_nontagged2.mat')


for iC = 1:size(CCR,1),
    nCCR(iC,:) = (CCR(iC,:) - MeanH0(iC)) / SDH0(iC); %        (nanmean(CCR(iC,:))));
%     nCCR(iC,:)=CCR(iC,:);
%     sCCR(iC,:)=filtfilt(b,a,nCCR(iC,:));
%     nCCR(iC,:) = sCCR(iC,:)./(nanmean(sCCR(iC,:)));
end

 KeepFigure=zeros(1,size(CCR,1));
% for iC = 219: 219 %size(CCR,1),
%    figure(iC)
%     cla
%     st_h = bar([-30:30],CCR(iC,:),1,'FaceColor','r','EdgeColor','none');
%     hold on
%     set(gca,'TickDir','out','xtick',[],'ytick',[], 'xlim',[-30 30],'box','off')
%     stairs([-30:30],LCCR(iC,:),'b')
%     stairs([-30:30],UCCR(iC,:),'b')
% %     if max(CCR(iC,:))>max(UCCR(iC,:)) | min(CCR(iC,:))<min(LCCR(iC,:))
% %        KeepFigure(1,iC)=1;
% %      end
% end
% ind=find(KeepFigure==1)
%     axis tight
%     ylim([0.5 1.2])
% end

subplot(222)
time_series = [-30:30];
SECCR=std(nCCR)/sqrt(size(nCCR,1));
errorbar(nanmean(nCCR),SECCR);
ylim([-2 4]);
%  stairs(time_series,nanmean(nCCR),'color','r');
%plot(time_series,nanmean(nCCR),'color','r');
% set(sub_h(2),'box','off','xlim',xlim,'xtick',xtick,'TickDir','out','xticklabel',[],'ylim',ylim,'ytick',ytick,'yticklabel',[],'fontsize',8)
clear nCCR CCR LCCR UCCR PairOfCells UniqueCellPairs MeanH0 SDH0


%som-som pairs
load('C:\My Dropbox\KepecsLab\_Duda\Analyses\newCCG\som_som_H0.mat')
load('C:\My Dropbox\KepecsLab\_Duda\Analyses\newCCG\som_som2.mat')
ifinvert = zeros(1,15);ifinvert(1,7)=1;
for iC =  [3,4,6,7,8,10,11,12,15];size(CCR,1), 
%     if ifinvert(iC) == 1,
%         nCCR(iC,:) = fliplr(CCR(iC,:)./(nanmean(CCR(iC,:))));
%     else
        nCCR(iC,:) = (CCR(iC,:) - MeanH0(iC)) / SDH0(iC);%        (nanmean(CCR(iC,:))));
%         nCCR(iC,:) = CCR(iC,:);
%         sCCR(iC,:)=filtfilt(b,a,nCCR(iC,:));
%         nCCR(iC,:) = sCCR(iC,:)./(nanmean(sCCR(iC,:)));
        %     end
end
 nCCR=nCCR(any(nCCR,2),:);
% for iC =[3,4,6,7,8,10,11,12,15]; % 1:size(CCR,1),
%     figure(iC+100)
%     cla
%     st_h = bar([-30:30],CCR(iC,:),1,'FaceColor','b','EdgeColor','none');
%     set(gca,'TickDir','out','xtick',[],'ytick',[], 'xlim',[-30 30],'box','off')
%     hold on
%     stairs([-30:30],LCCR(iC,:),'r')
%     stairs([-30:30],UCCR(iC,:),'r')
%     
% %     axis tight
% %     ylim([0.5 1.2])
% end

subplot(223)
time_series = [-30:30];
SECCR=std(nCCR)/sqrt(size(nCCR,1));
errorbar(nanmean(nCCR),SECCR);
ylim([-7 7.5]);
% stairs(time_series,nanmean(nCCR),'color','b');
% set(sub_h(4),'box','off','xlim',xlim,'xtick',xtick,'TickDir','out','xticklabel',xticklabel,'ylim',ylim,'ytick',ytick,'yticklabel',yticklabel,'fontsize',8)
clear CCR nCCR LCCR UCCR PairOfCells UniqueCellPairs MeanH0 SDH0

load('C:\My Dropbox\KepecsLab\_Duda\Analyses\newCCG\som_nontagged_H0.mat')
load('C:\My Dropbox\KepecsLab\_Duda\Analyses\newCCG\som_nontagged2.mat')
for iC = 1:size(CCR,1),
    nCCR(iC,:) = (CCR(iC,:) - MeanH0(iC)) / SDH0(iC);%        (nanmean(CCR(iC,:))));
%     nCCR(iC,:) = CCR(iC,:);
%     sCCR(iC,:)=filtfilt(b,a,nCCR(iC,:));
%     nCCR(iC,:) = sCCR(iC,:)./(nanmean(sCCR(iC,:)));
end
    
% for iC =201:size(CCR,1) 
%     figure(iC)
%     cla
%     st_h = bar([-30:30],CCR(iC,:),1,'FaceColor','b','EdgeColor','none');
%     hold on
%     stairs([-30:30],LCCR(iC,:),'r')
%     stairs([-30:30],UCCR(iC,:),'r')
%     set(gca,'TickDir','out','xtick',[],'ytick',[], 'xlim',[-30 30],'box','off')
%     
% %     axis tight
% %     ylim([0.5 1.2])
% end

subplot(224)
time_series = [-30:30];
SECCR=std(nCCR)/sqrt(size(nCCR,1));
errorbar(nanmean(nCCR),SECCR);
ylim([-2 4]);