% from Duda

%  function [bar_h1 bar_h2 bar_h3 bar_h4 bar_h5 bar_h6]=fig2c(non_tagged_pv,non_tagged_som)
figure(3)
[b,a]=butter(10,2/5);
sub_h = set_subplots(2,3,0.04,0.04);
xlim=[-30 30];xtick=[-30:10:30];xticklabel=[-30:10:30];ylim=[0.5 1.5]; ytick=[0.5:0.1:1.5];yticklabel=[0.5:0.1:1.5];
load('pv_pv2.mat')
ifinvert = [1 0 0 0 0 0 1 1 0];
for iC = 1:size(CCR,1),
    if ifinvert(iC) == 1,
        nCCR(iC,:) = fliplr(CCR(iC,:)./(nanmean(CCR(iC,:))));%        (nanmean(CCR(iC,:))));
%         nCCR(iC,:) = fliplr(CCR(iC,:));
%         sCCR(iC,:)=filtfilt(b,a,nCCR(iC,:));
%         nCCR(iC,:) = sCCR(iC,:)./(nanmean(sCCR(iC,:)));
    else
        nCCR(iC,:) =  CCR(iC,:)./(nanmean(CCR(iC,:)));
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
subplot(sub_h(1))
time_series = [-30:30];
SECCR=std(nCCR([1 2 4:9],:))/sqrt(size(nCCR,1));
errorbar(nanmean(nCCR([1 2 4:9],:)),SECCR)
% stairs(time_series,nanmean(nCCR([1 2 4:9],:)),'color','r');
% set(sub_h(1),'box','off','xlim',xlim,'xtick',xtick,'TickDir','out','xticklabel',[],'ylim',ylim,'ytick',ytick,'yticklabel',yticklabel,'fontsize',8)
clear CCR nCCR LCCR UCCR PairOfCells UniqueCellPairs

load('pv_nontagged2.mat')


for iC = 1:size(CCR,1),
    nCCR(iC,:) = CCR(iC,:)./(nanmean(CCR(iC,:))); %        (nanmean(CCR(iC,:))));
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

subplot(sub_h(2))
time_series = [-30:30];
SECCR=std(nCCR)/sqrt(size(nCCR,1));
errorbar(nanmean(nCCR),SECCR);
%  stairs(time_series,nanmean(nCCR),'color','r');
%plot(time_series,nanmean(nCCR),'color','r');
% set(sub_h(2),'box','off','xlim',xlim,'xtick',xtick,'TickDir','out','xticklabel',[],'ylim',ylim,'ytick',ytick,'yticklabel',[],'fontsize',8)
clear nCCR CCR LCCR UCCR PairOfCells UniqueCellPairs


%som-som pairs
load('som_som2.mat')
ifinvert = zeros(1,15);ifinvert(1,7)=1;
for iC =  [3,4,6,7,8,10,11,12,15];size(CCR,1), 
%     if ifinvert(iC) == 1,
%         nCCR(iC,:) = fliplr(CCR(iC,:)./(nanmean(CCR(iC,:))));
%     else
        nCCR(iC,:) = CCR(iC,:)./(nanmean(CCR(iC,:)));%        (nanmean(CCR(iC,:))));
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

subplot(sub_h(4))
time_series = [-30:30];
SECCR=std(nCCR)/sqrt(size(nCCR,1));
errorbar(nanmean(nCCR),SECCR);
% stairs(time_series,nanmean(nCCR),'color','b');
% set(sub_h(4),'box','off','xlim',xlim,'xtick',xtick,'TickDir','out','xticklabel',xticklabel,'ylim',ylim,'ytick',ytick,'yticklabel',yticklabel,'fontsize',8)
clear CCR nCCR LCCR UCCR PairOfCells UniqueCellPairs

load('som_nontagged2.mat')
for iC = 1:size(CCR,1),
    nCCR(iC,:) = CCR(iC,:)./(nanmean(CCR(iC,:)));%        (nanmean(CCR(iC,:))));
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

subplot(sub_h(5))
time_series = [-30:30];
SECCR=std(nCCR)/sqrt(size(nCCR,1));
errorbar(nanmean(nCCR),SECCR);
% stairs(time_series,nanmean(nCCR),'color','b');
% set(sub_h(5),'box','off','xlim',xlim,'xtick',xtick,'TickDir','out','xticklabel',xticklabel,'ylim',ylim,'ytick',ytick,'yticklabel',[],'fontsize',8)
clear CCR nCCR LCCR UCCR PairOfCells UniqueCellPairs
CCR_perc=zeros(4,4);
% for i=1:size(CCR_stat,1)
%     for k=1:size(CCR_stat,2)
%         CCR_perc(i,k)=(CCR(i,k)*100)/sum(CCR_stat(i,:))
%     end
% end
subplot(sub_h(3))
% plot the inhibition duration for Pyramydal cells from Pv vs Som animals
 allcells=listtag('cells');
 [cellids, cellids_ind]=intersect(allcells,non_tagged_pv)
 index= cellids_ind(find(TheMatrix(cellids_ind,10)<0.05 & TheMatrix(cellids_ind,13)>0));
 non_tagged_inh_pv=TheMatrix(index,13);
 N1= hist(non_tagged_inh_pv,[0:5:100]); %stairs([0:5:100],N,'color','r');

   % plotting cum distribution 
%  figure(10)
%  stairs(cumsum(N)/max(cumsum(N)),'color','r')
 [cellids, cellids_ind]=intersect(allcells,non_tagged_som)
 index= cellids_ind(find(TheMatrix(cellids_ind,10)<0.05 & TheMatrix(cellids_ind,13)>0));
 non_tagged_inh_som=TheMatrix(index,13);
 N2=hist(non_tagged_inh_som,[0:5:100]); %stairs([0:5:100],N,'color','b');
 bar([0:5:100],[N2' N1'],1.3);
 h = findobj(gca,'Type','patch');
 set(h(1),'FaceColor','r');
 set(h(2),'FaceColor','b');
 set(sub_h(3),'box','off','xlim',[-2.5 102.5],'xtick',[0:10:100],'TickDir','out','xticklabel',[],'ylim',[0 80],'ytick',[0:20:80],'yticklabel',[0:20:80],'fontsize',8)

 % plotting cum distribution 
%  figure(10)
%  hold on
%  stairs(cumsum(N)/max(cumsum(N)),'color','b')
%  
 subplot(sub_h(6))
% plot the inhibition duration for  Pv vs Som animals
 allcells=listtag('cells');
 [cellids, cellids_ind]=intersect(allcells,pv_cells_v2)
 index= cellids_ind(find(TheMatrix(cellids_ind,10)<0.05 & TheMatrix(cellids_ind,13)>0));
 inh_pv=TheMatrix(index,13);
 N1=hist(inh_pv,[0:5:100]);%stairs([0:5:100],N,'color','r');
 
%  figure(11)
%  stairs(cumsum(N)/max(cumsum(N)),'color','r')
%  
 [cellids, cellids_ind]=intersect(allcells,som_cells_v1)
 index= cellids_ind(find(TheMatrix(cellids_ind,10)<0.05 & TheMatrix(cellids_ind,13)>0));
 inh_som=TheMatrix(index,13);
 N2=hist(inh_som,[0:5:100]); %stairs([0:5:100],N,'color','b');
 bar([0:5:100],[N2' N1'], 1.3);
 h = findobj(gca,'Type','patch');
 set(h(1),'FaceColor','r');
 set(h(2),'FaceColor','b');
 set(sub_h(6),'box','off','xlim',[0 100],'xtick',[0:10:100],'TickDir','out','xticklabel',[0:10:100],'ylim',[0 10],'ytick',[0:2:10],'yticklabel',[0:2:10],'fontsize',8)

%  figure(11)
%  hold on
%  stairs(cumsum(N)/max(cumsum(N)),'color','b')
%  

 
  set(gcf, 'Renderer', 'painters');
 
 
 