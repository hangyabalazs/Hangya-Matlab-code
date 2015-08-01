% from Duda

figure

% pv-pv
clear CCR nCCR LCCR UCCR PairOfCells UniqueCellPairs MeanH0 SDH0
load('C:\Dropbox\KepecsLab\_Sachin\mPFC\CCG analysis\newCCG2\pv_pv_H0.mat')
% load('C:\Dropbox\KepecsLab\_Sachin\mPFC\CCG analysis\newCCG2\pv_pv2.mat')
ifinvert = [1 0 0 0 0 0 1 1 0];
for iC = 1:size(CCR,1),
    if ifinvert(iC) == 1,
        nCCR(iC,:) = (fliplr(CCR(iC,:)) - MeanH0(iC)) / SDH0(iC);%        (nanmean(CCR(iC,:))));
    else
        nCCR(iC,:) =  (CCR(iC,:) - MeanH0(iC)) / SDH0(iC);
    end
end

% inx = [];
% poc = {};
% for k=1:size(nCCR,1)
%     poc{k,1}=PairOfCells{k}{1};
%     poc{k,2}=PairOfCells{k}{2};
%     [r,s,t1,u] = cellid2tags(poc{k,1});
%     [r,s,t2,u] = cellid2tags(poc{k,2});
%     if isequal(t1,t2)
%         inx = [inx k];
%     end
% end
% nCCR(inx,:) = [];
% CCR(inx,:) = [];
% LCCR(inx,:) = [];
% UCCR(inx,:) = [];
% size(nCCR,1)

for k = 1:size(nCCR,1)
    if any(CCR(k,:)<LCCR(k,:))
        figure;plot(CCR(k,:),'k')
        hold on
        plot(LCCR(k,:),'r')
    end
end
% for k = 1:size(nCCR,1)
%     if any(CCR(k,:)>UCCR(k,:))
%         figure;plot(CCR(k,:),'k')
%         hold on
%         plot(UCCR(k,:),'r')
%     end
% end
% for k = 1:size(nCCR,1)
%     if any(CCR(k,:)>UCCR(k,:)) && any(CCR(k,:)<LCCR(k,:))
%         figure;plot(CCR(k,:),'k')
%         hold on
%         plot(UCCR(k,:),'r')
%         plot(LCCR(k,:),'r')
%     end
% end

subplot(221)
time_series = [-30:30];
SECCR=std(nCCR)/sqrt(size(nCCR,1));
errorbar(nanmean(nCCR),SECCR)
ylim([-7 9]);
clear CCR nCCR LCCR UCCR PairOfCells UniqueCellPairs MeanH0 SDH0

load('C:\Dropbox\KepecsLab\_Sachin\Analyses\newCCG2\pv_nontagged_H0.mat')
load('C:\Dropbox\KepecsLab\_Duda\Analyses\newCCG\pv_nontagged2.mat')


for iC = 1:size(CCR,1),
    nCCR(iC,:) = (CCR(iC,:) - MeanH0(iC)) / SDH0(iC); %        (nanmean(CCR(iC,:))));
end

nCCR(186:196,:) = [];

inx = [];
poc = {};
for k=1:size(nCCR,1)
    poc{k,1}=PairOfCells{k}{1};
    poc{k,2}=PairOfCells{k}{2};
    [r,s,t1,u] = cellid2tags(poc{k,1});
    [r,s,t2,u] = cellid2tags(poc{k,2});
    if isequal(t1,t2)
        inx = [inx k];
    end
end
nCCR(inx,:) = [];
CCR(inx,:) = [];
LCCR(inx,:) = [];
UCCR(inx,:) = [];
% nCCR = nCCR(inx,:);
% CCR = CCR(inx,:);
% LCCR = LCCR(inx,:);
% UCCR = UCCR(inx,:);
size(nCCR,1)

% for k = 1:size(nCCR,1)
%     if any(CCR(k,:)<LCCR(k,:))
%         figure;plot(CCR(k,:),'k')
%         hold on
%         plot(LCCR(k,:),'r')
%     end
% end
for k = 1:size(nCCR,1)
    if any(CCR(k,:)>UCCR(k,:))
        figure;plot(CCR(k,:),'k')
        hold on
        plot(UCCR(k,:),'r')
    end
end
% for k = 1:size(nCCR,1)
%     if any(CCR(k,:)>UCCR(k,:)) && any(CCR(k,:)<LCCR(k,:))
%         figure;plot(CCR(k,:),'k')
%         hold on
%         plot(UCCR(k,:),'r')
%         plot(LCCR(k,:),'r')
%     end
% end


subplot(222)
time_series = [-30:30];
SECCR=std(nCCR)/sqrt(size(nCCR,1));
errorbar(nanmean(nCCR),SECCR);
ylim([-2 2]);
%  stairs(time_series,nanmean(nCCR),'color','r');
%plot(time_series,nanmean(nCCR),'color','r');
% set(sub_h(2),'box','off','xlim',xlim,'xtick',xtick,'TickDir','out','xticklabel',[],'ylim',ylim,'ytick',ytick,'yticklabel',[],'fontsize',8)
clear nCCR CCR LCCR UCCR PairOfCells UniqueCellPairs MeanH0 SDH0


%som-som pairs
load('C:\Dropbox\KepecsLab\_Duda\Analyses\newCCG\som_som_H0.mat')
load('C:\Dropbox\KepecsLab\_Duda\Analyses\newCCG\som_som2.mat')
ifinvert = zeros(1,15);ifinvert(1,7)=1;
tinx = [3,4,6,7,8,10,11,12,15];
for iC =  1:size(CCR,1),
        nCCR(iC,:) = (CCR(iC,:) - MeanH0(iC)) / SDH0(iC);%        (nanmean(CCR(iC,:))));
%         nCCR(iC,:) = CCR(iC,:);
%         sCCR(iC,:)=filtfilt(b,a,nCCR(iC,:));
%         nCCR(iC,:) = sCCR(iC,:)./(nanmean(sCCR(iC,:)));
        %     end
end
% nCCR=nCCR(any(nCCR,2),:);

% inx = [];
% poc = {};
% for k=1:size(nCCR,1)
%     poc{k,1}=PairOfCells{k}{1};
%     poc{k,2}=PairOfCells{k}{2};
%     [r,s,t1,u] = cellid2tags(poc{k,1});
%     [r,s,t2,u] = cellid2tags(poc{k,2});
%     if isequal(t1,t2)
%         inx = [inx k];
%     end
% end
% nCCR(inx,:) = [];
% size(nCCR,1)

for k = 1:size(nCCR,1)
%     if any(CCR(k,:)>UCCR(k,:)) && any(CCR(k,:)<LCCR(k,:))
        figure;plot(CCR(k,:),'k')
        hold on
        plot(UCCR(k,:),'r')
        plot(LCCR(k,:),'r')
%     end
end

subplot(223)
time_series = [-30:30];
SECCR=std(nCCR)/sqrt(size(nCCR,1));
errorbar(nanmean(nCCR),SECCR);
ylim([-7 9]);
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
    
inx = [];
poc = {};
for k=1:size(nCCR,1)
    poc{k,1}=PairOfCells{k}{1};
    poc{k,2}=PairOfCells{k}{2};
    [r,s,t1,u] = cellid2tags(poc{k,1});
    [r,s,t2,u] = cellid2tags(poc{k,2});
    if isequal(t1,t2)
        inx = [inx k];
    end
end
nCCR(inx,:) = [];
CCR(inx,:) = [];
LCCR(inx,:) = [];
UCCR(inx,:) = [];
% nCCR = nCCR(inx,:);
% CCR = CCR(inx,:);
% LCCR = LCCR(inx,:);
% UCCR = UCCR(inx,:);
PairOfCells(inx) = [];
size(nCCR,1)

% for k = 1:size(nCCR,1)
%     if any(CCR(k,:)<LCCR(k,:))
%         figure;plot(CCR(k,:),'k')
%         hold on
%         plot(LCCR(k,:),'r')
%         k
%     end
% end
for k = 1:size(nCCR,1)
    if any(CCR(k,:)>UCCR(k,:))
        figure;plot(CCR(k,:),'k')
        hold on
        plot(UCCR(k,:),'r')
    end
end
% for k = 1:size(nCCR,1)
%     if any(CCR(k,:)>UCCR(k,:)) && any(CCR(k,:)<LCCR(k,:))
%         figure;plot(CCR(k,:),'k')
%         hold on
%         plot(UCCR(k,:),'r')
%         plot(LCCR(k,:),'r')
%     end
% end

subplot(224)
time_series = [-30:30];
SECCR=std(nCCR)/sqrt(size(nCCR,1));
errorbar(nanmean(nCCR),SECCR);
ylim([-2 2]);