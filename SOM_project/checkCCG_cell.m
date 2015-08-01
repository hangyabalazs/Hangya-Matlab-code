%% load

load('C:\Balazs\_analysis\SOM\nt_pairs\CCG_matrices1.mat')

%% all 0: CCG not calculated (not enough spikes or no segments)

length(find(all(CCR'==0)))

%% too low count

length(find(mean(CCR,2)<5))
length(find(mean(CCR,2)>5))

%% inhibition

% inh = find(any(CCR(:,27:35)<LCCR(:,27:35),2)&mean(LCCR(:,27:35),2)>1);
inh = find(any(CCR(:,27:35)<LCCR(:,27:35),2)&mean(CCR,2)>5);
inh = find(any(CCR(:,27:35)<LCCR(:,27:35),2)&mean(CCR,2)>5);
inh = find(any(CCR(:,27:35)<LCCR(:,27:35),2));

%%

for k = 1:11
    figure(inh(k))
    plot(CCR(inh(k),:))
    hold on
    plot(LCCR(inh(k),:),'Color',[0.7 0.7 0.7])
    plot(UCCR(inh(k),:),'Color',[0.7 0.7 0.7])
    title(num2str(mean(CCR(inh(k),:))))
end

%%

may_include = [240 311 365 673 926 944 945 1068 1078 1202 1545 2096];
looks_real = [103 107 541 585 642 643 671 672 677 679 682 685 838 846 849 889 890 925 933 942 946 948 961 975 986 987 991 ...
    1045 1069 1072 1073 1161 1199 1212 1216 1219 1220 1222 1263 1270 1287 1289 ...
    1343 1365 1371 1381 1587 1693 1699 1728 1738 2922 2933 2944 2948 2951 3035 3179 3201 3243 3334];
definitely_real = [103 541 643 672 677 682 685 1371 3334];

%% excitation

exc = find(any(CCR(:,27:35)>UCCR(:,27:35),2));

for k = 21:115
    figure(exc(k))
    plot(CCR(exc(k),:))
    hold on
    plot(LCCR(exc(k),:),'Color',[0.7 0.7 0.7])
    plot(UCCR(exc(k),:),'Color',[0.7 0.7 0.7])
    title(num2str(mean(CCR(exc(k),:))))
end

%% common contamination?

PairOfCells(679,:)
PairOfCells(671,:)
PairOfCells(642,:)
PairOfCells(585,:)
PairOfCells(107,:)
PairOfCells(3242,:)
PairOfCells(3243,:)
PairOfCells(3201,:)
PairOfCells(3179,:)
PairOfCells(3035,:)
PairOfCells(3007,:)
PairOfCells(1738,:)
PairOfCells(1728,:)
PairOfCells(1669,:)
PairOfCells(1693,:)
PairOfCells(1587,:)
PairOfCells(1381,:)
PairOfCells(1365,:)
PairOfCells(1343,:)
PairOfCells(1287,:)
PairOfCells(1220,:)
PairOfCells(1216,:)
PairOfCells(1161,:)
PairOfCells(1069,:)
PairOfCells(991,:)
PairOfCells(986,:)
PairOfCells(975,:)
PairOfCells(946,:)
PairOfCells(925,:)
PairOfCells(889,:)
PairOfCells(849,:)
PairOfCells(846,:)

%%

real_inh = [5 103 107 139 240 311 362 365 439 505];
symmetric_inh = [103 107];

%% inhibition

inh = find(any(CCR(:,20:41)<LCCR(:,20:41),2));

for k = 1:20
    figure(inh(k))
    plot(CCR(inh(k),:))
    hold on
    plot(LCCR(inh(k),:),'Color',[0.7 0.7 0.7])
    plot(UCCR(inh(k),:),'Color',[0.7 0.7 0.7])
end


%%

real_inh = [5 103 107 139 240 311 362 365 439 505];
symmetric_inh = [103 107];