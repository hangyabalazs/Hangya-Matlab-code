%% load

load('C:\Balazs\_analysis\SOM\nt_pairs\CCG_matrices1.mat')

%% all 0: CCG not calculated (not enough spikes or no segments)

length(find(all(CCR'==0))) % 1468 good pairs if excluding duplicates

%% inhibition

inh = find(any(CCR(:,27:35)<LCCR(:,27:35),2));

%%

for k = 461:480
    figure(inh(k))
    plot(CCR(inh(k),:))
    hold on
    plot(LCCR(inh(k),:),'Color',[0.7 0.7 0.7])
    plot(UCCR(inh(k),:),'Color',[0.7 0.7 0.7])
    title(num2str(mean(CCR(inh(k),:))))
end

%%

inhibited = [103 311 365 541 643 672 673 677 682 685 838 890 ...
    926 942 944 945 948 987 1045 1068 1072 1078 1199 1222 ...
    1260 1263 1270 1289 1332 1371 1545 ...
    2096 2105 2135 2750 2922 2933 2944 2948 2951 2963 ...
    3334];

%% excitation

exc = find(any(CCR(:,27:35)>UCCR(:,27:35),2));

%%

for k = 1:115
    figure(exc(k))
    plot(CCR(exc(k),:))
    hold on
    plot(LCCR(exc(k),:),'Color',[0.7 0.7 0.7])
    plot(UCCR(exc(k),:),'Color',[0.7 0.7 0.7])
    zoom on
    title(num2str(mean(CCR(exc(k),:))))
end

%%

excited = [139 140 397 535 540 545 553 565 674 677 682 685 688 ...
    696 705 785 800 807 876 893 906 937 948 ...
    1065 1117 1162 1260 ...
    1379 1382 1752 1769 1785 1830 1875 1939 1958 2080 2114 2132 ...
    2148 2532 2565 2571 2807 2831 2875 2944 3049 3111 3120 3163 ...
    3223 3368 3371 3407 3682 3757 3817];
sym_excited = [139 545 553 565 674 ...
    696 705 785 800 807 893 937 ...
    1065 1260 ...
    1752 ...
    2532 2807 2875 ...
    3371];


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