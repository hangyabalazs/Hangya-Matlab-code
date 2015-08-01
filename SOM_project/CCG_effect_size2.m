function CCG_effect_size2
%CCG_EFFECT_SIZE   Relative effect size for cross-correlograms.
%   CCG_EFFECT_SIZE quantifies negative/positive peak value (min/max) in
%   the -4 to +4 ms window with respect to the average count outside this
%   window (difference).
%
%   See also SOM_CCG_CONF_FILTER and CCG_EFFECT_SIZE2.

% Effect size for PV-PV
load ('C:\Dropbox\KepecsLab\_Sachin\mPFC\CCG analysis\CCG3\pv_pv_H02');   % load CCGs
excludeinx = [];
[pv_inh pv_exc inx] = efsize(CCR,LCCR,UCCR,excludeinx);
pv_inh = pv_inh(~isnan(pv_inh));   % remove NaNs
pv_exc = pv_exc(~isnan(pv_exc));   % remove NaNs

% Effect size for PV-NT
load ('C:\Dropbox\KepecsLab\_Sachin\mPFC\CCG analysis\CCG3\pv_nontagged_H02');   % load CCGs
excludeinx = [61 123];
[pvnt_inh pvnt_exc inx] = efsize(CCR,LCCR,UCCR,excludeinx);
pvnt_inh = pvnt_inh(~isnan(pvnt_inh));   % remove NaNs
pvnt_exc = pvnt_exc(~isnan(pvnt_exc));   % remove NaNs

% Effect size for SOM-NT
load('C:\Dropbox\KepecsLab\_Sachin\mPFC\CCG analysis\CCG3\som_nontagged_H02.mat')   % load CCGs
excludeinx = [1 2 12 39 61 73 89 92 129 142 143];
[somnt_inh somnt_exc inx] = efsize(CCR,LCCR,UCCR,excludeinx);
somnt_inh = somnt_inh(~isnan(somnt_inh));   % remove NaNs
somnt_exc = somnt_exc(~isnan(somnt_exc));   % remove NaNs

% Effect size for NT-NT
load('C:\Balazs\_analysis\SOM\nt_pairs\CCG_matrices3.mat')    % load CCGs
excited = [30 105 114 139 178 246 321 349 379 472 524 527 528 531 554 565 586 587 674 ...
    677 682 685 688 705 777 782 784 785 787 788 792 795 798 800 807 814 857 893 912 943 951 952 962 969 976 997 1043 1066 1975 1123 1125 1156 ...
    1202 1208 1212 1214 1217 1218 1222 1228 1256 1258 1260 1270 1271 1272 1273 1274 1282 1285 1286 1289 1297 1308 1335 1340 1342 1350 1352 1363 1366 1370 1372 1379 1380 ...
    1382 1383 1425 1476 1491 1577 1642 1648 1658 1787 1812 1830 1833 1869 1879 1892 ...
    1902 1939 1995 1996 1998 2026 2028 2029 2033 2042 2069 2083 2085 2086 2094 2095 2096 2102 2103 2107 2128 2132 2134 2135 2138 2140 2148 2203 2303 2306 2308 2326 2329 ...
    2444 2543 2545 2553 2719 2720 2737 2739 2793 2797 2801 2802 2807 2811 2814 2815 ...
    2821 2825 2830 2831 2875 2891 2914 2922 2924 2929 2931 2932 2933 2934 2935 2936 2942 2943 2944 2947 2948 2949 2951 2965 2973 2985 3022 3070 3121 3127 3205 ...
    3313 3318 3319 3323 3331 3357 3385 3395 3396 3409 3484 3500];
inhibited = [103 217 570 586 587 595 643 672 673 677 ...
     682 685 779 781 784 796 799 822 825 837 850 890 891 894 911 919 923 924 926 928 929 933 934 945 947 960 962 976 987 988 989 990 1004 1060 1067 ...
     1068 1073 1075 1080 1113 1116 1156 1160 1168 1188 1191 1192 1202 1205 1206 1214 1215 1222 1257 1259 1270 1274 1286 1289 1291 1297 1340 1348 1350 1371 1375 ...
     1586 1599 1604 1660 1760 1771 1789 1792 1815 1832 ...
     1879 1892 1939 1980 1992 2029 2033 2034 2064 2070 2131 2140 2240 2287 ...
     2596 2739 2744 2749 2750 2831 2888 2890 2914 2920 2922 2931 2934 2942 2944 2948 2951 2980 2998 3075 ...
     3083 3127 3145 3170 3175 3194 3195 3200 3206 3207 3276 3313 3319];
[ntnt_inh ntnt_exc] = efsize2(CCR,inhibited,excited);
ntnt_inh = ntnt_inh(~isnan(ntnt_inh));   % remove NaNs
ntnt_exc = ntnt_exc(~isnan(ntnt_exc));   % remove NaNs

% Plot bar graph with mean and SE
figure
hold on
bar(1,mean(pv_exc),'FaceColor','w','EdgeColor',[0 0 0],'LineWidth',2)  % PV-PV, excitation
bar(1,mean(pv_inh),'FaceColor','w','EdgeColor',[0 0 0],'LineWidth',2)  % PV-PV, inhibition
E = errorbar([1 1],[mean(pv_exc) mean(pv_inh)],...
    [std(pv_exc)/sqrt(length(pv_exc)) std(pv_inh)/sqrt(length(pv_inh))],...
    '+','Color',[0 0 0],'LineWidth',2);   % SE
errorbar_tick(E,0)   % eliminate horizontal line from errorbar

bar(2,mean(pvnt_exc),'FaceColor','w','EdgeColor',[0 0 0],'LineWidth',2)  % PV-NT, excitation
bar(2,mean(pvnt_inh),'FaceColor','w','EdgeColor',[0 0 0],'LineWidth',2)  % PV-NT, inhibition
E = errorbar([2 2],[mean(pvnt_exc) mean(pvnt_inh)],...
    [std(pvnt_exc)/sqrt(length(pvnt_exc)) std(pvnt_inh)/sqrt(length(pvnt_inh))],...
    '+','Color',[0 0 0],'LineWidth',2);   % SE
errorbar_tick(E,0)   % eliminate horizontal line from errorbar

bar(3,mean(somnt_exc),'FaceColor','w','EdgeColor',[0 0 0],'LineWidth',2)  % SOM-NT, excitation
bar(3,mean(somnt_inh),'FaceColor','w','EdgeColor',[0 0 0],'LineWidth',2)  % SOM-NT, inhibition
E = errorbar([3 3],[mean(somnt_exc) mean(somnt_inh)],...
    [std(somnt_exc)/sqrt(length(somnt_exc)) std(somnt_inh)/sqrt(length(somnt_inh))],...
    '+','Color',[0 0 0],'LineWidth',2);   % SE
errorbar_tick(E,0)   % eliminate horizontal line from errorbar

bar(4,mean(ntnt_exc),'FaceColor','w','EdgeColor',[0 0 0],'LineWidth',2)  % NT-NT, excitation
bar(4,mean(ntnt_inh),'FaceColor','w','EdgeColor',[0 0 0],'LineWidth',2)  % NT-NT, inhibition
E = errorbar([4 4],[mean(ntnt_exc) mean(ntnt_inh)],...
    [std(ntnt_exc)/sqrt(length(ntnt_exc)) std(ntnt_inh)/sqrt(length(ntnt_inh))],...
    '+','Color',[0 0 0],'LineWidth',2);   % SE
errorbar_tick(E,0)   % eliminate horizontal line from errorbar

keyboard

% Plot bar graph with median and IQR
figure
hold on
bar(1,median(pv_exc),'FaceColor','w','EdgeColor',[0 0 0],'LineWidth',2)  % PV-PV, excitation
bar(1,median(pv_inh),'FaceColor','w','EdgeColor',[0 0 0],'LineWidth',2)  % PV-PV, inhibition
E = errorbar([1 1],[median(pv_exc) median(pv_inh)],...
    [prctile(pv_exc,25) prctile(pv_inh,25)],...
    [prctile(pv_exc,75) prctile(pv_inh,75)],...
    '+','Color',[0 0 0],'LineWidth',2);   % IQR
errorbar_tick(E,0)   % eliminate horizontal line from errorbar

bar(2,median(pvnt_exc),'FaceColor','w','EdgeColor',[0 0 0],'LineWidth',2)  % PV-NT, excitation
bar(2,median(pvnt_inh),'FaceColor','w','EdgeColor',[0 0 0],'LineWidth',2)  % PV-NT, inhibition
E = errorbar([2 2],[median(pvnt_exc) median(pvnt_inh)],...
    [prctile(pvnt_exc,25) prctile(pvnt_inh,25)],...
    [prctile(pvnt_exc,75) prctile(pvnt_inh,75)],...
    '+','Color',[0 0 0],'LineWidth',2);   % IQR
errorbar_tick(E,0)   % eliminate horizontal line from errorbar

bar(3,median(somnt_exc),'FaceColor','w','EdgeColor',[0 0 0],'LineWidth',2)  % SOM-NT, excitation
bar(3,median(somnt_inh),'FaceColor','w','EdgeColor',[0 0 0],'LineWidth',2)  % SOM-NT, inhibition
E = errorbar([3 3],[median(somnt_exc) median(somnt_inh)],...
    [prctile(somnt_exc,25) prctile(somnt_inh,25)],...
    [prctile(somnt_exc,75) prctile(somnt_inh,75)],...
    '+','Color',[0 0 0],'LineWidth',2);   % IQR
errorbar_tick(E,0)   % eliminate horizontal line from errorbar

keyboard

% -------------------------------------------------------------------------
function [mins maxs inx] = efsize(CCR,LCCR,UCCR,excludeinx)

mins = [];  % effect size for inhibition
maxs = [];  % effect size for excitation
inx = [];
for k = 1:size(CCR,1)   % loop through all CCGs
    if ~ismember(k,excludeinx) && ...   % not excluded
            (any(CCR(k,27:35)<LCCR(k,27:35)) || any(CCR(k,27:35)>UCCR(k,27:35)))   % if there was an effect
        minpc = - mean(CCR(k,[1:26 36:61])) + min(CCR(k,27:35));   % diff. between mean and trough
%         minpc = min(CCR(k,27:35)) / mean(CCR(k,[1:26 36:61]));
        maxpc = max(CCR(k,27:35)) - mean(CCR(k,[1:26 36:61]));   % diff. between peak and mean
%         maxpc = max(CCR(k,27:35)) / mean(CCR(k,[1:26 36:61]));
        if any(CCR(k,27:35)<LCCR(k,27:35))   % inhibition
            mins(end+1) = minpc;  % effect size for inhibition
        else
            mins(end+1) = NaN;
        end
        if any(CCR(k,27:35)>UCCR(k,27:35))   % excitaion
            maxs(end+1) = maxpc;  % effect size for excitation
        else
            maxs(end+1) = NaN;
        end
        inx(end+1) = k;
    end
end

% -------------------------------------------------------------------------
function [mins maxs] = efsize2(CCR,inhinx,excinx)

mins = [];  % effect size for inhibition
maxs = [];  % effect size for excitation
for k = 1:size(CCR,1)   % loop through all CCGs
    if ismember(k,inhinx)   % inhibited
        minpc = - mean(CCR(k,[1:26 36:61])) + min(CCR(k,27:35));   % diff. between mean and trough
        mins(end+1) = minpc;  % effect size for inhibition
    end
    if ismember(k,excinx)   % excited
        maxpc = max(CCR(k,27:35)) - mean(CCR(k,[1:26 36:61]));   % diff. between peak and mean
        maxs(end+1) = maxpc;  % effect size for excitation
    end
end