function rjtp
%RJTP   Analysis of non-theta to theta transitions (raphe juxta project).
%   RJTP plots firing rate changes around transitions from non-theta to
%   theta segments evoked by sensory stimulation (tail-pinch, TP) or
%   occurring spontaneously. TPs that do not induce theta (within 3s) are
%   treated separately. Last 5 seconds of pre-transition non-theta segments
%   or pre-TP segments for TPs that do not evoke theta oscillation are used
%   as controls. First 3 seconds of theta segments or post-TP segments are
%   analyzed. Statistical testing is performed by using Wilcoxon's ranksum
%   test.
%
%   See also RJFRATE.

% Directories
global DATADIR
global DATAPATH
inpdir = [DATADIR 'raphe_matfiles\raphe_juxta_files_discriminated\temp\'];
inpdir_unit = [DATADIR 'raphe_matfiles\raphe_juxta_files\temp\'];
thetadir = [DATAPATH 'Raphe\raphe_juxta\Wavelet_temp\theta_segments\'];
nothdir = [DATAPATH 'Raphe\raphe_juxta\Wavelet_temp\nontheta_segments\'];
resdir = [DATAPATH 'Raphe\raphe_juxta\Frate_new2c\'];
tpdir = [DATAPATH 'Raphe\raphe_juxta\'];
mm = pwd;
cd(resdir)

% Filelist
[files files_short] = b_filelist(inpdir);
sf = length(files_short);

% Progress indicator
wb = waitbar(0,'Running RJTP...','Position',[360 250 275 50]);
global WB
WB(end+1) = wb;

% Tail pinch intervals
fn = [tpdir 'TPs.xls'];
headerrows = 1;
[mtx ntx atx] = xlsread(fn,'Sheet1');
ntx(1:headerrows,:) = [];
atx(1:headerrows,:) = [];

% Main
sr = 10000;     % sampling rate
for o = 1:sf
    fname = files(o).name;
    cmps = strread(fname,'%s','delimiter','_');
    titlestr = [cmps{1} ' ' cmps{2}];
    filenam = [cmps{1} '_' cmps{2}];
    ff = [inpdir fname];    % load EEG, unit, disc. unit
    load(ff)
    ff = [inpdir_unit fname(1:end-6) '.mat'];
    load(ff)
    unit = Unit.values;
    
    inx = find(strcmp(filenam,atx(:,1)));     % find TPs
    TPstart = mtx(inx,1) * sr;
    TPend = mtx(inx,2) * sr;
    
    ff = [thetadir 'THETA_SEGMENTS_' cmps{1} '_' cmps{2}(1:2)];  % phase; load theta segments
    load(ff)
    if ~isempty(ThetaSegments)
        ThetaSegments = unitesegt(ThetaSegments,sr,0.5);     % drop gaps < 0.5 s
        ThetaSegments = short_killer(ThetaSegments,sr,3);    % drop segments < 3 s
    end
    th_index = size(ThetaSegments,2);
    ff = [nothdir 'NONTHETA_SEGMENTS_' cmps{1} '_' cmps{2}(1:2)];  % load non-theta segments
    load(ff)
    if ~isempty(NonThetaSegments)
        NonThetaSegments = unitesegn(NonThetaSegments,sr);     % drop gaps < 0.5 s
        NonThetaSegments = short_killer(NonThetaSegments,sr,5);    % drop segments < 5 s
    end
    for t = 1:th_index      % theta segment cycle
        th1 = ThetaSegments(1,t);
        th2 = ThetaSegments(2,t);
        th2 = min(th2,length(eeg));
        
        if any(ThetaSegments(2,:)<th1&ThetaSegments(2,:)>th1-10*sr)  % check if theta is preceded by non-theta
            continue
        end
        tpinx = find(TPstart<th1+3*sr);     % TP before or within 3 s from theta start
        if ~isempty(tpinx) && TPend(tpinx(end)) > th1
            thtype = 'evoked';
        else
            thtype = 'spont';
        end
        nthinx = find(NonThetaSegments(1,:)<th1,1,'last');
        if isempty(nthinx)
            continue
        end
        nth1 = NonThetaSegments(1,nthinx);
        nth2 = NonThetaSegments(2,nthinx);
        ntfr = zeros(1,5);
        for k = 5:-1:1   % control segment: non-theta > 5s preceding transition
            lvd = vdisc(vdisc>nth2-k*sr&vdisc<nth2-(k-1)*sr);
            ntfr(5-k+1) = length(lvd);
        end
        tfr = zeros(1,3);
        for k = 1:3     % theta segment > 3s
            lvd = vdisc(vdisc>th1+(k-1)*sr&vdisc<th1+k*sr);
            tfr(k) = length(lvd);
        end
        [Wp Wh] = b_ranksum2(tfr,ntfr,'alpha',0.05); % testing for sign. activation
        H = figure;     % plot
        subplot(3,1,1)  % plot EEG
        plinx = nth1:th2;
        plntinx = nth2-5*sr:nth2;
        plthinx = th1:th1+3*sr;
        leeg = eeg(plinx);
        plot(plinx,leeg)
        hold on
        plot(plntinx,eeg(plntinx),'r')  % control: red
        plot(plthinx,eeg(plthinx),'g')  % theta: green
        line([th1 th1],[min(leeg) max(leeg)],'Color','k')   % line: transition
        subplot(3,1,2)  % plot unit
        lunit = unit(plinx);
        plot(plinx,lunit)
        hold on
        plot(plntinx,unit(plntinx),'r')  % control: red
        plot(plthinx,unit(plthinx),'g')  % theta: green
        line([th1 th1],[min(lunit) max(lunit)],'Color','k')   % line: transition
        subplot(3,1,3)  % plot firing rate
        plot([ntfr tfr])
        hold on
        line([5 5],[min([ntfr tfr]) max([ntfr tfr])],'Color','k')   % line: transition
        if Wh   % text for stat.
            clr = 'red';
        else
            clr = 'black';
        end
        y_lim = ylim;
        x_lim = xlim;
        tpos1 = x_lim(1) + (x_lim(2) - x_lim(1)) / 2;
        tpos2 = y_lim(1) + (y_lim(2)-y_lim(1)) * 3 / 5;
        text(tpos1,tpos2,num2str(Wp),'Color',clr,'Horizontalalignment','center')
        tpos1 = x_lim(1) + (x_lim(2) - x_lim(1)) / 2;
        tpos2 = y_lim(1) + (y_lim(2)-y_lim(1)) * 4 / 5;
        text(tpos1,tpos2,thtype)
        fn = [fname(1:end-4) '_NT2TH' num2str(t)];  % save
        saveas(H,fn);
    end
    for t = 1:length(TPstart)
        if isnan(TPstart)
            continue
        end
        tp1 = TPstart(t);
        tp2 = TPend(t);
        if ~isempty(ThetaSegments)
            if any(ThetaSegments(1,:)>tp1&ThetaSegments(1,:)<tp1+5*sr) || ...
                    any(ThetaSegments(2,:)>tp1&ThetaSegments(2,:)<tp1+5*sr) || ...
                    any(ThetaSegments(1,:)<tp1&ThetaSegments(2,:)>tp1+5*sr) % check if theta is present after TP
                continue
            end
        end
        thtype = 'TP';   % TPs with no evoked theta
        for k = 5:-1:1   % control: preTP 5 s
            lvd = vdisc(vdisc>tp1-k*sr&vdisc<tp1-(k-1)*sr);
            ntfr(5-k+1) = length(lvd);
        end
        tfr = zeros(1,3);
        for k = 1:3   % 3s after TP
            lvd = vdisc(vdisc>tp1+(k-1)*sr&vdisc<tp1+k*sr);
            tfr(k) = length(lvd);
        end
        [Wp Wh] = b_ranksum2(tfr,ntfr,'alpha',0.05);    % test for sign. activation
        H = figure;     % plot
        subplot(3,1,1)  % plot EEG
        plinx = tp1-10*sr:tp1+10*sr;
        plntinx = tp1-5*sr:tp1;
        plthinx = tp1:tp1+3*sr;
        leeg = eeg(plinx);
        plot(plinx,leeg)
        hold on
        plot(plntinx,eeg(plntinx),'r')  % control: red
        plot(plthinx,eeg(plthinx),'g')  % after TP: green
        line([tp1 tp1],[min(leeg) max(leeg)],'Color','k')   % line: TP
        subplot(3,1,2)  % plot unit
        lunit = unit(plinx);
        plot(plinx,lunit)
        hold on
        plot(plntinx,unit(plntinx),'r')  % control: red
        plot(plthinx,unit(plthinx),'g')  % after TP: green
        line([tp1 tp1],[min(lunit) max(lunit)],'Color','k')   % line: TP
        subplot(3,1,3)  % plot firing rate
        plot([ntfr tfr])
        hold on
        line([5 5],[min([ntfr tfr]) max([ntfr tfr])],'Color','k')   % line: TP
        if Wh   % text for stat
            clr = 'red';
        else
            clr = 'black';
        end
        y_lim = ylim;
        x_lim = xlim;
        tpos1 = x_lim(1) + (x_lim(2) - x_lim(1)) / 2;
        tpos2 = y_lim(1) + (y_lim(2)-y_lim(1)) * 3 / 5;
        text(tpos1,tpos2,num2str(Wp),'Color',clr,'Horizontalalignment','center')
        tpos1 = x_lim(1) + (x_lim(2) - x_lim(1)) / 2;
        tpos2 = y_lim(1) + (y_lim(2)-y_lim(1)) * 4 / 5;
        text(tpos1,tpos2,thtype)
        fn = [fname(1:end-4) '_TP' num2str(t)];  % save
        saveas(H,fn);
    end
    close all
    waitbar(o/sf)
end
close(wb)
cd(mm)



% -------------------------------------------------------------------------
function segments2 = unitesegt(segments,sr,gp)

len = size(segments,2);
segments2 = segments;
for k = 1:len-1
    la = segments(1,k+1);
    fi = segments(2,k);
    if (la-fi)/sr < gp
        [fnx fny] = find(segments2==fi);
        segments2(fnx,fny) = segments(2,k+1);
        segments2 = [segments2(1:2,1:fny) segments2(1:2,fny+2:end)];
    end
end



% -------------------------------------------------------------------------
function segments2 = unitesegn(segments,sr)

len = size(segments,2);
segments2 = segments;
for k = 1:len-1
    la = segments(1,k+1);
    fi = segments(2,k);
    pre = segments(1,k);
    post = segments(2,k+1);
    if (la-fi)/sr < 0.5 && (fi-pre)/sr > 0.5 && (post-la)/sr > 0.5
        [fnx fny] = find(segments2==fi);
        segments2(fnx,fny) = segments(2,k+1);
        segments2 = [segments2(1:2,1:fny) segments2(1:2,fny+2:end)];
    end
end



% ----------------------------------------------------------------------------------
function segments = short_killer(segments,sr,lm)

% Skip short segments
int = segments;
int1 = int(1,:);
int2 = int(2,:);
difint = int2 - int1;
fd = find(difint<lm*sr);         % leaving segments shorter than 3 sec.
int1(fd) = [];
int2(fd) = [];
segments = [int1; int2];



% -------------------------------------------------------------------------
function sacr = racorr(ncc,bno)
%RACORR   Autocorrelation.
%   RACORR2(VD,BNO) calculates autocorrelogram for discriminated unit VD, 
%   using a +-1000 ms time window and BNO bins.

% Calculate spike times in milliseconds
sr = 1000;
nc = ncc * sr;

% Autocorrelogram
zunit1 = zeros(1,length(round(nc))+5);
zunit1(round(nc)) = 1;
acr = xcorr(zunit1,1*sr);
acr(length(acr)/2+0.5) = [];
acr = reshape(acr,length(acr)/bno,bno);     % window: -200 ms - 200 ms
sacr = sum(acr);



% -------------------------------------------------------------------------
function H = stafig(sta,sta_index1,sta_index2,nn,wn,sr,titlestr)

time = linspace(-wn/sr/2,wn/sr/2,length(sta));
H = figure;
plot(time,sta,'LineWidth',1.5)
ach = allchild(H);     % figure title
ax = findobj(ach,'type','axes');
title(ax(end),titlestr)
x_lim = xlim;
y_lim = ylim;
str = ['\it{Max-mean: }' '\bf ' num2str(sta_index1)];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
str = ['\it{Max: }' '\bf ' num2str(sta_index2)];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
str = ['\it{n: }' '\bf ' num2str(nn)];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')

% -------------------------------------------------------------------------
function [sta sta_index1 sta_index2 sta_amp lenv] = astanormmod(vdisc,eeg,wn,nc)
%ASTANORMMOD    Normalized Spike Triggered Average.
%   [S O1 O2 A L] = ASTANORM(VD,EEG,WN,NC) calculates STA of EEG and
%   discriminated unit (VD) using WN windowsize. Each EEG window is
%   normalized by NC before average calculation. STA, maximum STA, maximum
%   STA minus mean STA and STA amplitude are returned in S, O1, O2 and A.
%   Number of spikes is returned in L.

% Standardize EEG
eeg = (eeg - mean(eeg)) / nc;

% Calculate STA
wn2 = round(wn/2);
vdisc = vdisc(find(vdisc-wn2>0&vdisc+wn2<=length(eeg)));
lenv = length(vdisc);
st = zeros(lenv,2*wn2+1);
for t = 1:lenv
    eeg2 = eeg(vdisc(t)-wn2:vdisc(t)+wn2);
    st(t,1:2*wn2+1) = eeg2;
end
sta = mean(st,1);

% Output
sta_index1 = max(sta) - mean(sta);
sta_index2 = max(sta);
sta_amp = max(sta) - min(sta);