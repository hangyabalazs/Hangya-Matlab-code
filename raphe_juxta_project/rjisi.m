function rjisi
%RJISI   Interspike interval analysis (raphe juxta project).
%   RJISI calculates and saves the following outputs for theta and 
%   non-theta segments:
%       minimal (rank = 1) interspike interval (ISI)
%       mean of ISI with ranks 1-10
%       mean of ISI with ranks 11-20
%       entropy of ISI distribution
%
%   ISI distribution plots are saved in fig format.
%
%   See also RJSTA.

% Directories
global DATADIR
global DATAPATH
inpdir = [DATADIR 'raphe_matfiles\raphe_juxta_files_discriminated\temp\'];
thetadir = [DATAPATH 'Raphe\raphe_juxta\Wavelet_temp\theta_segments\'];
nothdir = [DATAPATH 'Raphe\raphe_juxta\Wavelet_temp\nontheta_segments\'];
sharpdir = [DATAPATH 'Raphe\raphe_juxta\Wavelet_temp\sharpwave_segments\'];
resdir = [DATAPATH 'Raphe\raphe_juxta\ISI2c\'];
mm = pwd;
cd(resdir)

% Filelist
[files files_short] = b_filelist(inpdir);
sf = length(files_short);

% Progress indicator
wb = waitbar(0,'Running RJISI...','Position',[360 250 275 50]);
global WB
WB(end+1) = wb;

% Main
sr = 10000;     % sampling rate
nqf = sr / 2;   % Nyquist frequency
dsr = 1000;     % rate for downsampling
cnst = sr / dsr;
xlsout = cell(sf,7);    % initialize Excel output
for o = 1:sf
    fname = files(o).name;
    cmps = strread(fname,'%s','delimiter','_');
    titlestr = [cmps{1} ' ' cmps{2}];
    ff = [inpdir fname];
    load(ff)
    
    ff = [thetadir 'THETA_SEGMENTS_' cmps{1} '_' cmps{2}(1:2)];  % phase; load theta segments
    load(ff)
    if ~isempty(ThetaSegments)
        ThetaSegments = unitesegt(ThetaSegments,sr);     % drop gaps < 0.5 s
        ThetaSegments = short_killer(ThetaSegments);    % drop segments < 3 s
        lethe = ThetaSegments(2,:) - ThetaSegments(1,:);  % find the longest theta segment
        minx = find(lethe==max(lethe));
        minx1 = ThetaSegments(1,minx);
        minx2 = ThetaSegments(2,minx);
        minx2 = min(minx2,length(eeg));
        tnc = std(eeg(minx1:cnst:minx2));    % normalizing constant for STA
    end
    th_index = size(ThetaSegments,2);
    VdiscTheta = [];
    for t = 1:th_index      % theta segment cycle
        th1 = ThetaSegments(1,t);
        th2 = ThetaSegments(2,t);
        th2 = min(th2,length(eeg));
        VdiscTheta = [VdiscTheta vdisc(vdisc>th1&vdisc<th2)];
    end
    IsiTheta = diff(VdiscTheta);
    SVT = sort(IsiTheta,'ascend');    % ISI distribution
    edges = [0:500:10000];
    cnts = (edges(1:end-1) + edges(2:end)) / 2;
    HistT = histc(IsiTheta,edges);
    HistT = HistT(1:end-1);
    H = figure;
    if isempty(HistT)
        HistT = zeros(1,length(cnts));
        SVT = zeros(1,100);
    end
    bar(cnts,HistT)
    fnisi = [fname(1:end-4) '_theta_ISI'];
    saveas(H,fnisi);
    HistT = HistT / sum(HistT);
    ht = HistT(HistT~=0);
    EntT = - sum(ht.*log2(ht));
    Nt = length(IsiTheta);   % total number of trials
    R_bar = length(HistT);    % number of possible responses
    Bias_HRT = ((-1) / (2 * Nt * log(2))) * (R_bar - 1)    % bias correction
%     EntT = EntT - Bias_HRT;

    ff = [nothdir 'NONTHETA_SEGMENTS_' cmps{1} '_' cmps{2}(1:2)];  % load non-theta segments
    load(ff)
    if ~isempty(NonThetaSegments)
        NonThetaSegments = unitesegn(NonThetaSegments,sr);     % drop gaps < 0.5 s
        NonThetaSegments = short_killer(NonThetaSegments);    % drop segments < 3 s
        lenon = NonThetaSegments(2,:) - NonThetaSegments(1,:);  % find the longest non-theta segment
        minx = find(lenon==max(lenon));
        minx1 = NonThetaSegments(1,minx);
        minx2 = NonThetaSegments(2,minx);
        minx2 = min(minx2,length(eeg));
        nnc = std(eeg(minx1:cnst:minx2));    % normalizing constant for STA
    end
    no_index = size(NonThetaSegments,2);
    VdiscNonTheta = [];
    for t = 1:no_index      % non-theta segment cycle
        no1 = NonThetaSegments(1,t);
        no2 = NonThetaSegments(2,t);
        no2 = min(no2,length(eeg));
        vdisc_noth = vdisc(vdisc>no1&vdisc<no2) - no1;
        VdiscNonTheta = [VdiscNonTheta vdisc(vdisc>no1&vdisc<no2)];
    end
    IsiNonTheta = diff(VdiscNonTheta);
    SVNT = sort(IsiNonTheta,'ascend');    % ISI distribution
    edges = [0:500:10000];
    cnts = (edges(1:end-1) + edges(2:end)) / 2;
    HistN = histc(IsiNonTheta,edges);
    HistN = HistN(1:end-1);
    if isempty(HistN)
        HistN = zeros(1,length(cnts));
        SVNT = zeros(1,100);
    end
    H = figure;
    bar(cnts,HistN)
    fnisi = [fname(1:end-4) '_nontheta_ISI'];
    saveas(H,fnisi);
    HistN = HistN / sum(HistN);
    ht = HistN(HistN~=0);
    EntN = - sum(ht.*log2(ht));
    Nt = length(IsiNonTheta);   % total number of trials
    R_bar = length(HistN);    % number of possible responses
    Bias_HRN = ((-1) / (2 * Nt * log(2))) * (R_bar - 1)    % bias correction
%     EntN = EntN - Bias_HRN;
        
    xlsout{o,1} = fname(1:end-4);   % Excel output
    xlsout{o,2} = SVT(1);       % min ISI
    if length(SVT) >= 10
        xlsout{o,3} = mean(SVT(1:10));  % ISI 1-10
    else
        xlsout{o,3} = NaN;
    end
    if length(SVT) >= 20
        xlsout{o,4} = mean(SVT(11:20)); % ISI 11-20
    else
        xlsout{o,4} = NaN;
    end
    xlsout{o,5} = EntT;         % entropy
    xlsout{o,6} = Bias_HRT;     % entropy bias
    xlsout{o,7} = SVNT(1);      % min ISI
    if length(SVNT) >= 10
        xlsout{o,8} = mean(SVNT(1:10));  % ISI 1-10
    else
        xlsout{o,8} = NaN;
    end
    if length(SVNT) >= 20
        xlsout{o,9} = mean(SVNT(11:20)); % ISI 11-20
    else
        xlsout{o,9} = NaN;
    end
    xlsout{o,10} = EntN;         % entropy
    xlsout{o,11} = Bias_HRN;     % entropy bias
    
    close all
    waitbar(o/sf)
end
xlswrite('isisummary',xlsout)
close(wb)
cd(mm)



% -------------------------------------------------------------------------
function segments2 = unitesegt(segments,sr)

len = size(segments,2);
segments2 = segments;
for k = 1:len-1
    la = segments(1,k+1);
    fi = segments(2,k);
    if (la-fi)/sr < 0.5
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
function segments = short_killer(segments)

% Skip short segments
int = segments;
int1 = int(1,:);
int2 = int(2,:);
difint = int2 - int1;
fd = find(difint<30000);         % leaving segments shorter than 3 sec.
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
function H = stafig(sta,sta_amp,sta_index2,nn,wn,sr,titlestr)

time = linspace(-wn/sr/2,wn/sr/2,length(sta));
H = figure;
plot(time,sta,'LineWidth',1.5)
ach = allchild(H);     % figure title
ax = findobj(ach,'type','axes');
title(ax(end),titlestr)
x_lim = xlim;
y_lim = ylim;
str = ['\it{Amp: }' '\bf ' num2str(sta_amp)];
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