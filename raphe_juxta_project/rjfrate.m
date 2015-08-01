function rjfrate
%RJFRATE   Analysis of non-theta to theta transitions (raphe juxta project).
%   RJFRATE plots firing rate changes around transitions from non-theta to
%   theta segments.
%
%   It should be noted that a more precise determination of the transition
%   points would be necessary.
%
%   See also RJSTA and RJSTA2.

% Directories
global DATADIR
global DATAPATH
inpdir = [DATADIR 'raphe_matfiles\raphe_juxta_files_discriminated\'];
inpdir_unit = [DATADIR 'raphe_matfiles\raphe_juxta_files\'];
thetadir = [DATAPATH 'Raphe\raphe_juxta\Wavelet\theta_segments\'];
% nothdir = [DATAPATH 'Raphe\raphe_juxta\Wavelet\nontheta_segments\'];
% sharpdir = [DATAPATH 'Raphe\raphe_juxta\Wavelet\sharpwave_segments\'];
resdir = [DATAPATH 'Raphe\raphe_juxta\Frate\'];
mm = pwd;
cd(resdir)

% Filelist
[files files_short] = b_filelist(inpdir);
sf = length(files_short);

% Progress indicator
wb = waitbar(0,'Running RJFRATE...','Position',[360 250 275 50]);
global WB
WB(end+1) = wb;

% Main
sr = 10000;     % sampling rate
nqf = sr / 2;   % Nyquist frequency
% dsr = 1000;     % rate for downsampling
% cnst = sr / dsr;
% xlsout = cell(sf,7);    % initialize Excel output
for o = 1:sf
    fname = files(o).name;
    cmps = strread(fname,'%s','delimiter','_');
    titlestr = [cmps{1} ' ' cmps{2}];
    ff = [inpdir fname];
    load(ff)
    ff = [inpdir_unit fname(1:end-6) '.mat'];
    load(ff)
    unit = Unit.values;
    
    ff = [thetadir 'THETA_SEGMENTS_' cmps{1} '_' cmps{2}(1:2)];  % phase; load theta segments
    load(ff)
    if ~isempty(ThetaSegments)
        ThetaSegments = unitesegt(ThetaSegments,sr);     % drop gaps < 0.5 s
        ThetaSegments = short_killer(ThetaSegments);    % drop segments < 3 s
%         lethe = ThetaSegments(2,:) - ThetaSegments(1,:);  % find the longest theta segment
%         minx = find(lethe==max(lethe));
%         minx1 = ThetaSegments(1,minx);
%         minx2 = ThetaSegments(2,minx);
%         minx2 = min(minx2,length(eeg));
%         tnc = std(eeg(minx1:cnst:minx2));    % normalizing constant for STA
    end
%     flt = fir1(4096,[2/nqf 8/nqf]);  % filtering EEG
%     seeg = (eeg - mean(eeg)) / std(eeg);
%     feeg = filtfilt(flt,1,seeg);
    th_index = size(ThetaSegments,2);
%     VdiscTheta = [];
%     ThetaFiringRateC1 = zeros(1,th_index);
%     ThetaFiringRateC2 = zeros(1,th_index);
%     ThetaThetamod = zeros(1,th_index);
    figure
    wn = 10000;
    z = zeros(1,2*wn);
    for t = 1:th_index      % theta segment cycle
        th1 = ThetaSegments(1,t);
%         th2 = ThetaSegments(2,t);
%         th2 = min(th2,length(eeg));
%         eeg_theta = eeg(th1:th2);
%         feeg_theta = feeg(th1:th2);
%         vdisc_theta = vdisc(vdisc>th1&vdisc<th2) - th1;
%         VdiscTheta = [VdiscTheta vdisc(vdisc>th1&vdisc<th2)];
                
%         if ~isempty(vdisc_theta)        % firing rate and theta-modulation
%             efflen = (vdisc_theta(end) - vdisc_theta(1)) / sr;
%             ThetaFiringRateC1(t) = (length(vdisc_theta) - 1);
%             ThetaFiringRateC2(t) = efflen;
%             ac = racorr(vdisc_theta/10000,500);    % autocorrelation
%             [yac wac] = b_fft2(ac,250);
%             thbpower = yac(wac>=2.5&wac<=6);        % theta band: 2.5-6 Hz
%             allpower = yac(wac>=0.5&wac<=100);      % all power between 0.5 and 100 Hz
%             ThetaThetamod(t) = sum(thbpower) / sum(allpower);
%         else
%             ThetaFiringRateC1(t) = 0;
%             ThetaFiringRateC2(t) = 0;
%             ThetaThetamod(t) = NaN;
%         end
        cnt = th1;
        loceeg = eeg(cnt-wn:cnt+wn);
        subplot(3,1,1)
        hold on
        plot(loceeg)
        subplot(3,1,2)
        hold on
        locunit = unit(cnt-wn:cnt+wn);
        plot(locunit)
        locvdisc = (vdisc(vdisc>cnt-wn&vdisc<cnt+wn)) - cnt + wn;
        z(locvdisc) = 1;
    end
    z2 = reshape(z,2*wn/20,20);
    subplot(3,1,3)
    bar(1:20,sum(z2))
    title(titlestr)
    fn = [fname(1:end-4) '_TRANSITION'];
    saveas(gcf,fn);
    
%     ff = [nothdir 'NONTHETA_SEGMENTS_' cmps{1} '_' cmps{2}(1:2)];  % load non-theta segments
%     load(ff)
%     if ~isempty(NonThetaSegments)
%         NonThetaSegments = unitesegn(NonThetaSegments,sr);     % drop gaps < 0.5 s
%         NonThetaSegments = short_killer(NonThetaSegments);    % drop segments < 3 s
%         lenon = NonThetaSegments(2,:) - NonThetaSegments(1,:);  % find the longest non-theta segment
%         minx = find(lenon==max(lenon));
%         minx1 = NonThetaSegments(1,minx);
%         minx2 = NonThetaSegments(2,minx);
%         minx2 = min(minx2,length(eeg));
%         nnc = std(eeg(minx1:cnst:minx2));    % normalizing constant for STA
%     end
%     no_index = size(NonThetaSegments,2);
%     VdiscNonTheta = [];
%     NonThetaFiringRateC1 = zeros(1,no_index);
%     NonThetaFiringRateC2 = zeros(1,no_index);
%     NonThetaThetamod = zeros(1,th_index);
%     for t = 1:no_index      % non-theta segment cycle
%         no1 = NonThetaSegments(1,t);
%         no2 = NonThetaSegments(2,t);
%         no2 = min(no2,length(eeg));
%         vdisc_noth = vdisc(vdisc>no1&vdisc<no2) - no1;
%         VdiscNonTheta = [VdiscNonTheta vdisc(vdisc>no1&vdisc<no2)];
%         
%         if ~isempty(vdisc_noth)        % firing rate and theta-modulation
%             efflen = (vdisc_noth(end) - vdisc_noth(1)) / sr;
%             NonThetaFiringRateC1(t) = (length(vdisc_noth) - 1);
%             NonThetaFiringRateC2(t) = efflen;
%             ac = racorr(vdisc_noth/10000,500);    % autocorrelation
%             [yac wac] = b_fft2(ac,250);
%             thbpower = yac(wac>=2.5&wac<=6);        % theta band: 2.5-6 Hz
%             allpower = yac(wac>=0.5&wac<=100);      % all power between 0.5 and 100 Hz
%             NonThetaThetamod(t) = sum(thbpower) / sum(allpower);
%         else
%             NonThetaFiringRateC1(t) = 0;
%             NonThetaFiringRateC2(t) = 0;
%             NonThetaThetamod(t) = NaN;
%         end
%     end
%     n = length(ThetaHilbert);
%     ftm = sum(exp(1).^(i*ThetaHilbert)) / n;    % first trigonometric moment
%     ThetaAng = angle(ftm);   % mean angle
%     ThetaMvl = abs(ftm);     % mean resultant length
%     
%     VdiscTheta2 = round(VdiscTheta/cnst);   % downsample on 1000 Hz
%     VdiscNonTheta2 = round(VdiscNonTheta/cnst);
%     eeg2 = eeg(1:cnst:end);
%     wn = 2 * dsr;    % 2 sec. window
%     [StaTheta StaIndexTheta1 StaIndexTheta2 StaAmpTheta nn] = astanormmod(VdiscTheta2,eeg2,wn,tnc);    % STA
%     H = stafig(StaTheta,StaIndexTheta1,StaIndexTheta2,nn,wn,sr,titlestr);
%     fnsta = [fname(1:end-4) '_theta_STA'];
%     saveas(H,fnsta);
%     [StaNonTheta StaIndexNonTheta1 StaIndexNonTheta2 StaAmpNonTheta nn] = astanormmod(VdiscNonTheta2,eeg2,wn,nnc);
%     H = stafig(StaNonTheta,StaIndexNonTheta1,StaIndexNonTheta2,nn,wn,sr,titlestr);
%     fnsta = [fname(1:end-4) '_nontheta_STA'];
%     saveas(H,fnsta);
    
%     xlsout{o,1} = fname(1:end-4);   % Excel output
%     xlsout{o,2} = ThetaAng * 180 / pi;
%     xlsout{o,3} = ThetaMvl;
%     xlsout{o,4} = sum(ThetaFiringRateC1) / sum(ThetaFiringRateC2);
%     xlsout{o,5} = b_mean_nonnan(ThetaThetamod);
%     xlsout{o,6} = StaIndexTheta1;
%     xlsout{o,7} = StaIndexTheta2;
%     xlsout{o,8} = StaAmpTheta;
%     xlsout{o,9} = sum(NonThetaFiringRateC1) / sum(NonThetaFiringRateC2);
%     xlsout{o,10} = b_mean_nonnan(NonThetaThetamod);
%     xlsout{o,11} = StaIndexNonTheta1;
%     xlsout{o,12} = StaIndexNonTheta2;
%     xlsout{o,13} = StaAmpNonTheta;
    
    close all
    waitbar(o/sf)
end
% xlswrite('summary',xlsout)
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