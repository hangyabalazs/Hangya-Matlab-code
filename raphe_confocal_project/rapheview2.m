function rapheview2
%RAPHEVIEW   Analysis for MRN stimulation project.
%   RAPHEVIEW calculates and saves the following outputs for raphe
%   stimulation baseline mat files:
%       ISI histogram
%       instantenous unit frequency
%       EEG and unit wavelet
%       theta/delta band power
%       unit autocorrelation
%       firing rate
%       firing phase during theta segments greater than 3 s.
%       unit histogram centered to ripple onset/center/end
%
%   See also RAPHEPREVIEW.

% Directories
global DATADIR
global DATAPATH
inpdir = [DATADIR 'raphe_matfiles\temp\discriminated\'];
inpdir2 = [DATADIR 'raphe_matfiles\temp\'];
thetadir = [DATAPATH 'Wavelet_raphe\theta_segments\'];
nothdir = [DATAPATH 'Wavelet_raphe\nontheta_segments\'];
sharpdir = [DATAPATH 'Wavelet_raphe\sharpwave_segments\'];
resdir = [DATAPATH 'Raphe\View2\'];
mm = pwd;
cd(resdir)

% Filelist
[files files_short] = b_filelist(inpdir);
sf = length(files_short);

% Progress indicator
wb = waitbar(0,'Running RAPHEVIEW...','Position',[360 250 275 50]);
global WB
WB(end+1) = wb;

% Main
sr = 10000;     % sampling rate
nqf = sr / 2;   % Nyquist frequency
xlsout = cell(sf,7);    % initialize Excel output
for o = 1:sf
    fname = files(o).name;
    cmps = strread(fname,'%s','delimiter','_');
    titlestr = [cmps{1} ' ' cmps{2}];
    ff = [inpdir fname];
    load(ff)
    
    ff = [thetadir 'THETA_SEGMENTS_' cmps{1} '_' cmps{2}(1:2)];  % phase; load theta segments
    load(ff)
    ThetaSegments = uniteseg(ThetaSegments,sr);     % drop gaps < 0.5 s
    ThetaSegments = short_killer(ThetaSegments);    % drop segments < 3 s
    th_index = size(ThetaSegments,2);
%     thlen = ThetaSegments(2,:) - ThetaSegments(1,:);    % longest theta segment
%     ml = find(thlen==max(thlen));
    flt = fir1(4096,[2/nqf 8/nqf]);  % filtering EEG
    seeg = (eeg - mean(eeg)) / std(eeg);
    feeg = filtfilt(flt,1,seeg);
    ThetaHilbert = [];
    ThetaFiringRateC1 = zeros(1,th_index);
    ThetaFiringRateC2 = zeros(1,th_index);
    ThetaThetamod = zeros(1,th_index);
    for t = 1:th_index      % theta segment cycle
        th1 = ThetaSegments(1,t);
        th2 = ThetaSegments(2,t);
        th2 = min(th2,length(eeg));
        eeg_theta = eeg(th1:th2);
        feeg_theta = feeg(th1:th2);
        vdisc_theta = vdisc(vdisc>th1&vdisc<th2) - th1;
        ahee = angle(hilbert(feeg_theta));    % Hilbert-transformation of the EEG
        fn = find(-diff(ahee)>2*pi-0.3);
        sd = std(feeg_theta);
        inx = find(vdisc_theta<fn(1));
        for k = 1:length(fn)-1
            seeg = feeg_theta(fn(k):fn(k+1));
            axs = max(seeg) - min(seeg);
            sahee = ahee(fn(k):fn(k+1));
            if (axs < 1 * sd)  || any(diff(sahee(2:end))<0) || (fn(k+1) - fn(k) < 1 / 8 * sr)
                inx = [inx find(vdisc_theta>fn(k)&vdisc_theta<fn(k+1))];
            end    % discard segments > 8 Hz freq. | Hilbert-transform < mean + 2SD
        end
        inx = [inx find(vdisc_theta>fn(end))];
        vdisc_theta2 = vdisc_theta;
        vdisc_theta2(inx) = [];
        bahee = ahee(vdisc_theta2(vdisc_theta2>0&vdisc_theta2<...
            length(eeg_theta)));        % phase angles - Hilbert
        n = length(bahee);
        ftm = sum(exp(1).^(i*bahee)) / n;    % first trigonometric moment
        hang = angle(ftm);   % mean angle
        hmvl = abs(ftm);     % mean resultant length
        z = n * (hmvl ^ 2);  % Rayleigh's Z statistic
        p = exp(1) ^ (-1 * z) * (1 + (2 * z - z ^ 2) / ...
            (4 * n) - (24 * z - 132 * z ^ 2 + 76 * z ^ 3 - 9 * z ^ 4) / (288 * n ^ 2));
        if p < 0.001
            ThetaHilbert = [ThetaHilbert bahee];
        end
        
        if ~isempty(vdisc_theta)        % firing rate and theta-modulation
            efflen = (vdisc_theta(end) - vdisc_theta(1)) / sr;
            ThetaFiringRateC1(t) = (length(vdisc_theta) - 1)
            ThetaFiringRateC2(t) = efflen;
            ac = racorr(vdisc_theta/10000,500);    % autocorrelation
            [yac wac] = b_fft2(ac,250);
            thbpower = yac(wac>=2.5&wac<=6);        % theta band: 2.5-6 Hz
            allpower = yac(wac>=0.5&wac<=100);      % all power between 0.5 and 100 Hz
            ThetaThetamod(t) = sum(thbpower) / sum(allpower);
        else
            ThetaFiringRateC1(t) = 0;
            ThetaFiringRateC2(t) = 0;
            ThetaThetamod(t) = NaN;
        end
    end
    ff = [nothdir 'NONTHETA_SEGMENTS_' cmps{1} '_' cmps{2}(1:2)];  % load non-theta segments
    load(ff)
    NonThetaSegments = uniteseg(NonThetaSegments,sr);     % drop gaps < 0.5 s
    NonThetaSegments = short_killer(NonThetaSegments);    % drop segments < 3 s
    no_index = size(NonThetaSegments,2);
%     nthlen = NonThetaSegments(2,:) - NonThetaSegments(1,:);    % longest non-theta segment
%     nml = find(nthlen==max(nthlen));
    NonThetaFiringRateC1 = zeros(1,no_index);
    NonThetaFiringRateC2 = zeros(1,no_index);
    NonThetaThetamod = zeros(1,th_index);
    for t = 1:no_index      % non-theta segment cycle
        no1 = NonThetaSegments(1,t);
        no2 = NonThetaSegments(2,t);
        no2 = min(no2,length(eeg));
        vdisc_noth = vdisc(vdisc>no1&vdisc<no2) - no1;
        
        if ~isempty(vdisc_noth)        % firing rate and theta-modulation
            efflen = (vdisc_noth(end) - vdisc_noth(1)) / sr;
            NonThetaFiringRateC1(t) = (length(vdisc_noth) - 1);
            NonThetaFiringRateC2(t) = efflen;
            ac = racorr(vdisc_noth/10000,500);    % autocorrelation
            [yac wac] = b_fft2(ac,250);
            thbpower = yac(wac>=2.5&wac<=6);        % theta band: 2.5-6 Hz
            allpower = yac(wac>=0.5&wac<=100);      % all power between 0.5 and 100 Hz
            NonThetaThetamod(t) = sum(thbpower) / sum(allpower);
        else
            NonThetaFiringRateC1(t) = 0;
            NonThetaFiringRateC2(t) = 0;
            NonThetaThetamod(t) = NaN;
        end
    end
    n = length(ThetaHilbert);
    ftm = sum(exp(1).^(i*ThetaHilbert)) / n;    % first trigonometric moment
    ThetaAng = angle(ftm);   % mean angle
    ThetaMvl = abs(ftm);     % mean resultant length
    
    xlsout{o,1} = fname(1:end-4);   % Excel output
    xlsout{o,2} = ThetaAng * 180 / pi;
    xlsout{o,3} = ThetaMvl;
    xlsout{o,4} = sum(ThetaFiringRateC1) / sum(ThetaFiringRateC2);
    xlsout{o,5} = b_mean_nonnan(ThetaThetamod);
    xlsout{o,6} = sum(NonThetaFiringRateC1) / sum(NonThetaFiringRateC2);
    xlsout{o,7} = b_mean_nonnan(NonThetaThetamod);
    
%     ff = [sharpdir 'SHARPWAVE_SEGMENTS_' cmps{1} '_' ...
%         cmps{2}(1:2)];  % sharpwaves; load sharpwave segments
%     load(ff)
%     flt = fir1(2048,[90/nqf 140/nqf]);
%     feeg = filtfilt(flt,1,eeg);
%     figure
%     wn = 10000;
%     z = zeros(1,2*wn);
%     for k = 1:size(SharpWaveSegments,2)
%         cnt = round(mean([SharpWaveSegments(1,k),SharpWaveSegments(2,k)]));
%         loceeg = feeg(cnt-wn:cnt+wn);
%         subplot(3,1,1)
%         hold on
%         plot(loceeg)
%         subplot(3,1,2)
%         hold on
%         locunit = unit(cnt-wn:cnt+wn);
%         plot(locunit)
%         locvdisc = (vdisc(vdisc>cnt-wn&vdisc<cnt+wn)) - cnt + wn;
%         z(locvdisc) = 1;
%     end
%     z2 = reshape(z,2*wn/20,20);
%     subplot(3,1,3)
%     bar(1:20,sum(z2))
%     title(titlestr)
%     fn = [fname(1:end-4) '_RIPPLEb'];
%     saveas(gcf,fn);
%     figure
%     wn = 10000;
%     z = zeros(1,2*wn);
%     for k = 1:size(SharpWaveSegments,2)
%         cnt = SharpWaveSegments(1,k);
%         loceeg = feeg(cnt-wn:cnt+wn);
%         subplot(3,1,1)
%         hold on
%         plot(loceeg)
%         subplot(3,1,2)
%         hold on
%         locunit = unit(cnt-wn:cnt+wn);
%         plot(locunit)
%         locvdisc = (vdisc(vdisc>cnt-wn&vdisc<cnt+wn)) - cnt + wn;
%         z(locvdisc) = 1;
%     end
%     z2 = reshape(z,2*wn/20,20);
%     subplot(3,1,3)
%     bar(1:20,sum(z2))
%     title(titlestr)
%     fn = [fname(1:end-4) '_RIPPLEa'];
%     saveas(gcf,fn);
%     figure
%     wn = 10000;
%     z = zeros(1,2*wn);
%     for k = 1:size(SharpWaveSegments,2)
%         cnt = SharpWaveSegments(2,k);
%         loceeg = feeg(cnt-wn:cnt+wn);
%         subplot(3,1,1)
%         hold on
%         plot(loceeg)
%         subplot(3,1,2)
%         hold on
%         locunit = unit(cnt-wn:cnt+wn);
%         plot(locunit)
%         locvdisc = (vdisc(vdisc>cnt-wn&vdisc<cnt+wn)) - cnt + wn;
%         z(locvdisc) = 1;
%     end
%     z2 = reshape(z,2*wn/20,20);
%     subplot(3,1,3)
%     bar(1:20,sum(z2))
%     title(titlestr)
%     fn = [fname(1:end-4) '_RIPPLEc'];
%     saveas(gcf,fn);
    
    close all
    waitbar(o/sf)
end
xlswrite('summary',xlsout)
close(wb)
cd(mm)

% -------------------------------------------------------------------------
function segments2 = uniteseg(segments,sr)

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