function raphepreview2
%RAPHEPREVIEW2   Preliminary analysis for MRN stimulation project.
%   RAPHEPREVIEW2 calculates and saves the following outputs for raphe
%   stimulation baseline mat files:
%       ISI histogram
%       instantenous unit frequency
%       EEG and unit wavelet
%       theta/delta band power
%       firing rate
%       firing phase during theta segments greater than 3 s.
%       unit histogram centered to ripple onset/center/end

% Directories
global DATADIR
global DATAPATH
inpdir = [DATADIR 'raphe_matfiles\temp\discriminated\'];
inpdir2 = [DATADIR 'raphe_matfiles\temp\'];
thetadir = [DATAPATH 'Wavelet_raphe\theta_segments\'];
nothdir = [DATAPATH 'Wavelet_raphe\nontheta_segments\'];
sharpdir = [DATAPATH 'Wavelet_raphe\sharpwave_segments\'];
resdir = [DATAPATH 'Raphe\View\'];
mm = pwd;
cd(resdir)

% Filelist
[files files_short] = b_filelist(inpdir);
sf = length(files_short);

% Progress indicator
wb = waitbar(0,'Running RAPHEPREVIEW2...','Position',[360 250 275 50]);
global WB
WB(end+1) = wb;

% Main
sr = 10000;     % sampling rate
nqf = sr / 2;   % Nyquist frequency
fl = fir1(512,200/sr,'high');   % highpass FIR filter
for o = 2:sf
    fname = files(o).name;
    cmps = strread(fname,'%s','delimiter','_');
    titlestr = [cmps{1} ' ' cmps{2}];
    ff = [inpdir fname];
    load(ff)
    ff2 = [inpdir2 fname(1:end-6) '.mat'];
    load(ff2)
    if exist('data')
        unitr = data(:,2);
    else
        unitr = Unit.values';
    end
    unit = filtfilt(fl,1,unitr)';   % unit may come along with field!
    unit = unit(1:length(eeg));     % some EEGs are truncated!
    
    isi = diff(vdisc);
    instfrek = [];      % instantenous frequency
    lenu = length(eeg);
    for k = 1:length(vdisc)-1
        instfrek(vdisc(k):vdisc(k+1)) = 1 / isi(k);
    end
    instfrek(1:vdisc(1)-1) = 1 / vdisc(1);
    instfrek(vdisc(end):lenu) = 1 / (lenu - vdisc(end));
    instfrek2 = instfrek * sr;
    figure
    plot(instfrek2)
    title(titlestr)
    fn = [fname(1:end-4) '_INSTFREQ'];
    saveas(gcf,fn);
    
    [powunit,phaseunit,f] = unitwavelet(vdisc,lenu,sr);     % wavelet
    figure
    imagesc(powunit)
    b_rescaleaxis('Y',f)
    title(titlestr)
    fn = [fname(1:end-4) '_UNITWAVELET.jpg'];
    saveas(gcf,fn);
    [poweeg,phaseeeg,f] = eegwavelet(eeg(1:10:end),sr/10);
    figure
    imagesc(poweeg)
    b_rescaleaxis('Y',f)
    title(titlestr)
    fn = [fname(1:end-4) '_EEGWAVELET.jpg'];
    saveas(gcf,fn);
    
    seglen = 5 * sr;        % firing rate; 5 sec. long segments
    len = length(eeg);
    olp = 5;    % 80% overlapping windows
    ind1 = 1:seglen/olp:len-seglen+1;
    ind2 = ind1 + seglen - 1;
    lenr = length(ind1);
    frate = zeros(1,lenr);
    for k = 1:lenr
        vd = vdisc(vdisc>ind1(k)&vdisc<ind2(k)) - ind1(k);
        if ~isempty(vd)
            efflen = (vd(end) - vd(1)) / sr;
            frate(k) = length(vd) / efflen;
        else
            frate(k) = NaN;
        end
    end
    figure
    plot(frate)
    ylabel('firing rate')
    title(titlestr)
    fn = [fname(1:end-4) '_FR'];
    saveas(gcf,fn);
    
    ff = [thetadir 'THETA_SEGMENTS_' cmps{1} '_' cmps{2}(1:2)];  % phase; load theta segments
    load(ff)
    ThetaSegments = uniteseg(ThetaSegments,sr);     % drop gaps < 0.5 s
    thlen = ThetaSegments(2,:) - ThetaSegments(1,:);    % longest theta segment
    ml = find(thlen==max(thlen));
    flt = fir1(4096,[2/nqf 8/nqf]);  % filtering EEG
    seeg = (eeg - mean(eeg)) / std(eeg);
    feeg = filtfilt(flt,1,seeg);
    th1 = ThetaSegments(1,ml);
    th2 = ThetaSegments(2,ml);
    eeg_theta = eeg(th1:th2);
    feeg_theta = feeg(th1:th2);
    vdisc_theta = vdisc(vdisc>th1&vdisc<th2) - th1;
    ahee = angle(hilbert(feeg_theta));    % Hilbert-transformation of the EEG
    fn = find(-diff(ahee)>2*pi-0.1);
    sd = std(feeg);
    inx = find(vdisc_theta<fn(1));
    for k = 1:length(fn)-1
        seeg = feeg(fn(k):fn(k+1));
        axs = max(seeg) - min(seeg);
        sahee = ahee(fn(k):fn(k+1));
        if (axs < 2 * sd)  || any(diff(sahee(2:end))<0) || (fn(k+1) - fn(k) < 1 / 12 * sr)
            inx = [inx find(vdisc_theta>fn(k)&vdisc_theta<fn(k+1))];
        end    % discard segments > 12 Hz freq. | Hilbert-transform < mean + 2SD
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
    figure
    [nm xout] = hist(bahee,20);
    bar([xout xout+2*pi],[nm nm])
    x_lim = xlim;
    y_lim = ylim;
    text(x_lim(1)+(x_lim(2)-x_lim(1))/2,y_lim(1)+(y_lim(2)-y_lim(1))*0.9,...
        ['mean angle: ' num2str(hang)])
    text(x_lim(1)+(x_lim(2)-x_lim(1))/2,y_lim(1)+(y_lim(2)-y_lim(1))*0.75,...
        ['mean vector length: ' num2str(hmvl)])
    title(titlestr)
    fn = [fname(1:end-4) '_PHASE'];
    saveas(gcf,fn);
    
    isi = diff(vdisc);      % ISI histogram
    fn1 = [fname(1:end-4) '_ISIh'];
    fn2 = [fname(1:end-4) '_ISIh2'];
    isihist(isi,fn1,fn2,titlestr,sr)
    
    ff = [nothdir 'NONTHETA_SEGMENTS_' cmps{1} '_' cmps{2}(1:2)];  % load non-theta segments
    load(ff)
    NonThetaSegments = uniteseg(NonThetaSegments,sr);     % drop gaps < 0.5 s
    nthlen = NonThetaSegments(2,:) - NonThetaSegments(1,:);    % longest non-theta segment
    nml = find(nthlen==max(nthlen));
    nth1 = NonThetaSegments(1,nml);
    nth2 = NonThetaSegments(2,nml);
    isi_theta = diff(vdisc_theta);
    fn1 = [fname(1:end-4) '_ISIh_theta'];
    fn2 = [fname(1:end-4) '_ISIh2_theta'];
    isihist(isi_theta,fn1,fn2,[titlestr ' theta'],sr)
    vdisc_noth = vdisc(vdisc>nth1&vdisc<nth2) - nth1;
    isi_noth = diff(vdisc_noth);
    fn1 = [fname(1:end-4) '_ISIh_noth'];
    fn2 = [fname(1:end-4) '_ISIh2_noth'];
    isihist(isi_noth,fn1,fn2,[titlestr ' noth'],sr)
    close all
    
    czacorr2(vdisc_theta/1000,50)    % autocorrelation
    title([titlestr ' theta'])
    fn = [fname(1:end-4) '_ACORR_theta'];
    saveas(gcf,fn);
    czacorr2(vdisc_noth/1000,50)
    title([titlestr ' noth'])
    fn = [fname(1:end-4) '_ACORR_noth'];
    saveas(gcf,fn);
    
    ff = [sharpdir 'SHARPWAVE_SEGMENTS_' cmps{1} '_' ...
        cmps{2}(1:2)];  % sharpwaves; load sharpwave segments
    load(ff)
    flt = fir1(2048,[90/nqf 140/nqf]);
    feeg = filtfilt(flt,1,eeg);
    figure
    wn = 10000;
    z = zeros(1,2*wn);
    for k = 1:size(SharpWaveSegments,2)
        cnt = round(mean([SharpWaveSegments(1,k),SharpWaveSegments(2,k)]));
        loceeg = feeg(cnt-wn:cnt+wn);
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
    fn = [fname(1:end-4) '_RIPPLEb'];
    saveas(gcf,fn);
    figure
    wn = 10000;
    z = zeros(1,2*wn);
    for k = 1:size(SharpWaveSegments,2)
        cnt = SharpWaveSegments(1,k);
        loceeg = feeg(cnt-wn:cnt+wn);
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
    fn = [fname(1:end-4) '_RIPPLEa'];
    saveas(gcf,fn);
    figure
    wn = 10000;
    z = zeros(1,2*wn);
    for k = 1:size(SharpWaveSegments,2)
        cnt = SharpWaveSegments(2,k);
        loceeg = feeg(cnt-wn:cnt+wn);
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
    fn = [fname(1:end-4) '_RIPPLEc'];
    saveas(gcf,fn);
    
    close all
    waitbar(o/sf)
end
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

% -------------------------------------------------------------------------
function isihist(isi,fn1,fn2,str,sr)

figure
hist(isi,50)
title(str)
saveas(gcf,fn1);
figure
hist(isi(isi<1*sr),50)
title(str)
saveas(gcf,fn2);    