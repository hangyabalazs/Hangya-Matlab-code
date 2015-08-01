function raphepreview
%RAPHEPREVIEW   Preliminary analysis for MRN stimulation project.
%   RAPHEPREVIEW calculates and saves the following outputs for raphe
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
sharpdir = [DATAPATH 'Wavelet_raphe\sharpwave_segments\'];
resdir = [DATAPATH 'Raphe\Preview\'];
mm = pwd;
cd(resdir)

% Filelist
[files files_short] = b_filelist(inpdir);
sf = length(files_short);

% Progress indicator
wb = waitbar(0,'Running RAPHEPREVIEW...','Position',[360 250 275 50]);
global WB
WB(end+1) = wb;

% Main
sr = 10000;     % sampling rate
nqf = sr / 2;   % Nyquist frequency
fl = fir1(512,200/nqf,'high');   % highpass FIR filter
for o = 1:sf
    fname = files(o).name;
    cmps = strread(fname,'%s','delimiter','_');
    titlestr = [cmps{1} ' ' cmps{2}];
    ff = [inpdir fname];
    load(ff)
    ff2 = [inpdir2 fname(1:end-6) '.mat'];
    load(ff2)
    unitr = data(:,2);
    unit = filtfilt(fl,1,unitr)';   % unit may come along with field!
    unit = unit(1:length(eeg));     % some EEGs are truncated!
    
    isi = diff(vdisc);      % ISI histogram
    figure
    hist(isi,50)
    fn = [fname(1:end-4) '_ISIh'];
    saveas(gcf,fn);
    figure
    hist(isi(isi<1*sr),50)
    fn = [fname(1:end-4) '_ISIh2'];
    saveas(gcf,fn);
    
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
    
%     [powunit,phaseunit,f] = unitwavelet(vdisc,lenu,sr);     % wavelet
%     figure
%     imagesc(powunit)
%     b_rescaleaxis('Y',f)
%     title(titlestr)
%     fn = [fname(1:end-4) '_UNITWAVELET.jpg'];
%     saveas(gcf,fn);
%     [poweeg,phaseeeg,f] = eegwavelet(eeg(1:10:end),sr/10);
%     figure
%     imagesc(poweeg)
%     b_rescaleaxis('Y',f)
%     title(titlestr)
%     fn = [fname(1:end-4) '_EEGWAVELET.jpg'];
%     saveas(gcf,fn);
%     
%     pwind1 = find(f<6,1,'first');       % theta/delta power; theta band
%     pwind2 = find(f>2.5,1,'last');
%     pwind3 = find(f<2.5,1,'first');      % delta band
%     pwind4 = find(f>0.5,1,'last');
%     tpd = sum(powunit(pwind1:pwind2,:)) ./ sum(powunit(pwind3:pwind4,:));
%     figure
%     plot(tpd)
%     title(titlestr)
%     fn = [fname(1:end-4) '_TPD'];
%     saveas(gcf,fn);
    
    seglen = 5 * sr;        % firing rate; 5 sec. long segments
    len = length(eeg);
    olp = 5;    % 80% overlapping windows
    ind1 = 1:seglen/olp:len-seglen+1;
    ind2 = ind1 + seglen - 1;
    lenr = length(ind1);
    frate = zeros(1,lenr);
    for k = 1:lenr
        vd = vdisc(vdisc>ind1(k)&vdisc<ind2(k)) - ind1(k);
        efflen = (vd(end) - vd(1)) / sr;
        frate(k) = length(vd) / efflen;
    end
    figure
    plot(frate)
    ylabel('firing rate')
    title(titlestr)
    fn = [fname(1:end-4) '_FR'];
    saveas(gcf,fn);
    
    ff = [thetadir 'THETA_SEGMENTS_' cmps{1} '_' cmps{2}(1:2)];  % phase; load theta segments
    load(ff)
    flt = fir1(512,[2/nqf 8/nqf]);  % filtering EEG
    feeg = filtfilt(flt,1,eeg);
    bahee = [];
    for t = 1:size(ThetaSegments,2)
        th1 = ThetaSegments(1,t);
        th2 = ThetaSegments(2,t);
        if (th2-th1)/sr < 3
            continue
        end
        eeg_theta = eeg(th1:th2);
        feeg_theta = feeg(th1:th2);
        vdisc_theta = vdisc(vdisc>th1&vdisc<th2)-th1;
        ahee = angle(hilbert(feeg_theta));    % Hilbert-transformation of the EEG
        angs = ahee(vdisc_theta(vdisc_theta>0&vdisc_theta<...
            length(eeg_theta)));        % phase angles - Hilbert
        bahee = [bahee angs];
    end
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
    
    th1 = 60 * sr;
    th2 = 90 * sr;
    eeg_theta = eeg(th1:th2);
    feeg_theta = feeg(th1:th2);
    vdisc_theta = vdisc(vdisc>th1&vdisc<th2)-th1;
    ahee = angle(hilbert(feeg_theta));    % Hilbert-transformation of the EEG
    bahee = ahee(vdisc_theta(vdisc_theta>0&vdisc_theta<...
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
    fn = [fname(1:end-4) '_PHASEtp'];
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