function czphase_resubmission
%CZPHASE   Theta phase for hippocampal interneurons.
%   CZPHASE plots and saves theta phase histograms for hippocampal
%   interneurons. Theta phase is calculated via Hilbert-transform of the EEG
%   filtered between 4 and 12 Hz. First trigonometric moments are also
%   saved for all cells.
%
%   See also CZRIPPLE.

% Directories
global DATAPATH
inpdir_eeg = [DATAPATH 'Czurko\czurko_EEG\'];
inpdir_unit = [DATAPATH 'Czurko\discriminated2\'];
resdir = [DATAPATH 'Czurko\czurko_EEG\phase\resubmission\'];
mm = pwd;
cd(resdir)

% Load
xlsname = [inpdir_eeg 'EEG3.xls'];
[ntz mtz] = xlsread(xlsname,'all_unique');
sf = size(mtz,1);   % number of pairs
ftms = zeros(1,sf);
ftms_rayleigh = [];
ftms_rao = [];
PRayleigh = [];
PRao = [];
for o = 1:1
    disp(o)
    dr = [inpdir_unit mtz{o,4} '\' mtz{o,1} '\'];
    fn = [dr mtz{o,2} '.mat'];
    load(fn)        % load unit
    vdisc = data;
    
    fn = [inpdir_eeg 'EEG_' mtz{o,3} '_' mtz{o,1} '.mat'];
    load(fn)        % load EEG
    eval(['Eeg = ' mtz{o,3} ';']);
    eval(['clear ' mtz{o,3}]);
    eeg = Eeg.values;
    sr = 1 / Eeg.interval;
    eeg_start = Eeg.start;
%     eeg_start = 0;
    eeg_end = eeg_start + (Eeg.length - 1) * Eeg.interval;
    eeg_times = eeg_start:Eeg.interval:eeg_end;
    
% Downsample EEG
    eeg2 = eeg(1:5:end);
    eeg_times2 = eeg_times(1:5:end);
    sr2 = sr / 5;
        
% Phase calculation
    [hang hmvl ftm angs isexcluded exclratio exclline excltime excllinetime] = ...
        thetaphase(eeg2,vdisc,sr2,eeg_times2);
    ftms(o) = ftm;
    [Z,p_rayleigh,U,p_rao] = b_rao(angs);
    if p_rayleigh < 0.001
        ftms_rayleigh = [ftms_rayleigh ftm];
    end
    PRayleigh = [PRayleigh p_rayleigh];
    if p_rao(2) <= 0.001
        ftms_rao = [ftms_rao ftm];
    end
    PRao = [PRao; p_rao];
    
    edges = -180:20:180;     % edges for phase histogram
    cnts = (edges(1:end-1) + edges(2:end)) / 2;
    [nm,xout] = histc(angs*180/pi,edges);   % phase histogram
    nm = nm(1:end-1);
%     H = figure;
%     B = bar(cnts,nm'/length(angs));
%     set(B,'FaceColor',[0.16 0.38 0.27])
%     y_lim = ylim;
%     axis([-200 200 y_lim(1) y_lim(2)])
%     str = ['\it{Mean angle: }' '\bf ' num2str(hang*180/pi)];
%     text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
%     str = ['\it{Mean resultant length: }' '\bf ' num2str(hmvl)];
%     text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
%     str = ['\it{n: }' '\bf ' num2str(length(angs))];
%     text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    ts = [mtz{o,1} ' ' mtz{o,2}];
    ts(ts=='_') = ' ';
%     title(ts)
    fns = [ts '_exclude_ratio.mat'];
%     save(fns,'exclratio','exclline')
    
% Wavelet
    dsr = 200;
    const = round(sr/dsr);
    dsr = sr / const;
    eegw = eeg(1:const:end);     % downsampling EEG
    len = length(eegw);
    [wave_eeg,f] = eegwavelet(eegw,dsr);        % EEG wavelet
    im_eeg = abs(wave_eeg) .^ 2;
    clear wave_eeg
    
    Hw = figure;
    S1 = subplot('Position',[0.13 0.25+0.1 0.775 0.70-0.1]);
    imagesc(linspace(0,length(eeg2)/sr2/60,size(im_eeg,2)),0:size(im_eeg,1),im_eeg)
    set(gca,'CLim',[0 70])
    b_rescaleaxis('Y',f)
    S2 = subplot('Position',[0.13 0.09 0.775 0.18]);
    plot(excllinetime/sr2/60,exclline)
    xlim([0 length(eeg2)/sr2/60])
    linkaxes([S1,S2],'x')
    
% Save
    fns = [ts '_WAVELET.fig'];
    saveas(Hw,fns)
    close(Hw)
end
cd(mm)



% -------------------------------------------------------------------------
function [hang, hmvl, ftm, bahee isexcluded exclratio exclline excltime excllinetime] = ...
    thetaphase(eeg,vdisc,sr,eeg_times)

% Filtering EEG
nqf = sr / 2;
flt = fir1(4096,[4 12]/nqf,'band');      % bandpass filtering
feeg = filtfilt(flt,1,eeg);

% Hilbert transformation of the EEG
ahee = angle(hilbert(feeg));
aheedeg = ahee * (180 / pi);

% Check criteria:
% 1. discard cicles with EEG amp. lower then 2SD
% 2. discard cicles shorter then 83 ms
% 3. discard cicles longer then 250 ms
fn = find(-diff(ahee)>2*pi-0.3);
sd = std(feeg);
inx = [];
isexcluded = zeros(1,length(fn)-1);
excltime = zeros(1,length(fn)-1);
for k = 1:length(fn)-1
    seeg = feeg(fn(k):fn(k+1));
    axs = max(seeg) - min(seeg);
    sahee = ahee(fn(k):fn(k+1));
    stim = eeg_times(fn(k):fn(k+1));
    if (axs < 2 * sd)  || (fn(k+1) - fn(k) < 0.083 * sr) || (fn(k+1) - fn(k) > 0.25 * sr)
        inx = [inx; find(vdisc>stim(1)&vdisc<stim(end))];
        isexcluded(k) = 1;
    end
    excltime(k) = fn(k) + (fn(k+1) - fn(k)) / 2;
end
vdisc(inx) = [];
exclratio = length(find(isexcluded)) / length(isexcluded);
exclline = zeros(1,length(fn)-50);
excllinetime = zeros(1,length(fn)-50);
for k = 50:length(fn)-1
    exclline(k-49) = length(find(isexcluded(k-49:k))) / 50;
    excllinetime(k-49) = mean(excltime(k-49:k));
end

% Phase angles - Hilbert
vdisc(vdisc<eeg_times(1)|vdisc>eeg_times(end)) = [];
lvd = length(vdisc);
bahee = zeros(1,lvd);
for k = 1:lvd    % interpolate
    inx1 = find(eeg_times<vdisc(k),1,'last');
    inx2 = inx1 + 1;
    rto = (vdisc(k) - eeg_times(inx1)) / (eeg_times(inx2) - eeg_times(inx1));
    if rto > 1
        error('Technical error 55!')
    end
    bahee(k) = ahee(inx1) + (ahee(inx2) - ahee(inx1)) * rto;
end
n = length(bahee);
ftm = sum(exp(1).^(i*bahee)) / n;    % first trigonometric moment
hang = angle(ftm);   % mean angle
hmvl = abs(ftm);     % mean resultant length



% -------------------------------------------------------------------------
function [wave,f] = eegwavelet(dat,sr)

% Prepare for wavelet transformation
variance = std(dat) ^ 2;
dat = (dat - mean(dat)) / sqrt(variance) ;
n = length(dat);
dt = 1 / sr;
pad = 1;
dj = 0.02;    
j1 = ceil((1/dj) * log2(n/2));
j1 = ceil(j1);
j = (0:j1);
s0 = 2 * dt; 
s = s0 .* 2 .^ (j * dj);
omega0 = 6;
c = 4 * pi / (omega0 + sqrt(2+omega0^2));
fperiod = c .* s;
f = 1 ./ fperiod;
lag1 = 0.72;
param = -1;
mif = 0.5;          %minimal intresting frequency
mis = find(f>mif);
mis = mis(end);     %maximal intristing scale
mother = 'Morlet';

% Wavelet transformation
[wave,period,scale,coi] = b_wavelet_new3(dat,dt,pad,dj,s0,j1,mother,param,mis);