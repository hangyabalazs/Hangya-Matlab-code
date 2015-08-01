function cz2phase
%CZ2PHASE   Wavelet phase.
%   CZ2PHASE calculates theta phase of hippocampal interneurons using
%   wavelet and Hilbert methods. Theta selection is based on amplitude mean
%   + 2*SD threshold and time limitations of the filtered theta cycles in
%   the Hilbert method and theta band wavelet power mean + SD threshold
%   applied on the power values corresponding to time-wise maximum along
%   theta-band scales in the wavelet method. Wavelet phase is assessed as
%   the angle of the mean complex theta band wavelet coefficients.
%
%   See also CZPHASE.

% Directories
global DATAPATH
inpdir_unit = [DATAPATH 'Czurko2\units\'];
inpdir_EEG = [DATAPATH 'Czurko2\EEG\'];
EEGtable = [DATAPATH 'Czurko2\inter_anal_02d_L_2SD_depth_2jmp_DG_OSC2_2blzs_bh2b.xls'];
resdir = [DATAPATH 'Czurko2\Phase2\'];
files = b_filelist(inpdir_unit);
sf = length(files);
mm = pwd;
cd(resdir)

% Progress indicator
wb = waitbar(0,'Running CZ2PHASE...','Position',[360 250 315 50]);    %Progress indicator
global WB
WB(end+1) = wb;

% Load EEG table
[tbl0 tbl00 tbl] = xlsread(EEGtable,'OSC_2');

% Main: calculate phase distribution
ftms = zeros(1,sf);
ftms_rayleigh = [];
ftms_rao = [];
PRayleigh = [];
PRao = [];
wftms = zeros(1,sf);
wftms_rayleigh = [];
wftms_rao = [];
wPRayleigh = [];
wPRao = [];
O = [1:10 12:16 19:28 30:33 36:40 45 47 49:53 59:61 68:sf];
for o = O
    disp(o)
    fname = files(o).name
    fs = findstr(fname,'_');
    fn = fname(1:fs(1)-1);
    ff = [inpdir_unit fname];    % load unit
    [unit0 unit] = xlsread(ff);
    if isempty(str2num(unit{1}))
        unit(1,:)=[];
        unit = unit(:);
        iu = cellfun(@isempty,unit);
        unit(iu) = [];
    end
    unit = cellfun(@str2num,unit);
    
    inx = find(strcmp(cellfun(@eegname,{tbl{1:end,2}},'UniformOutput',0),fn));
    EEGfile = tbl{inx,19};  % load EEG
    ff = [inpdir_EEG fn '_' EEGfile '.mat'];
    load(ff)
    eval(['Eeg = ' EEGfile ';']);
    eval(['clear ' EEGfile]);
    eeg = (-1) * Eeg.values;        % Neuralynx reverses the polarity!
    sr = 1 / Eeg.interval;
    eeg_start = Eeg.start;
    eeg_end = eeg_start + (Eeg.length - 1) * Eeg.interval;
    eeg_times = eeg_start:Eeg.interval:eeg_end;
    
    sp = find(abs(eeg)>eps*10,1,'first');   % cut spikes with no corresponding EEG
    unit(unit<eeg_times(sp)) = [];
    
% Downsample EEG
    eeg2 = eeg(1:5:end);
    eeg_times2 = eeg_times(1:5:end);
    sr2 = sr / 5;
    vdisc = unit;   % in seconds
        
% Phase calculation
    [hang, hmvl, ftm, angs] = thetaphase(eeg2,vdisc,sr2,eeg_times2);
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
    H = figure;
    B = bar(cnts,nm'/length(angs));
    set(B,'FaceColor',[0.16 0.38 0.27])
    y_lim = ylim;
    axis([-200 200 y_lim(1) y_lim(2)])
    str = ['\it{Mean angle: }' '\bf ' num2str(hang*180/pi)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = ['\it{Mean resultant length: }' '\bf ' num2str(hmvl)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
    str = ['\it{n: }' '\bf ' num2str(length(angs))];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    ts = fn;
    title(ts)
    
% Save
    cd Hilbert2
    fns = [fname(1:end-4) '_PHASEHIST.fig'];
    saveas(H,fns)
    fns = [fname(1:end-4) '_STAT.mat'];
    save(fns,'p_rayleigh','p_rao')
    fns = [fname(1:end-4) '_ANGS.mat'];
    save(fns,'angs','hang','hmvl','ftm','p_rayleigh')
    cd ..
    
% Wavelet
    dsr = 200;
    const = round(sr/dsr);
    dsr = sr / const;
    eegw = eeg(1:const:end);     % downsampling EEG
    vdiscw = round(unit*dsr);
    len = length(eegw);
    [wave_eeg,f] = eegwavelet(eegw,dsr);        % EEG wavelet
    pow_eeg = abs(wave_eeg) .^ 2;
    
% Wavelet phase
    fnd = find(f>12);    % find theta band
    pwind1 = fnd(end);
    fnd = find(f<4);
    pwind2 = fnd(1);
    
    thp = pow_eeg(pwind1:pwind2,:);
    ththres = mean(thp(:)) + std(thp(:));     % theta threshold: mean + SD in theta band
        
    sw2 = size(pow_eeg,2);    % maximum localizations
    maxes = max(thp);         % theta segments: max. above the threshold
%     inxs = find(maxes>ththres);
%     ThetaSegments = thresholder(inxs);
%     inxs = find(maxes<=ththres);
%     NonThetaSegments = thresholder(inxs);
    istheta = zeros(1,len);
    seglen = 0.5 * dsr;
    lenr = floor(len/seglen);
    ind1 = 1:seglen:len;
    ind2 = ind1 + seglen - 1;
    for k = 1:lenr
        ist = mean(maxes(ind1(k):ind2(k))) > ththres;
        if ist
            istheta(ind1(k):ind2(k)) = 1;
        end
    end
    vdiscw(vdiscw>length(istheta)) = [];
    vdiscw(vdiscw<1) = [];
    thetaspikes = vdiscw(find(istheta(vdiscw)));
    nonthetaspikes = vdiscw(find(~istheta(vdiscw)));
    
%     maxloc = zeros(1,sw2);
%     for k = 1:sw2
%         maxloc(k) = find(pow_eeg(:,k)==maxes(k));
%     end
%     wph = angle(wave_eeg(maxloc,:));  % angles
    pwph = mean(wave_eeg(pwind1:pwind2,:));
    wph = angle(pwph);
%     wph = circular_mean(angle(wave_eeg(pwind1:pwind2,:)),'rad');  % angles
    wphh = wph(thetaspikes);
%     figure      % phase definition plot
%     plot(eegw*10000)
%     hold on
%     nqf = dsr / 2;
%     flt = fir1(1024,[4 12]/nqf,'band');      % bandpass filtering: theta frequency band bounderies: 8-10 Hz
%     feeg = filtfilt(flt,1,eegw*10000);
%     plot(feeg,'k')
%     fis = find(istheta);
%     cm = colormap('HSV');
%     for k = 1:length(fis)
%         ph = wph(fis(k));
%         cph = (ph - (-pi)) / (2 * pi) * 64;
%         cph = round(cph);
%         cph = max(cph,1);
%         cph = min(cph,64);
%         clr = cm(cph,:);
% %         plot(fis(k):fis(k+1),feeg(fis(k):fis(k+1)),'Color',clr,'Marker','.','MarkerSize',16)
%         plot(fis(k):fis(k+1),feeg(fis(k):fis(k+1)),'Color',clr,'LineWidth',3)
%     end
        
    n = length(wphh);
    wftm = sum(exp(1).^(i*wphh)) / n;    % first trigonometric moment
    wang = angle(wftm);   % mean angle
    wmvl = abs(wftm);     % mean resultant length
    wangs = wphh;
    
    wftms(o) = wftm;     % significance testing
    [Z,wp_rayleigh,U,wp_rao] = b_rao(wangs);
    if wp_rayleigh < 0.001
        wftms_rayleigh = [wftms_rayleigh wftm];
    end
    wPRayleigh = [wPRayleigh wp_rayleigh];
    if wp_rao(2) <= 0.001
        wftms_rao = [wftms_rao wftm];
    end
    wPRao = [wPRao; wp_rao];
    
    [nm,xout] = histc(wangs*180/pi,edges);   % phase histogram
    nm = nm(1:end-1);
    H = figure;
    B = bar(cnts,nm'/length(wangs));
    set(B,'FaceColor',[0.16 0.38 0.27])
    y_lim = ylim;
    axis([-200 200 y_lim(1) y_lim(2)])
    str = ['\it{Mean angle: }' '\bf ' num2str(wang*180/pi)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = ['\it{Mean resultant length: }' '\bf ' num2str(wmvl)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
    str = ['\it{n: }' '\bf ' num2str(length(wangs))];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    title(fn)
    
% Save
    cd wavelet2
    fns = [fname(1:end-4) '_WPHASEHIST.fig'];
    saveas(H,fns)
    fns = [fname(1:end-4) '_WSTAT.mat'];
    save(fns,'wp_rayleigh','wp_rao')
    fns = [fname(1:end-4) '_WANGS.mat'];
    save(fns,'wangs','wang','wmvl','wftm','wp_rayleigh')
    cd ..
    
    close all
    waitbar(o/sf)
end
fns = 'FTMS_pos_rayleigh.mat';
save(fns,'ftms_rayleigh')
fns = 'FTMS_pos_rao.mat';
save(fns,'ftms_rao')
fns = 'PVAL.mat';
save(fns,'PRayleigh','PRao')
fns = 'WFTMS_pos_rayleigh.mat';
save(fns,'wftms_rayleigh')
fns = 'WFTMS_pos_rao.mat';
save(fns,'wftms_rao')
fns = 'WPVAL.mat';
save(fns,'wPRayleigh','wPRao')

close(wb)
cd(mm)

% -------------------------------------------------------------------------
function A = eegname(S)

fs = findstr(S,'/');
if isempty(fs)
    A = '';
    return
end
A = S(fs(end)+1:end);

% -------------------------------------------------------------------------
function [hang, hmvl, ftm, bahee] = thetaphase(eeg,vdisc,sr,eeg_times)

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
for k = 1:length(fn)-1
    seeg = feeg(fn(k):fn(k+1));
    axs = max(seeg) - min(seeg);
    sahee = ahee(fn(k):fn(k+1));
    stim = eeg_times(fn(k):fn(k+1));
    if (axs < 2 * sd)  || (fn(k+1) - fn(k) < 0.083 * sr) || (fn(k+1) - fn(k) > 0.25 * sr)
        inx = [inx; find(vdisc>stim(1)&vdisc<stim(end))];
    end
end
vdisc(inx) = [];

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

% -------------------------------------------------------------------------
function th = thresholder(ip)

drml = diff(ip);
fdr = find(drml>1);
lenfdr = length(fdr);
preth = zeros(2,lenfdr+1);
preth(1,1) = ip(1);
for t = 1:lenfdr
    preth(2,t) = ip(fdr(t));
    preth(1,t+1) = ip(fdr(t)+1);
end
preth(2,end) = ip(end);
th = preth;

% ----------------------------------------------------------------------------------
function segments = short_killer(segments,sr)

% Skip short segments
int = segments;
int1 = int(1,:);
int2 = int(2,:);
difint = int2 - int1;
fd = find(difint<0.5*sr);         % leaving segments shorter than 0.5 sec.
int1(fd) = [];
int2(fd) = [];
segments = [int1; int2];