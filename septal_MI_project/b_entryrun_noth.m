function b_entryrun_noth
%ENTRYRUN_NOTH   Runs entropy and PDC on a sequence of non-theta segments.
%   ENTRYRUN_NOTH calculates wavelet power and wavelet phase entropy for
%   non-theta segments using 1000-point non-overlapping and 50% overlapping
%   time windows. It save line and time-frequency image visualization and
%   entropy matrices as well.
%
%   It also calculates Partial Directed Coherence and saves PDC images, 
%   lines and PDC matrices.
%
%   See also DTF, 2DTF, ENTRY, ENTRYRUN3 and ENTRYRUN_THETA.

% Input argument check
error(nargchk(0,0,nargin))

% Directories
global DATAPATH
global DATADIR2
inpdir = [DATAPATH 'Burst\Cluster\Noth\'];   % input directory
resdir = [DATAPATH 'Entry\Noth\'];
datadir = DATADIR2;

mm = pwd;
cd(resdir)
create_subdir

% Import
files = b_filelist(inpdir);
sf = length(files);
[datalist, dlist_short] = b_filelist(datadir);

% Progress indicator
wb = waitbar(0,'Running ENTRYRUN...','Position',[360 250 275 50]);
global WB
WB(end+1) = wb;

% Load
for o = 1:sf
    fname = files(o).name;
    ffnm = [inpdir fname];
    filenam = fname(1:6);
    fnm = fname(1:end-12);
    pont = findstr(fname,'.');
    fn_noext = fname(1:pont(1)-1);
    load(ffnm);     % load burst analysis results
    
    cmps = strread(fname,'%s','delimiter','_');
    seglen = str2num(cmps{5}) - str2num(cmps{4});
    i_first = str2num(cmps{4});
    i_second = str2num(cmps{5});
    global FILENAME
    FILENAME = fnm;
    global DATINX1
    DATINX1 = i_first;
    global DATINX2
    DATINX2 = i_second;
    
    inx = find(strcmp(filenam,dlist_short));
    fn = datalist(inx).name;
    ffn = [datadir fn];
    load(ffn);      % load raw data
    
% Create random unit
    lenu = length(eeg);
    r = random('exp',1000,1,10000);
    s = cumsum(r);
    psvd = unique(ceil(s));     % 'pseudo vdisc'
    psvd = psvd(find(psvd<lenu));
    
% Downsample & wavelet transformation of eeg
    [wavea_abs,wavea_phase,f] = eeg_wavelet(eeg(i_first:10:i_second)); % eeg (downsample on 1000 Hz)

% Sinc convolution & wavelet transformation of unit
    [waveb_abs,waveb_phase,f] = unit_wavelet(vdisc(find(vdisc>=i_first&vdisc<=i_second))-i_first,seglen);   % unit
    [wavec_abs,wavec_phase,f] = unit_wavelet(psvd(find(psvd>=i_first&psvd<=i_second))-i_first,seglen);    % random unit
    
    if size(wavea_abs,2) > size(waveb_abs,2)    % adjust size
        wavea_abs = wavea_abs(:,1:end-1);
        wavea_phase = wavea_phase(:,1:end-1);
    elseif size(waveb_abs,2) > size(wavea_abs,2)
        waveb_abs = waveb_abs(:,1:end-1);
        waveb_phase = waveb_phase(:,1:end-1);
        wavec_abs = wavec_abs(:,1:end-1);
        wavec_phase = wavec_phase(:,1:end-1);
    end

% Entropy
    fnd = find(f>6);    % frequency band bounderies
    pwind1 = fnd(end);
    fnd = find(f<2.5);
    pwind2 = fnd(1);
    
    abs1 = wavea_abs(pwind1:pwind2,:);  % eeg
    phase1 = wavea_phase(pwind1:pwind2,:);
    abs2 = waveb_abs(pwind1:pwind2,:);  % unit
    phase2 = waveb_phase(pwind1:pwind2,:);
    abs3 = wavec_abs(pwind1:pwind2,:);  % random unit
    phase3 = wavec_phase(pwind1:pwind2,:);
    
    cd([resdir 'entropy\power\line\windowsize1000_overlap1\real\'])
    entropy_line(abs1,abs2,f,1000,1)
    cd([resdir 'entropy\phase\line\windowsize1000_overlap1\real\'])
    entropy_line(phase1,phase2,f,1000,1)
    cd([resdir 'entropy\power\line\windowsize1000_overlap1\control\'])
    entropy_line(abs1,abs3,f,1000,1)
    cd([resdir 'entropy\phase\line\windowsize1000_overlap1\control\'])
    entropy_line(phase1,phase3,f,1000,1)
    cd([resdir 'entropy\power\line\windowsize1000_overlap2\real\'])
    entropy_line(abs1,abs2,f,1000,2)
    cd([resdir 'entropy\phase\line\windowsize1000_overlap2\real\'])
    entropy_line(phase1,phase2,f,1000,2)
    cd([resdir 'entropy\power\line\windowsize1000_overlap2\control\'])
    entropy_line(abs1,abs3,f,1000,2)
    cd([resdir 'entropy\phase\line\windowsize1000_overlap2\control\'])
    entropy_line(phase1,phase3,f,1000,2)
    
    fnd = find(f>6);    % frequency band bounderies
    pwind1 = fnd(end);
    fnd = find(f<2.5);
    pwind2 = fnd(1);
    
    abs1 = wavea_abs(pwind1:pwind2,:);  % eeg
    phase1 = wavea_phase(pwind1:pwind2,:);
    abs2 = waveb_abs(pwind1:pwind2,:);  % unit
    phase2 = waveb_phase(pwind1:pwind2,:);
    abs3 = wavec_abs(pwind1:pwind2,:);  % random unit
    phase3 = wavec_phase(pwind1:pwind2,:);
    
    cd([resdir 'entropy\power\image\windowsize1000_overlap1\real\'])
    entropy_image(abs1,abs2,pwind1,pwind2,f,1000,1)
    cd([resdir 'entropy\phase\image\windowsize1000_overlap1\real\'])
    entropy_image(phase1,phase2,pwind1,pwind2,f,1000,1)
    cd([resdir 'entropy\power\image\windowsize1000_overlap1\control\'])
    entropy_image(abs1,abs3,pwind1,pwind2,f,1000,1)
    cd([resdir 'entropy\phase\image\windowsize1000_overlap1\control\'])
    entropy_image(phase1,phase3,pwind1,pwind2,f,1000,1)
    cd([resdir 'entropy\power\image\windowsize1000_overlap2\real\'])
    entropy_image(abs1,abs2,pwind1,pwind2,f,1000,2)
    cd([resdir 'entropy\phase\image\windowsize1000_overlap2\real\'])
    entropy_image(phase1,phase2,pwind1,pwind2,f,1000,2)
    cd([resdir 'entropy\power\image\windowsize1000_overlap2\control\'])
    entropy_image(abs1,abs3,pwind1,pwind2,f,1000,2)
    cd([resdir 'entropy\phase\image\windowsize1000_overlap2\control\'])
    entropy_image(phase1,phase3,pwind1,pwind2,f,1000,2)
    cd(resdir)
    
    clear wavea_abs wavea_phase waveb_abs waveb_phase wavec_abs wavec_phase r s
    
% Directed Transfer Function
    cd([resdir 'dtf\halfsec\real\'])
    dtf(eeg,vdisc,0.5)
    cd([resdir 'dtf\halfsec\control\'])
    dtf(eeg,psvd,0.5)
    cd([resdir 'dtf\onesec\real\'])
    dtf(eeg,vdisc,1)
    cd([resdir 'dtf\onesec\control\'])
    dtf(eeg,psvd,1)
    
    waitbar(o/sf)
end
close(wb)



% -------------------------------------------------------------------------
% WAVELET
% -------------------------------------------------------------------------
function [pow,phase,f] = eeg_wavelet(dat)

% Prepare for wavelet transformation
variance = std(dat) ^ 2;
dat = (dat - mean(dat)) / sqrt(variance) ;
n = length(dat);
dt = 1 / 1000;
pad = 1;
dj = 0.08;    
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
pow = abs(wave).^2;
phase = angle(wave);



% -------------------------------------------------------------------------
function [pow,phase,f] = unit_wavelet(vdisc,lenu)

% Sinc convolution
fs = 10000;     % unit
dto = 1 / fs;
ts = zeros(1,lenu);
ts(vdisc) = 1;
du = diff(vdisc);
fdu = 1 ./ du;
fdu = [fdu 0.0001];
fcut = 100; 
fsnew = 1000;
dtnew = 1 / 1000;
fsold = 10000;
fsratio = fsnew / fsold;
told = vdisc * dto * fcut;
tnew = (1:lenu*fsratio) * dtnew * fcut;
lentold = length(told);
zint = 0;
for i = 1:lentold
    zint = zint + sinc(tnew-told(i));
end

% Prepare for wavelet transformation
variance = std(zint) ^ 2;
zint = (zint - mean(zint)) / sqrt(variance) ;
n = length(zint);
dt = 1 / 1000;
pad = 1;
dj = 0.08;    
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
[wave,period,scale,coi] = b_wavelet_new3(zint,dt,pad,dj,s0,j1,mother,param,mis);
pow = abs(wave).^2;
phase = angle(wave);



% -------------------------------------------------------------------------
% ENTROPY
% -------------------------------------------------------------------------
function entropy_line(data1,data2,f,WindowSize,Overlap)

[k1 k2] = size(data1);

winlen = WindowSize;   % window size
maxi = floor(k2/winlen);

aHx = []; aHy = []; aHxy = [];     % preallocating output variables
aHycx = []; aHxcy = [];
aIxy = []; aUxy = []; aUyx = [];
aIxynorm = [];
aRelHy = []; aRelHx = [];

% Entropy calculation
ovlp = Overlap;
for i = 1:maxi*ovlp-ovlp+1        % ABS LOOP
    inx1 = (i - 1) * winlen / ovlp + 1;  % Note: overlaping windows!
    inx1 = round(inx1);
    inx2 = inx1 + winlen - 1;

    y1 = data1(:,inx1:inx2);   % eeg wavelet "cells"
    y2 = data2(:,inx1:inx2);   % unit wavelet "cells"
    numy = numel(y1);
    bno = fix(exp(0.626+0.4*log(numy-1)));   % number of bins
    [hx,hy,jh,Hx,Hy,Hxy,Hxcy,Hycx,Ixy,Uxy,Uyx] = b_entry_sub1(data1,data2,y1,bno,y2,bno); % MAIN
    Ixynorm = Ixy / log(numy);
    RelHy = Hy / log(bno);      % relative Shannon entropy
    RelHx = Hx / log(bno);

    aHx = [aHx Hx];
    aHy = [aHy Hy];
    aHxy = [aHxy Hxy];
    aHycx = [aHycx Hycx];
    aHxcy = [aHxcy Hxcy];
    aIxy = [aIxy Ixy];
    aUxy = [aUxy Uxy];
    aUyx = [aUyx Uyx];
    aIxynorm = [aIxynorm Ixynorm];
    aRelHy = [aRelHy RelHy];
    aRelHx = [aRelHx RelHx];
end     % end of abs loop

% Save
global FILENAME
fnm = FILENAME;
global DATINX1
datinx1 = DATINX1;
global DATINX2
datinx2 = DATINX2;

fnts = [fnm '_ENTROPYline'];
save(fnts,'aHx','aHy','aHxy','aIxy','aIxynorm','aRelHy','aRelHx',...
    'aHxcy','aHycx','aUxy','aUyx');

% Plot
sr = 10000;
dsr = 1000;
const = sr / dsr;
x1 = datinx1 / sr;
x2 = (datinx1 + maxi * winlen * const) / sr;
xa = linspace(x1,x2,length(aHx));
x11 = datinx1 / sr;
x22 = datinx2 / sr;
xb = linspace(x11,x22,datinx2-datinx1);

figure
subplot(3,1,1)
plot(xa,aUyx)
y_lim = ylim;
axis([x1 x2 y_lim(1) y_lim(2)]);
title('unit->eeg')

subplot(3,1,2)
plot(xa,aUxy)
y_lim = ylim;
axis([x1 x2 y_lim(1) y_lim(2)]);
title('eeg->unit')

subplot(3,1,3)
plot(xa,aUyx-aUxy)
y_lim = ylim;
axis([x1 x2 y_lim(1) y_lim(2)]);
title('difference')

% Save
saveas(gcf,[fnm '_ENTROPY1line'],'fig');

% Plot
figure
subplot(9,1,1)
plot(xa,aHx)
y_lim = ylim;
axis([x1 x2 y_lim(1) y_lim(2)]);
title('H(unit)')

subplot(9,1,2)
plot(xa,aHy)
y_lim = ylim;
axis([x1 x2 y_lim(1) y_lim(2)]);
title('H(eeg)')

subplot(9,1,3)
plot(xa,aHxy)
y_lim = ylim;
axis([x1 x2 y_lim(1) y_lim(2)]);
title('H(unit,eeg)')

subplot(9,1,4)
plot(xa,aIxy)
y_lim = ylim;
axis([x1 x2 y_lim(1) y_lim(2)]);
title('I(unit&eeg)')

subplot(9,1,5)
plot(xa,aIxynorm)
y_lim = ylim;
axis([x1 x2 y_lim(1) y_lim(2)]);
title('I(unit&eeg) normalized')

subplot(9,1,6)
plot(xa,aRelHx)
y_lim = ylim;
axis([x1 x2 y_lim(1) y_lim(2)]);
title('H(unit) relative')

subplot(9,1,7)
plot(xa,aRelHy)
y_lim = ylim;
axis([x1 x2 y_lim(1) y_lim(2)]);
title('H(eeg) relative')

subplot(9,1,8)
plot(xa,aHxcy)
y_lim = ylim;
axis([x1 x2 y_lim(1) y_lim(2)]);
title('H(unit|eeg)')

subplot(9,1,9)
plot(xa,aHycx)
y_lim = ylim;
axis([x1 x2 y_lim(1) y_lim(2)]);
title('H(eeg|unit)')

% Save
saveas(gcf,[fnm '_ENTROPY2line'],'fig');
close all



% -------------------------------------------------------------------------
function entropy_image(data1,data2,pwind1,pwind2,f,WindowSize,Overlap,fnm)

[k1 k2] = size(data1);

winlen = WindowSize;   % window size
maxi = floor(k2/winlen);

aHx = []; aHy = []; aHxy = [];     % preallocating output variables
aHycx = []; aHxcy = [];
aIxy = []; aUxy = []; aUyx = [];
aIxynorm = [];
aRelHy = []; aRelHx = [];

% Entropy calculation
ovlp = Overlap;
for j = 1:k1
    tHx = []; tHy = []; tHxy = [];  % preallocating temporal variables
    tHycx = []; tHxcy = [];
    tIxy = []; tUxy = []; tUyx = [];
    tIxynorm = [];
    tRelHy = []; tRelHx = [];
    for i = 1:maxi*ovlp-ovlp+1        % ABS LOOP
        inx1 = (i - 1) * winlen / ovlp + 1;  % Note: overlaping windows!
        inx1 = round(inx1);
        inx2 = inx1 + winlen - 1;

        y1 = data1(j,inx1:inx2);   % eeg wavelet "cells"
        y2 = data2(j,inx1:inx2);   % unit wavelet "cells"
        numy = numel(y1);
        bno = fix(exp(0.626+0.4*log(numy-1)));   % number of bins
        [hx,hy,jh,Hx,Hy,Hxy,Hxcy,Hycx,Ixy,Uxy,Uyx] = b_entry_sub1(data1,data2,y1,bno,y2,bno); % MAIN
        Ixynorm = Ixy / log(numy);
        RelHy = Hy / log(bno);      % relative Shannon entropy
        RelHx = Hx / log(bno);

        tHx = [tHx Hx];
        tHy = [tHy Hy];
        tHxy = [tHxy Hxy];
        tHycx = [tHycx Hycx];
        tHxcy = [tHxcy Hxcy];
        tIxy = [tIxy Ixy];
        tUxy = [tUxy Uxy];
        tUyx = [tUyx Uyx];
        tIxynorm = [tIxynorm Ixynorm];
        tRelHy = [tRelHy RelHy];
        tRelHx = [tRelHx RelHx];
    end
    aHx(j,:) = tHx;
    aHy(j,:) = tHy;
    aHxy(j,:) = tHxy;
    aHycx(j,:) = tHycx;
    aHxcy(j,:) = tHxcy;
    aIxy(j,:) = tIxy;
    aUxy(j,:) = tUxy;
    aUyx(j,:) = tUyx;
    aIxynorm(j,:) = tIxynorm;
    aRelHy(j,:) = tRelHy;
    aRelHx(j,:) = tRelHx;
end     % end of abs loop

% Plot
sr = 10000;
dsr = 1000;
const = sr / dsr;
global FILENAME
fnm = FILENAME;
global DATINX1
datinx1 = DATINX1;
global DATINX2
datinx2 = DATINX2;

x1 = datinx1 / sr;
x2 = (datinx1 + maxi * winlen * const) / sr;
xa = linspace(x1,x2,length(aHx));
x11 = datinx1 / sr;
x22 = datinx2 / sr;
xb = linspace(x11,x22,datinx2-datinx1);

figure
subplot(3,1,1)
imagesc([x1 x2],[pwind1 pwind2],aUyx,[0 1])
b_rescaleaxis('Y',f)
title('unit->eeg')

subplot(3,1,2)
imagesc([x1 x2],[pwind1 pwind2],aUxy,[0 1])
b_rescaleaxis('Y',f)
title('eeg->unit')

subplot(3,1,3)
imagesc([x1 x2],[pwind1 pwind2],aUyx-aUxy,[-1 1])
b_rescaleaxis('Y',f)
title('difference')

% Save
fnts = [fnm '_ENTROPY1image'];
saveas(gcf,fnts,'jpg');

figure
subplot(9,1,1)
imagesc([x1 x2],[pwind1 pwind2],aHx,[0 1])
b_rescaleaxis('Y',f)
title('H(unit)')

subplot(9,1,2)
imagesc([x1 x2],[pwind1 pwind2],aHy,[0 1])
b_rescaleaxis('Y',f)
title('H(eeg)')

subplot(9,1,3)
imagesc([x1 x2],[pwind1 pwind2],aHxy,[0 1])
b_rescaleaxis('Y',f)
title('H(unit,eeg)')

subplot(9,1,4)
imagesc([x1 x2],[pwind1 pwind2],aIxy,[0 1])
b_rescaleaxis('Y',f)
title('I(unit&eeg)')

subplot(9,1,5)
imagesc([x1 x2],[pwind1 pwind2],aIxynorm,[0 1])
b_rescaleaxis('Y',f)
title('I(unit&eeg) normalized')

subplot(9,1,6)
imagesc([x1 x2],[pwind1 pwind2],aRelHx,[0 1])
b_rescaleaxis('Y',f)
title('H(unit) relative')

subplot(9,1,7)
imagesc([x1 x2],[pwind1 pwind2],aRelHy,[0 1])
b_rescaleaxis('Y',f)
title('H(eeg) relative')

subplot(9,1,8)
imagesc([x1 x2],[pwind1 pwind2],aHxcy,[0 1])
b_rescaleaxis('Y',f)
title('H(unit|eeg)')

subplot(9,1,9)
imagesc([x1 x2],[pwind1 pwind2],aHycx,[0 1])
b_rescaleaxis('Y',f)
title('H(eeg|unit)')

% Save
fnts = [fnm '_ENTROPY2image'];
saveas(gcf,fnts,'jpg');
close all



% -------------------------------------------------------------------------
% DTF
% -------------------------------------------------------------------------
function dtf(eeg,vdisc,c)

% Set number of frequencies
global N
N = 256;
dsr = 500;

% Sampling rate of raw data
sr = 10000;

% Initialize global LIN
global LIN1
global LIN2
LIN1 = [];
LIN2 = [];

% Get the analysed segment
global FILENAME
fnm = FILENAME;
global DATINX1
datinx1 = DATINX1;
global DATINX2
datinx2 = DATINX2;
eeg = eeg(datinx1:datinx2);
vdisc = vdisc(find(vdisc>=datinx1&vdisc<=datinx2)) - datinx1;

% Main
seglen = sr * c;
maxi = floor(length(eeg)/seglen);
PDC = zeros(2,2,N,maxi);
for i = 1:maxi
    ind1 = (i - 1) * seglen + 1;
    ind2 = ind1 + seglen -1;
    eeg2 = eeg(ind1:ind2);
    vdisc2 = vdisc(find(vdisc>=ind1&vdisc<=ind2)) - ind1;
    
    if ~isempty(vdisc2)
        [PDC(:,:,:,i),f] = dircoh(eeg2,vdisc2,sr,dsr);  % partial directed coherence
    else
        PDC(:,:,:,i) = 0;
        LIN1(end+1) = NaN;
        LIN2(end+1) = NaN;
    end
end

% Plot LIN
global LIN1
global LIN2
figure
time = linspace(datinx1/10000,datinx2/10000,length(LIN1));
plot(time,LIN1,'r')  % unit->eeg
hold on
plot(time,LIN2,'b')  % eeg->unit

% Save
saveas(gcf,[fnm '_DTFline'],'fig');

% Create image (time-frequency)
psu = zeros(1,length(eeg));     % pseudounit
psu(vdisc) = 1;
psu(vdisc+1) = -1;

PDC_ue = squeeze(PDC(1,2,:,:));    % unit->eeg
PDC_eu = squeeze(PDC(2,1,:,:));    % eeg->unit
clow = 0;       % common scaling for the images!
clow2 = -1;
chigh = 1;
x11 = datinx1 / sr;
x22 = floor(datinx2/seglen) * seglen / sr;
x1 = (datinx1 + (seglen / 2)) / sr;
x2 = (floor(datinx2/seglen) * seglen - (seglen / 2)) / sr;
x = [x1 x2];
ls = 10;   % last plotted scale
y = [f(1) f(ls)];
xb = linspace(x11,datinx2/sr,length(eeg));
figure
subplot(4,1,1)
imagesc(x,y,PDC_ue(1:ls,:),[clow chigh])
title('unit->eeg')
subplot(4,1,2)
imagesc(x,y,PDC_eu(1:ls,:),[clow chigh])
title('eeg->unit')
subplot(4,1,3)
imagesc(x,y,PDC_ue(1:ls,:)-PDC_eu(1:ls,:),[clow2 chigh])
title('difference')
subplot(4,1,4)
plot(xb,eeg)
hold on
plot(xb,psu,'r')
y_lim = ylim;
axis([x11,x22,y_lim(1),y_lim(2)]);
title('unit:red')

% Save
fnts = [fnm '_DTFimage'];
saveas(gcf,fnts,'jpg');
close all

fnts = [fnm '_DTF'];
save(fnts,'PDC_eu','PDC_ue');



% -------------------------------------------------------------------------
function [PDC,f] = dircoh(eeg,vdisc,sr,dsr)

% Get number of frequencies
global N

% Instantanous frequency
if vdisc(1) == 0
    vdisc = vdisc(2:end);
end
zint = ifreq(vdisc,length(eeg));
zint = zint(1:sr/dsr:end);      % downsampling instantenous frequency on 'dsr'

% Estimate optimal model order
en = eeg(1:sr/dsr:end)';        % downsampling EEG on 'dsr'
un = zint';
if length(en) > length(un)
    en = en(1:end-1);
elseif length(un) > length(en)
    un = un(1:end-1);
end
% en = (en - mean(en)) / std(en);      % normalization
% un = (un - mean(un)) / std(un);
en = en ./ max(en);      % normalization
un = un ./ max(un);
s = [en un];
p = arorder(s,1,32);

% Calculate DTF and PDC
global N;    % number of frequencies
Fs = dsr;  % sampling rate
[AR,RC,PE] = mvar(s,p,5);
X.A = [eye(size(AR,1)),-AR]; 
X.B = eye(size(X.A,1));
X.C = eye(size(X.A,1));
X.datatype = 'MVAR';
[S,h,PDC,COH,DTF,f] = main(X,'DTF',N,Fs);



% -------------------------------------------------------------------------
function instfrek = ifreq(vdisc,lenu)
instfrek = [];
isi = diff(vdisc);
for i = 1:length(vdisc)-1
    instfrek(vdisc(i):vdisc(i+1)) = 1 / isi(i);
end
instfrek(1:vdisc(1)-1) = 1 / vdisc(1);
instfrek(vdisc(end):lenu) = 1 / (lenu - vdisc(end));



% -------------------------------------------------------------------------
function popt = arorder(v,pmin,pmax)
% See ARFIT for a detailed help.

% Input argument check
error(nargchk(3,3,nargin))
if (pmin ~= round(pmin) | pmax ~= round(pmax))
    error('Order must be integer.');
end
if (pmax < pmin)
    error('PMAX must be greater than or equal to PMIN.')
end

% n: number of observations; m: dimension of state vectors
[n,m] = size(v);

mcor = 1;       % fit intercept vector
selector = 'sbc';       % use SBC as order selection criterion

ne = n - pmax;      % number of block equations of size m
npmax = m * pmax + mcor;        % maximum number of parameter vectors of length m

if (ne <= npmax)
    error('Time series too short.')
end

% Compute QR factorization for model of order pmax
[R, scale] = arqr(v, pmax, mcor);

% Compute approximate order selection criteria for models
% of order pmin:pmax
[sbc, fpe] = arord(R, m, mcor, ne, pmin, pmax);

% Get index iopt of order that minimizes the order selection
% criterion specified by the variable selector
[val, iopt] = min(eval(selector));

% Select order of model
popt = pmin + iopt - 1;     % estimated optimum order
np = m * popt + mcor;       % number of parameter vectors of length m



% -------------------------------------------------------------------------
function [S,h,PDC,COH,DTF,f] = main(X,Mode,N,Fs)
% See PLOTA for a detailed help.

% Initialize
[K1,K2] = size(X.A);
p = K2 / K1 - 1;
[K1,K2] = size(X.B);
q = K2 / K1 - 1;
f = (1:N) / N / 2 * Fs;
z = i * 2 * pi / Fs;

h = zeros(K1,K1,N);
S = zeros(K1,K1,N);
DTF = zeros(K1,K1,N);
COH = zeros(K1,K1,N);
PDC = zeros(K1,K1,N);
PDCF = zeros(K1,K1,N);
invC = inv(X.C);
tmp1 = zeros(1,K1);
tmp2 = zeros(1,K1);

% Ask global LIN
global LIN1
global LIN2

% Calculate PDC
for n = 1:N
    atmp = zeros(K1,K1);
    for k = 1:p+1,
        atmp = atmp + X.A(:,k*K1+(1-K1:0))*exp(z*(k-1)*f(n));
    end
    btmp = zeros(K1,K2);
    for k = 1:q+1,
        btmp = btmp + X.B(:,k*K1+(1-K1:0))*exp(z*(k-1)*f(n));
    end
    h(:,:,n) = atmp \ btmp;
    S(:,:,n) = h(:,:,n)*X.C*h(:,:,n)';

    for k1 = 1:K1
        tmp = squeeze(atmp(:,k1));
        tmp1(k1) = sqrt(tmp'*tmp);
        tmp2(k1) = sqrt(tmp'*invC*tmp);
    end

    PDCF(:,:,n) = abs(atmp) ./ tmp2(ones(1,K1),:);
    PDC(:,:,n)  = abs(atmp) ./ tmp1(ones(1,K1),:);
end

% Calculate DTF
for k1 = 1:K1
    DEN = sqrt(sum(abs(h(k1,:,:)).^2,2));
    for k2 = 1:K2
        COH(k1,k2,:) = abs(S(k1,k2,:)) ./ sqrt(abs(S(k1,k1,:).*S(k2,k2,:)));
        DTF(k1,k2,:) = abs(h(k1,k2,:)) ./ DEN;
    end
end

% Plot
for k1 = 1:K1
    for k2 = 1:K2
        sqPDC = squeeze(PDC(k1,k2,:));
        if k1 == 1 & k2 == 2
            LIN1(end+1) = mean(sqPDC(3:6));  % unit -> EEG
        end
        if k1 == 2 & k2 == 1
            LIN2(end+1) = mean(sqPDC(3:6));  % EEG -> unit
        end
    end
end



% -------------------------------------------------------------------------
% SUBDIRECTORIES
% -------------------------------------------------------------------------
function create_subdir

if ~b_isdir2('entropy')
    mkdir entropy
end
cd entropy
if ~b_isdir2('power')
    mkdir power
end
if ~b_isdir2('phase')
    mkdir phase
end
cd power
if ~b_isdir2('line')
    mkdir line
end
cd line
if ~b_isdir2('windowsize1000_overlap1')
    mkdir windowsize1000_overlap1
end
cd windowsize1000_overlap1
if ~b_isdir2('real')
    mkdir real
end
if ~b_isdir2('control')
    mkdir control
end
cd ..
if ~b_isdir2('windowsize1000_overlap2')
    mkdir windowsize1000_overlap2
end
cd windowsize1000_overlap2
if ~b_isdir2('real')
    mkdir real
end
if ~b_isdir2('control')
    mkdir control
end
cd ..
cd ..
if ~b_isdir2('image')
    mkdir image
end
cd image
if ~b_isdir2('windowsize1000_overlap1')
    mkdir windowsize1000_overlap1
end
cd windowsize1000_overlap1
if ~b_isdir2('real')
    mkdir real
end
if ~b_isdir2('control')
    mkdir control
end
cd ..
if ~b_isdir2('windowsize1000_overlap2')
    mkdir windowsize1000_overlap2
end
cd windowsize1000_overlap2
if ~b_isdir2('real')
    mkdir real
end
if ~b_isdir2('control')
    mkdir control
end
cd ..
cd ..
cd ..
cd phase
if ~b_isdir2('line')
    mkdir line
end
cd line
if ~b_isdir2('windowsize1000_overlap1')
    mkdir windowsize1000_overlap1
end
cd windowsize1000_overlap1
if ~b_isdir2('real')
    mkdir real
end
if ~b_isdir2('control')
    mkdir control
end
cd ..
if ~b_isdir2('windowsize1000_overlap2')
    mkdir windowsize1000_overlap2
end
cd windowsize1000_overlap2
if ~b_isdir2('real')
    mkdir real
end
if ~b_isdir2('control')
    mkdir control
end
cd ..
cd ..
if ~b_isdir2('image')
    mkdir image
end
cd image
if ~b_isdir2('windowsize1000_overlap1')
    mkdir windowsize1000_overlap1
end
cd windowsize1000_overlap1
if ~b_isdir2('real')
    mkdir real
end
if ~b_isdir2('control')
    mkdir control
end
cd ..
if ~b_isdir2('windowsize1000_overlap2')
    mkdir windowsize1000_overlap2
end
cd windowsize1000_overlap2
if ~b_isdir2('real')
    mkdir real
end
if ~b_isdir2('control')
    mkdir control
end
cd ..
cd ..
cd ..
cd ..

if ~b_isdir2('dtf')
    mkdir dtf
end
cd dtf
if ~b_isdir2('halfsec')
    mkdir halfsec
end
cd halfsec
if ~b_isdir2('real')
    mkdir real
end
if ~b_isdir2('control')
    mkdir control
end
cd ..
if ~b_isdir2('onesec')
    mkdir onesec
end
cd onesec
if ~b_isdir2('real')
    mkdir real
end
if ~b_isdir2('control')
    mkdir control
end
cd ..
cd ..