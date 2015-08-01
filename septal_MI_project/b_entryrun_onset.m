function b_entryrun_onset
%ENTRYRUN_ONSET   Runs entropy on a sequence of theta segments.
%   ENTRYRUN_ONSET calculates wavelet power and wavelet phase entropy for
%   theta segments using 1000-point non-overlapping time windows. It saves
%   line and time-frequency image visualization and entropy matrices as well.
%
%   ENTRYRUN_ONSET calculates entropy for the first sec. of theta segment,
%   without downsampling (resolution: 100 ms).
%
%   See also ENTRY, ENTRYRUN3 and ENTRYRUN_THETA2.

% Input argument check
error(nargchk(0,0,nargin))

% Directories
global DATAPATH
global DATADIR
global DATADIR2
inpdir = [DATAPATH 'Burst\Cluster\Theta\'];   % input directory
resdir = [DATAPATH 'Entry_onset\Theta\'];
datadir = DATADIR2;

mm = pwd;
cd(resdir)
create_subdir

% Import
files = b_filelist(inpdir);
sf = length(files);
[datalist, dlist_short] = b_filelist(datadir);

% Progress indicator
wb = waitbar(0,'Running ENTRYRUN ONSET...','Position',[360 250 275 50]);
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
    
    cmps = strread(fname,'%s','delimiter','_MH');   % M, H: delimiters for 3ch data
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
    [wavea_abs,wavea_phase,f] = eeg_wavelet(eeg(i_first:min(i_second,i_first+10000))); % eeg (without downsampling)

% Sinc convolution & wavelet transformation of unit
    [waveb_abs,waveb_phase,f] = unit_wavelet(vdisc(find(vdisc>=i_first&vdisc<=min(i_second,i_first+10000)))-i_first,10001);   % unit
    [wavec_abs,wavec_phase,f] = unit_wavelet(psvd(find(psvd>=i_first&psvd<=min(i_second,i_first+10000)))-i_first,10001 );    % random unit
    
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
    entropy_line(abs1,abs2,f,1000,1,'power')
    cd([resdir 'entropy\phase\line\windowsize1000_overlap1\real\'])
    entropy_line(phase1,phase2,f,1000,1,'phase')
    cd([resdir 'entropy\power\line\windowsize1000_overlap1\control\'])
    entropy_line(abs1,abs3,f,1000,1,'power')
    cd([resdir 'entropy\phase\line\windowsize1000_overlap1\control\'])
    entropy_line(phase1,phase3,f,1000,1,'phase')
        
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
    cd(resdir)
    
    clear wavea_abs wavea_phase waveb_abs waveb_phase wavec_abs wavec_phase r s
    
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
dt = 1 / 10000;
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
fsnew = 10000;
dtnew = 1 / fsnew;
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
dt = 1 / 10000;
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
function entropy_line(data1,data2,f,WindowSize,Overlap,porp)

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
    if isequal(porp,'power')
        [hx,hy,jh,Hx,Hy,Hxy,Hxcy,Hycx,Ixy,Uxy,Uyx] = b_entry_sub1(data1,data2,y1,bno,y2,bno); % MAIN
    elseif isequal(porp,'phase')
        [hx,hy,jh,Hx,Hy,Hxy,Hxcy,Hycx,Ixy,Uxy,Uyx] = b_entry_sub1(data1,data2,y1,bno,y2,bno); % MAIN
    else
        error('technical error 287')
    end
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
if ~b_isdir2('windowsize5000_overlap5')
    mkdir windowsize5000_overlap5
end
cd windowsize5000_overlap5
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
if ~b_isdir2('windowsize5000_overlap5')
    mkdir windowsize5000_overlap5
end
cd windowsize5000_overlap5
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
if ~b_isdir2('windowsize5000_overlap5')
    mkdir windowsize5000_overlap5
end
cd windowsize5000_overlap5
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
if ~b_isdir2('windowsize5000_overlap5')
    mkdir windowsize5000_overlap5
end
cd windowsize5000_overlap5
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