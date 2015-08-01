function b_entryrun3c_inv
%ENTRYRUN3C_INV   Runs entropy and PDC on a sequence of recordings.
%   ENTRYRUN3C_INV calculates wavelet power and wavelet phase entropy for
%   theta segments using 1000-point non-overlapping and 50% overlapping
%   time windows. It saves line and time-frequency image visualization and
%   entropy matrices as well.
%
%   It also calculates Partial Directed Coherence and saves PDC images, 
%   lines and PDC matrices.
%
%   ENTRYRUN3C_INV attaches non-theta and theta segment together and
%   completes the calculation on the attached segments. (Note that the
%   difference between ENTRYRUN3C and ENTRYRUN3C_INV is the order of theta
%   and non-theta segment.)
%
%   It computes frequency adjusted random Poisson controls for entropy
%   and PDC calculations. Frequency is adjusted for theta and non-theta
%   segment separately, in contrast to ENTRYRUN3.
%
%   Panzeri-Treves bias correction is implemented for entropy and mutual
%   information calculation.
%
%   The program is implemented for 3-channel  data files as well.
%
%   Reference: Panzeri S, Senatore R, Montemurro MA, Petersen RS (2007)
%   Correcting for the sampling bias problem in spike train information
%   measures. Journal of neurophysiology 98:1064-1072
%
%   See also DTF, 2DTF, ENTRY, ENTRYRUN3, ENTRYRUN_THETA and ENTRYRUN_NOTH.

% Input argument check
error(nargchk(0,0,nargin))

% Directories
global DATAPATH
global DATADIR
global DATADIR2
inpdir1 = [DATAPATH 'Burst\Cluster\Theta_newcells\'];   % input directory
inpdir2 = [DATAPATH 'Burst\Cluster\Noth_newcells\'];   % input directory
resdir = [DATAPATH 'Entry3c_inv\Theta_newcells\'];
datadir = [DATADIR 'mi_zshift_data\discriminated\'];

mm = pwd;
cd(resdir)
create_subdir

% Import
[files1 files1_short] = b_filelist(inpdir1);
[files2 files2_short] = b_filelist(inpdir2);
files_short = intersect(files1_short,files2_short);

sf = length(files_short);
[datalist, dlist_short] = b_filelist(datadir);

% Progress indicator
wb = waitbar(0,'Running ENTRYRUN3C INV...','Position',[360 250 275 50]);
global WB
WB(end+1) = wb;

% Load
for o = 1:sf
    ftheta = files1(find(strcmp(files_short{o},files1_short))).name;
    fnoth = files2(find(strcmp(files_short{o},files2_short))).name;
    filenam = ftheta(1:6);
    if ~isequal(filenam,fnoth(1:6))
        error('Technical error 47')
    end
    cmps = strread(ftheta,'%s','delimiter','_MH');   % M, H: delimiters for 3ch data
    seglen_theta = str2num(cmps{5}) - str2num(cmps{4});
    i_first_theta = str2num(cmps{4});
    i_second_theta = str2num(cmps{5});
    cmps = strread(ftheta,'%s','delimiter','_');
    ftheta = cmps{1};
    for kt = 2:length(cmps)-1
        ftheta = [ftheta '_' cmps{kt}];
    end
    cmps = strread(fnoth,'%s','delimiter','_MH');   % M, H: delimiters for 3ch data
    seglen_noth = str2num(cmps{5}) - str2num(cmps{4});
    i_first_noth = str2num(cmps{4});
    i_second_noth = str2num(cmps{5});
    cmps = strread(fnoth,'%s','delimiter','_');
    fnoth = cmps{1};
    for kn = 2:length(cmps)-1
        fnoth = [fnoth '_' cmps{kn}];
    end
    
    global FTHETA
    FTHETA = ftheta;
    global FNOTH
    FNOTH = fnoth;
    global DATINX1_THETA
    DATINX1_THETA = i_first_theta;
    global DATINX2_THETA
    DATINX2_THETA = i_second_theta;
    global DATINX1_NOTH
    DATINX1_NOTH = i_first_noth;
    global DATINX2_NOTH
    DATINX2_NOTH = i_second_noth;
    
    inx = find(strcmp(filenam,dlist_short));
    fn = datalist(inx).name;
    ffn = [datadir fn];
    load(ffn);      % load raw data
    
% Create new data (attach theta and non-theta segment)
    eeg_theta = eeg(i_first_theta:i_second_theta);
    eeg_noth = eeg(i_first_noth:i_second_noth);
    eeg_new = [eeg_noth eeg_theta];
    vdisc_theta = vdisc(find(vdisc>=i_first_theta&vdisc<=i_second_theta));
    vdisc_noth = vdisc(find(vdisc>=i_first_noth&vdisc<=i_second_noth));
    vdisc_new = [vdisc_noth-i_first_noth vdisc_theta-i_first_theta+i_second_noth-i_first_noth];
    seglen_new = seglen_theta + seglen_noth;
    
% Create random unit (frequency-adjusted Poisson process)
    leneeg = length(eeg_new);
    
    vdisc_theta = vdisc(find(vdisc>=i_first_theta&vdisc<=i_second_theta));
    vdisc_noth = vdisc(find(vdisc>=i_first_noth&vdisc<=i_second_noth));
    frq_theta = length(vdisc_theta) / (i_second_theta - i_first_theta);
    frq_noth = length(vdisc_noth) / (i_second_noth - i_first_noth);
    
    r = random('exp',1/frq_theta,1,10000);     % 1/lambda param. exp.!!!! (MATLAB-ban forditva van...)
    s = cumsum(r);
    psvd1 = unique(ceil(s));     % 'pseudo vdisc'
    psvd_theta = seglen_noth + psvd1(find(psvd1<seglen_theta));   % expected value of length(psvd): leneeg*lambda
    if isequal(length(psvd1),length(psvd_theta))
        error('Technical error 118')
    end
    
    r = random('exp',1/frq_noth,1,10000);
    s = cumsum(r);
    psvd2 = unique(ceil(s));     % 'pseudo vdisc'
    psvd_noth = psvd2(find(psvd2<seglen_noth));   % expected value of length(psvd): leneeg*lambda
    if isequal(length(psvd2),length(psvd_noth))
        error('Technical error 126')
    end
    psvd_new = [psvd_noth psvd_theta];
    
% Downsample & wavelet transformation of eeg
    [wavea_abs,wavea_phase,f] = eeg_wavelet(eeg_new(1:10:end)); % eeg (downsample on 1000 Hz)

% Sinc convolution & wavelet transformation of unit
    [waveb_abs,waveb_phase,f] = unit_wavelet(vdisc_new,seglen_new);   % unit
    [wavec_abs,wavec_phase,f] = unit_wavelet(psvd_new,seglen_new);    % random unit
    
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
%     cd([resdir 'entropy\phase\line\windowsize1000_overlap1\real\'])
%     entropy_line(phase1,phase2,f,1000,1,'phase')
    cd([resdir 'entropy\power\line\windowsize1000_overlap1\control\'])
    entropy_line(abs1,abs3,f,1000,1,'power')
%     cd([resdir 'entropy\phase\line\windowsize1000_overlap1\control\'])
%     entropy_line(phase1,phase3,f,1000,1,'phase')
    
%     fnd = find(f>6);    % frequency band bounderies
%     pwind1 = fnd(end);
%     fnd = find(f<2.5);
%     pwind2 = fnd(1);
%     
%     abs1 = wavea_abs(pwind1:pwind2,:);  % eeg
%     phase1 = wavea_phase(pwind1:pwind2,:);
%     abs2 = waveb_abs(pwind1:pwind2,:);  % unit
%     phase2 = waveb_phase(pwind1:pwind2,:);
%     abs3 = wavec_abs(pwind1:pwind2,:);  % random unit
%     phase3 = wavec_phase(pwind1:pwind2,:);
%     
%     cd([resdir 'entropy\power\image\windowsize1000_overlap1\real\'])
%     entropy_image(abs1,abs2,pwind1,pwind2,f,1000,1)
%     cd([resdir 'entropy\phase\image\windowsize1000_overlap1\real\'])
%     entropy_image(phase1,phase2,pwind1,pwind2,f,1000,1)
%     cd([resdir 'entropy\power\image\windowsize1000_overlap1\control\'])
%     entropy_image(abs1,abs3,pwind1,pwind2,f,1000,1)
%     cd([resdir 'entropy\phase\image\windowsize1000_overlap1\control\'])
%     entropy_image(phase1,phase3,pwind1,pwind2,f,1000,1)
%     cd(resdir)
    
    clear wavea_abs wavea_phase waveb_abs waveb_phase wavec_abs wavec_phase r s
    
% Directed Transfer Function
%     cd([resdir 'dtf\halfsec\real\'])
%     dtf(eeg_new,vdisc_new,0.5)
%     cd([resdir 'dtf\halfsec\control\'])
%     dtf(eeg_new,psvd_new,0.5)
%     cd([resdir 'dtf\onesec\real\'])
%     dtf(eeg_new,vdisc_new,1)
%     cd([resdir 'dtf\onesec\control\'])
%     dtf(eeg_new,psvd_new,1)
    
    waitbar(o/sf)
end
close(wb)
clear global FTHETA FNOTH DATINX1_THETA DATINX2_THETA DATINX1_NOTH DATINX2_NOTH



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
function entropy_line(data1,data2,f,WindowSize,Overlap,porp)

[k1 k2] = size(data1);

winlen = WindowSize;   % window size
maxi = floor(k2/winlen);

aHx = []; aHy = []; aHxy = [];     % preallocating output variables
aHycx = []; aHxcy = [];
aIxy = []; aUxy = []; aUyx = [];
aIxynorm = [];
aRelHy = []; aRelHx = [];
aBias_HR = []; aBias_I = [];

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
        [hx,hy,jh,Hx,Hy,Hxy,Hxcy,Hycx,Ixy,Uxy,Uyx,Bias_HR,Bias_I] =...
            entry_sub1(data1,data2,y1,bno,y2,bno); % MAIN
    elseif isequal(porp,'phase')
        [hx,hy,jh,Hx,Hy,Hxy,Hxcy,Hycx,Ixy,Uxy,Uyx,Bias_HR,Bias_I] =...
            entry_sub1(data1,data2,y1,bno,y2,bno); % MAIN
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
    aBias_HR = [aBias_HR Bias_HR];
    aBias_I = [aBias_I Bias_I];
end     % end of abs loop

% Save
global FTHETA
ftheta = FTHETA;
global FNOTH
fnoth = FNOTH;
global DATINX1_THETA
datinx1_theta = DATINX1_THETA;
global DATINX2_THETA
datinx2_theta = DATINX2_THETA;
global DATINX1_NOTH
datinx1_noth = DATINX1_NOTH;
global DATINX2_NOTH
datinx2_noth = DATINX2_NOTH;

sr = 10000;
dsr = 1000;
const = sr / dsr;

% Plot
x1 = 0;
x2 = (maxi * winlen * const) / sr;
xa = linspace(x1,x2,length(aHx));
lna = floor((datinx2_noth-datinx1_noth)/winlen) * winlen / sr;
ld = (datinx2_theta - datinx1_theta) + (datinx2_noth - datinx1_noth);
x11 = 0;
x22 = ld / sr;
xb = linspace(x11,x22,ld);
lnb = (datinx2_noth - datinx1_noth) / sr;

figure
subplot(3,1,1)
plot(xa,aUyx)
y_lim = ylim;
axis([x1 x2 y_lim(1) y_lim(2)]);
line([lna lna],[y_lim(1) y_lim(2)],'Color','green');
title('unit->eeg')

subplot(3,1,2)
plot(xa,aUxy)
y_lim = ylim;
axis([x1 x2 y_lim(1) y_lim(2)]);
line([lna lna],[y_lim(1) y_lim(2)],'Color','green');
title('eeg->unit')

subplot(3,1,3)
plot(xa,aUyx-aUxy)
y_lim = ylim;
axis([x1 x2 y_lim(1) y_lim(2)]);
line([lna lna],[y_lim(1) y_lim(2)],'Color','green');
title('difference')

% Save
fnm = ftheta(1:6);
saveas(gcf,[fnm '_ENTROPY1line'],'fig');

% Plot
figure
subplot(9,1,1)
plot(xa,aHx)
y_lim = ylim;
axis([x1 x2 y_lim(1) y_lim(2)]);
line([lna lna],[y_lim(1) y_lim(2)],'Color','green');
title('H(unit)')

subplot(9,1,2)
plot(xa,aHy)
y_lim = ylim;
axis([x1 x2 y_lim(1) y_lim(2)]);
line([lna lna],[y_lim(1) y_lim(2)],'Color','green');
title('H(eeg)')

subplot(9,1,3)
plot(xa,aHxy)
y_lim = ylim;
axis([x1 x2 y_lim(1) y_lim(2)]);
line([lna lna],[y_lim(1) y_lim(2)],'Color','green');
title('H(unit,eeg)')

subplot(9,1,4)
plot(xa,aIxy)
y_lim = ylim;
axis([x1 x2 y_lim(1) y_lim(2)]);
line([lna lna],[y_lim(1) y_lim(2)],'Color','green');
title('I(unit&eeg)')

subplot(9,1,5)
plot(xa,aIxynorm)
y_lim = ylim;
axis([x1 x2 y_lim(1) y_lim(2)]);
line([lna lna],[y_lim(1) y_lim(2)],'Color','green');
title('I(unit&eeg) normalized')

subplot(9,1,6)
plot(xa,aRelHx)
y_lim = ylim;
axis([x1 x2 y_lim(1) y_lim(2)]);
line([lna lna],[y_lim(1) y_lim(2)],'Color','green');
title('H(unit) relative')

subplot(9,1,7)
plot(xa,aRelHy)
y_lim = ylim;
axis([x1 x2 y_lim(1) y_lim(2)]);
line([lna lna],[y_lim(1) y_lim(2)],'Color','green');
title('H(eeg) relative')

subplot(9,1,8)
plot(xa,aHxcy)
y_lim = ylim;
axis([x1 x2 y_lim(1) y_lim(2)]);
line([lna lna],[y_lim(1) y_lim(2)],'Color','green');
title('H(unit|eeg)')

subplot(9,1,9)
plot(xa,aHycx)
y_lim = ylim;
axis([x1 x2 y_lim(1) y_lim(2)]);
line([lna lna],[y_lim(1) y_lim(2)],'Color','green');
title('H(eeg|unit)')

% Save
saveas(gcf,[fnm '_ENTROPY2line'],'fig');

figure
plot(xa,aUxy)
hold on
plot(xa,aUyx,'r')
y_lim = ylim;
axis([x1 x2 y_lim(1) y_lim(2)]);
line([lna lna],[y_lim(1) y_lim(2)],'Color','green');

% Save
saveas(gcf,[fnm '_ENTROPY3line'],'fig');
close all

% Save variables
lm = (floor((datinx2_theta-datinx1_theta)/winlen) * ovlp - ovlp + 1) / const;
lm = floor(lm);
aHx_theta = aHx(1:lm);
aHx_noth = aHx(lm+2:end);
aHy_theta = aHy(1:lm);
aHy_noth = aHy(lm+2:end);
aHxy_theta = aHxy(1:lm);
aHxy_noth = aHxy(lm+2:end);
aIxy_theta = aIxy(1:lm);
aIxy_noth = aIxy(lm+2:end);
aIxynorm_theta = aIxynorm(1:lm);
aIxynorm_noth = aIxynorm(lm+2:end);
aRelHy_theta = aRelHy(1:lm);
aRelHy_noth = aRelHy(lm+2:end);
aRelHx_theta = aRelHx(1:lm);
aRelHx_noth = aRelHx(lm+2:end);
aHxcy_theta = aHxcy(1:lm);
aHxcy_noth = aHxcy(lm+2:end);
aHycx_theta = aHycx(1:lm);
aHycx_noth = aHycx(lm+2:end);
aUxy_theta = aUxy(1:lm);
aUxy_noth = aUxy(lm+2:end);
aUyx_theta = aUyx(1:lm);
aUyx_noth = aUyx(lm+2:end);
aBias_HR_theta = aBias_HR(1:lm);
aBias_HR_noth = aBias_HR(lm+2:end);
aBias_I_theta = aBias_I(1:lm);
aBias_I_noth = aBias_I(lm+2:end);

aHx = aHx_theta;
aHy = aHy_theta;
aHxy = aHxy_theta;
aIxy = aIxy_theta;
aIxynorm = aIxynorm_theta;
aRelHy = aRelHy_theta;
aRelHx = aRelHx_theta;
aHxcy = aHxcy_theta;
aHycx = aHycx_theta;
aUxy = aUxy_theta;
aUyx = aUyx_theta;
aBias_HR = aBias_HR_theta;
aBias_I = aBias_I_theta;
fnts = [ftheta '_ENTROPYline'];
save(fnts,'aHx','aHy','aHxy','aIxy','aIxynorm','aRelHy','aRelHx',...
    'aHxcy','aHycx','aUxy','aUyx','aBias_HR','aBias_I');

aHx = aHx_noth;
aHy = aHy_noth;
aHxy = aHxy_noth;
aIxy = aIxy_noth;
aIxynorm = aIxynorm_noth;
aRelHy = aRelHy_noth;
aRelHx = aRelHx_noth;
aHxcy = aHxcy_noth;
aHycx = aHycx_noth;
aUxy = aUxy_noth;
aUyx = aUyx_noth;
aBias_HR = aBias_HR_noth;
aBias_I = aBias_I_noth;
fnts = [fnoth '_ENTROPYline'];
save(fnts,'aHx','aHy','aHxy','aIxy','aIxynorm','aRelHy','aRelHx',...
    'aHxcy','aHycx','aUxy','aUyx','aBias_HR','aBias_I');



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
global FTHETA
ftheta = FTHETA;
global FNOTH
fnoth = FNOTH;
global DATINX1_THETA
datinx1_theta = DATINX1_THETA;
global DATINX2_THETA
datinx2_theta = DATINX2_THETA;
global DATINX1_NOTH
datinx1_noth = DATINX1_NOTH;
global DATINX2_NOTH
datinx2_noth = DATINX2_NOTH;

x1 = 0;
x2 = (maxi * winlen * const) / sr;
xa = linspace(x1,x2,length(aHx));
lna = floor((datinx2_noth-datinx1_noth)/winlen) * winlen / sr;
ld = (datinx2_theta - datinx1_theta) + (datinx2_noth - datinx1_noth);
x11 = 0;
x22 = ld / sr;
xb = linspace(x11,x22,ld);
lnb = (datinx2_noth - datinx1_noth) / sr;

figure
subplot(3,1,1)
imagesc([x1 x2],[pwind1 pwind2],aUyx,[0 1])
b_rescaleaxis('Y',f)
y_lim = ylim;
line([lna lna],[y_lim(1) y_lim(2)],'Color','white');
title('unit->eeg')

subplot(3,1,2)
imagesc([x1 x2],[pwind1 pwind2],aUxy,[0 1])
b_rescaleaxis('Y',f)
y_lim = ylim;
line([lna lna],[y_lim(1) y_lim(2)],'Color','white');
title('eeg->unit')

subplot(3,1,3)
imagesc([x1 x2],[pwind1 pwind2],aUyx-aUxy,[-1 1])
b_rescaleaxis('Y',f)
y_lim = ylim;
line([lna lna],[y_lim(1) y_lim(2)],'Color','white');
title('difference')

% Save
fnm = ftheta(1:6);
fnts = [fnm '_ENTROPY1image'];
saveas(gcf,fnts,'jpg');

figure
subplot(9,1,1)
imagesc([x1 x2],[pwind1 pwind2],aHx)
b_rescaleaxis('Y',f)
y_lim = ylim;
line([lna lna],[y_lim(1) y_lim(2)],'Color','white');
title('H(unit)')

subplot(9,1,2)
imagesc([x1 x2],[pwind1 pwind2],aHy)
b_rescaleaxis('Y',f)
y_lim = ylim;
line([lna lna],[y_lim(1) y_lim(2)],'Color','white');
title('H(eeg)')

subplot(9,1,3)
imagesc([x1 x2],[pwind1 pwind2],aHxy)
b_rescaleaxis('Y',f)
y_lim = ylim;
line([lna lna],[y_lim(1) y_lim(2)],'Color','white');
title('H(unit,eeg)')

subplot(9,1,4)
imagesc([x1 x2],[pwind1 pwind2],aIxy)
b_rescaleaxis('Y',f)
y_lim = ylim;
line([lna lna],[y_lim(1) y_lim(2)],'Color','white');
title('I(unit&eeg)')

subplot(9,1,5)
imagesc([x1 x2],[pwind1 pwind2],aIxynorm)
b_rescaleaxis('Y',f)
y_lim = ylim;
line([lna lna],[y_lim(1) y_lim(2)],'Color','white');
title('I(unit&eeg) normalized')

subplot(9,1,6)
imagesc([x1 x2],[pwind1 pwind2],aRelHx)
b_rescaleaxis('Y',f)
y_lim = ylim;
line([lna lna],[y_lim(1) y_lim(2)],'Color','white');
title('H(unit) relative')

subplot(9,1,7)
imagesc([x1 x2],[pwind1 pwind2],aRelHy)
b_rescaleaxis('Y',f)
y_lim = ylim;
line([lna lna],[y_lim(1) y_lim(2)],'Color','white');
title('H(eeg) relative')

subplot(9,1,8)
imagesc([x1 x2],[pwind1 pwind2],aHxcy)
b_rescaleaxis('Y',f)
y_lim = ylim;
line([lna lna],[y_lim(1) y_lim(2)],'Color','white');
title('H(unit|eeg)')

subplot(9,1,9)
imagesc([x1 x2],[pwind1 pwind2],aHycx)
b_rescaleaxis('Y',f)
y_lim = ylim;
line([lna lna],[y_lim(1) y_lim(2)],'Color','white');
title('H(eeg|unit)')
set(gcf,'Position',[1 31 1280 927])     % maximize figure window

% Save
fnts = [fnm '_ENTROPY2image'];
saveas(gcf,fnts,'jpg');
close all



% -------------------------------------------------------------------------
function [hx,hy,jh,Hx,Hy,Hxy,Hxcy,Hycx,Ixy,Uxy,Uyx,Bias_HR,Bias_I] = ...
    entry_sub1(tfv1,tfv2,y1,h1,y2,h2)
%ENTRY_SUB1   Wavelet entropy calculation.
%   [hx,hy,jh,Hx,Hy,Hxy,Hxcy,Hycx,Ixy,Uxy,Uyx] = ENTRY_SUB1(tfv1,tfv2,y1,h1,y2,h2)
%   needs six inputs:
%       tvf1: eeg wavelet "row"
%       tvf2: unit wavelet "row"
%       y1: eeg wavelet "cell"
%       y2: unit wavelet "cell"
%       h1: number of bins for eeg
%       h2: number of bins for unit
%   It produces 13 outputs:
%       hx: normalized distribution of unit
%       hy: normalized distribution of eeg
%       jh: common distribution of eeg and unit
%       Hx: unit wavelet magnitude entropy
%       Hy: eeg wavelet magnitude entropy
%       Hxy: combined entropy
%       Hxcy: conditional entropy (H(unit|eeg))
%       Hycx: conditional entropy (H(eeg|unit))
%       Ixy: mutual information
%       Uxy: uncertainity coefficient (eeg->unit)
%       Uyx: uncertainity coefficient (unit->eeg)
%       Bias_HR: bias of entropy estimation
%       Bias_I: bias of mutual information estimation
%
%   ENTROPY_SUB1 uses zeros when deviding with zero entropy values.
%
%   See also ENTRY_SUB2, ENTRY and ENTRY_BETA1.

% Input argument check
error(nargchk(6,6,nargin))

% Histogram estimation
n1 = numel(y1);    % size of wavelet "cell"
n2 = numel(y2);

miny1 = min(min(tfv1));     % minimum of eeg wavelet "row"
maxy1 = max(max(tfv1));     % maximum of eeg wavelet "row"
binwidth1 = (maxy1 - miny1) ./ h1;
xx1 = miny1 + binwidth1 * (0:h1);   % bin edges
xx1(length(xx1)) = maxy1;
xx1(1) = -inf;
x1 = xx1(1:length(xx1)-1) + binwidth1 / 2;     % bin halves

miny2 = min(min(tfv2));     % minimum of unit wavelet "row"
maxy2 = max(max(tfv2));     % maximum of unit wavelet "row"
binwidth2 = (maxy2 - miny2) ./ h2;
xx2 = miny2 + binwidth2 * (0:h2);   % bin edges
xx2(length(xx2)) = maxy2;
xx2(1) = -inf;
x2 = xx2(1:length(xx2)-1) + binwidth2 / 2;     % bin halves

nbin1 = length(xx1);
nbin2 = length(xx2);
jh = zeros(nbin1-1,nbin2-1);

ty1 = y1(:) - miny1;
ty2 = y2(:) - miny2;
p = fix((ty1-1000000*eps)/binwidth1) + 1;
q = fix((ty2-1000000*eps)/binwidth2) + 1;
for i = 1:n1
    ind1 = min(p(i),size(jh,1));
    ind2 = min(q(i),size(jh,2));
    jh(ind1,ind2) = jh(ind1,ind2) + 1;
end

% Calculation of entropies & uncertainity coefficients
[m,n] = size(jh); 
N = sum(sum(jh));
hxx = sum(jh);      % marginal distribution: unit
hyy = sum(jh');     % marginal distribution: eeg
N1 = sum(hxx);
N2 = sum(hyy);
hx = hxx / N1;      % normalized distribution: unit
hy = hyy / N2;      % normalized distribution: eeg

Hx = 0;
for i = 1:m
    if hx(i) ~= 0
        a = hxx(i) / N;     % normalization
        Hx = Hx - (a * log(a));   % UNIT ENTROPY
    end
end

Hy = 0;
for k = 1:n
    if hy(k) ~= 0
        a = hyy(k) / N;     % normalization
        Hy = Hy - (a * log(a));   % EEG ENTROPY
    end
end

Hxy = 0;
for i = 1:m
    for k = 1:n 
        if jh(i,k) ~= 0
            a = jh(i,k) / N;     % normalization
            Hxy = Hxy - a * log(a);     % COMMON ENTROPY
        end
    end
end

Hycx = Hxy - Hx;    % conditional entropy
Hxcy = Hxy - Hy;

Uyx = (Hy - Hycx) / (Hy + eps);     % uncertainity coefficient
Uxy = (Hx - Hxcy) / (Hx + eps);
Ux2y = 2 * ((Hy + Hx - Hxy) / (Hx + Hy + eps));

Ixy = Hx + Hy - Hxy;    % mutual information

% Bias correction
if ~isequal(n1,n2) | ~isequal(h1,h2)
    error('Repair bias correction!')
end

Nt = n2;   % total number of trials
Ns = hyy;
R_bar = h2;    % number of possible responses
for k = 1:size(jh,1)
    Rs_bar(k) = bayescount(Nt,jh(k,:)/N);
    Rs_bar2(k) = length(find(jh(k,:)));
end
Bias_HR = ((-1) / (2 * Nt)) * (R_bar - 1);
Bias_HRS = ((-1) / (2 * Nt)) * sum(Rs_bar-1);
Bias_I = (1 / (2 * Nt)) * (sum(Rs_bar-1) - (R_bar - 1));

Ixy = Ixy - Bias_I;
Hx = Hx - Bias_HR;
Hy = Hy - Bias_HR;
Uyx = Ixy / (Hy + eps);
Uxy = Ixy / (Hx + eps);



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

% Get global variables
global FTHETA
ftheta = FTHETA;
global FNOTH
fnoth = FNOTH;
global DATINX1_THETA
datinx1_theta = DATINX1_THETA;
global DATINX2_THETA
datinx2_theta = DATINX2_THETA;
global DATINX1_NOTH
datinx1_noth = DATINX1_NOTH;
global DATINX2_NOTH
datinx2_noth = DATINX2_NOTH;


% Main
seglen = sr * c;
maxi = floor(length(eeg)/seglen);
PDC = zeros(2,2,N,maxi);
for i = 1:maxi
    ind1 = (i - 1) * seglen + 1;
    ind2 = ind1 + seglen - 1;
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
ld = (datinx2_theta - datinx1_theta) + (datinx2_noth - datinx1_noth);
time = linspace(0,ld/10000,length(LIN1));
lnb = (datinx2_noth - datinx1_noth) / 10000;
plot(time,LIN1,'r')  % unit->eeg
hold on
plot(time,LIN2,'b')  % eeg->unit
y_lim = ylim;
line([lnb lnb],[y_lim(1) y_lim(2)],'Color','green');

% Save
fnm = ftheta(1:6);
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
x11 = 0;
x22 = floor(ld/seglen) * seglen / sr;
lna = floor((datinx2_noth-datinx1_noth)/seglen) * seglen / sr;
x1 = (seglen / 2) / sr;
x2 = (floor(ld/seglen) * seglen - (seglen / 2)) / sr;
x = [x1 x2];
ls = 10;   % last plotted frequency (10 Hz)
y = [f(1) f(ls)];
xb = linspace(x11,ld/sr,length(eeg));
figure
subplot(4,1,1)
imagesc(x,y,PDC_ue(1:ls,:),[clow chigh])
y_lim = ylim;
line([lna lna],[y_lim(1) y_lim(2)],'Color','white');
title('unit->eeg')
subplot(4,1,2)
imagesc(x,y,PDC_eu(1:ls,:),[clow chigh])
y_lim = ylim;
line([lna lna],[y_lim(1) y_lim(2)],'Color','white');
title('eeg->unit')
subplot(4,1,3)
imagesc(x,y,PDC_ue(1:ls,:)-PDC_eu(1:ls,:),[clow2 chigh])
y_lim = ylim;
line([lna lna],[y_lim(1) y_lim(2)],'Color','white');
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

lim = floor((datinx2_theta-datinx1_theta)/seglen);
PDC_eu_theta = PDC_eu(:,1:lim);
PDC_ue_theta = PDC_ue(:,1:lim);
PDC_eu_noth = PDC_eu(:,lim+2:end);
PDC_ue_noth = PDC_ue(:,lim+2:end);

PDC_eu = PDC_eu_theta;
PDC_ue = PDC_ue_theta;
fnts = [ftheta '_DTF'];
save(fnts,'PDC_eu','PDC_ue');

PDC_eu = PDC_eu_noth;
PDC_ue = PDC_ue_noth;
fnts = [fnoth '_DTF'];
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