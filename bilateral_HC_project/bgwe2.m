function bgwe2(data,trials)
%BGWE2   Wavelet entropy on Buzsaki's EEG files.
%   BGWE2 calculates wavelet, FFT and filtered signals for low and high
%   gamma band and wavelet entropy.
%
%   See also BGWE.

% Results directory
global DATAPATH
resdir = [DATAPATH 'Wheel\WE\'];
mm = pwd;
cd(resdir)

% Get trial EEG
lt = length(trials);
% lt = 1;
error_WE = [];
futureL_WE = [];
futureR_WE = [];
for k = 1:lt
    tS = round(trials{k}.tStart*4/5);  % indices are transformed to match new samp. rate
    tE = round(trials{k}.tEnd*4/5);
    eeg = data(10,tS:tE);   % shank: 4 or 10
    
% Filter
    sr = 1000;
    nqf = sr / 2;
    flt = fir1(1024,[30 40]/nqf,'band');      % bandpass filtering
    feeg_low = filtfilt(flt,1,eeg);
    feeg_low = (feeg_low - mean(feeg_low)) / std(feeg_low);
    flt = fir1(1024,[80 100]/nqf,'band');      % bandpass filtering
    feeg_high = filtfilt(flt,1,eeg);
    feeg_high = (feeg_high - mean(feeg_high)) / std(feeg_high);
    
% Wavelet
    [wavea_abs wavea_phase f spf] = eeg_wavelet(eeg); % eeg sampled on 1000 Hz
    H1 = figure;
    imagesc(wavea_abs)
    b_rescaleaxis('Y',f)
        
    fnd = find(f>30);    % frequency band bounderies
    pwind1 = fnd(end);
    fnd = find(f<100);
    pwind2 = fnd(1);
    H2 = figure;
    S1 = subplot(3,1,1);
    imagesc(wavea_abs(1:pwind1,:))
    b_rescaleaxis('Y',f)
    S2 = subplot(3,1,2);
    plot(eeg)
    hold on
    plot(feeg_low*1000,'Color',[204 204 204]/255)
    x_lim = xlim;
    xlim([x_lim(1) x_lim(1)+length(eeg)])
    S2 = subplot(3,1,3);
    plot(eeg)
    hold on
    plot(feeg_high*1000,'Color',[204 204 204]/255)
    x_lim = xlim;
    xlim([x_lim(1) x_lim(1)+length(eeg)])
       
% Wavelet entropy
    aHy = entropy_image(wavea_abs,f,1000,4,spf);
    H3 = figure;
    imagesc(aHy)
    
% FFT
    eeg2 = eeg(1:5:end);        % downsample on 200 Hz
    [y,w] = b_fft2(eeg2,200);     % FFT
    inxs_low = w > 30 & w < 40;
    inxs_high = w > 80 & w < 100;
    inxs2 = w > 5 & w < 100;
    integrated_power_low = sum(y(inxs_low));
    integrated_power_high = sum(y(inxs_high));
    H4 = figure;
    plot(w(inxs2),y(inxs2))
    x_lim = xlim;
    y_lim = ylim;
    text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.75+y_lim(1),...
        ['30-40 Hz: ' num2str(integrated_power_low/1000000) ' * 10^6'])
    text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.70+y_lim(1),...
        ['80-100 Hz: ' num2str(integrated_power_high/1000000) ' * 10^6'])
    
% Save
    if trials{k}.errorTrial
        cd rightHC_error
        error_WE = [error_WE sum(aHy,2)];
    else
        if isequal(trials{k}.futureLR,'L')
            cd rightHC_futureL
            futureL_WE = [futureL_WE sum(aHy,2)];
        elseif isequal(trials{k}.futureLR,'R')
            cd rightHC_futureR
            futureR_WE = [futureR_WE sum(aHy,2)];
        else
            keyboard
        end
    end
    fn = ['g01_m13_t' num2str(k) '_wavelet.jpg'];
    saveas(H1,fn)
    fn = ['g01_m13_t' num2str(k) '_gamma.jpg'];
    saveas(H2,fn)
    fn = ['g01_m13_t' num2str(k) '_WE.jpg'];
    saveas(H3,fn)
    fn = ['g01_m13_t' num2str(k) '_FFT.fig'];
    saveas(H4,fn)
    cd ..
    close all
end

% Plot and save grouped wavelet entropy results
H5 = figure;
imagesc(futureR_WE)
fn = ['g01_m13_futureR_WE.jpg'];
saveas(H5,fn)

H6 = figure;
imagesc(futureL_WE)
fn = ['g01_m13_futureL_WE.jpg'];
saveas(H6,fn)

H7 = figure;
imagesc(error_WE)
fn = ['g01_m13_error_WE.jpg'];
saveas(H7,fn)
save('g01_m13_grouped_WE','error_WE','futureR_WE','futureL_WE')
cd(mm)


% -------------------------------------------------------------------------
% WAVELET
% -------------------------------------------------------------------------
function [pow,phase,f,spf] = eeg_wavelet(dat)

% Prepare for wavelet transformation
variance = std(dat) ^ 2;
dat = (dat - mean(dat)) / sqrt(variance) ;
dt = 1 / 1000;
pad = 1;
spf = 5;
omega0 = 6;
c = 4 * pi / (omega0 + sqrt(2+omega0^2));
s = 1 ./ (120:-(1/spf):1) ./ c;
fperiod = c .* s;
f = 1 ./ fperiod;
param = -1;
mother = 'Morlet';

% Wavelet transformation
[wave,period,scale,coi] = wavelet(dat,dt,pad,s,mother,param);
pow = abs(wave).^2;
phase = angle(wave);

% -------------------------------------------------------------------------
function [wave,period,scale,coi] = wavelet(Y,dt,pad,s,mother,param)
%WAVELET   Wavelet with scales for every frequency.
%
%   Copyright (C) 1995-1998, Christopher Torrence and Gilbert P. Compo
%   University of Colorado, Program in Atmospheric and Oceanic Sciences.
%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made.  This
%   routine is provided as is without any express or implied warranties
%   whatsoever.

n1 = length(Y);

% Construct time series to analyze, pad if necessary
x(1:n1) = Y - mean(Y);
if (pad == 1)
	base2 = fix(log(n1)/log(2)+0.4999);   % power of 2 nearest to N
	x = [x,zeros(1,2^(base2+1)-n1)];
end
n = length(x);

% Construct wavenumber array used in transform [Eqn(5)]
k = 1:fix(n/2);
k = k .* ((2 .* pi) / (n * dt));
k = [0., k, -k(fix((n-1)/2):-1:1)];

% Compute FFT of the (padded) time series
f = fft(x);    % [Eqn(3)]

% Construct SCALE array & empty PERIOD & WAVE arrays
scale = s;

% Loop through all scales and compute transform
wave = zeros(length(scale),n1);  % define the wavelet array
for a1 = 1:length(scale)
	[daughter,fourier_factor,coi,dofmin] = wave_bases(mother,k,scale(a1),param);	
	ifd = ifft(f.*daughter);  % wavelet transform[Eqn(4)]
    ifd = ifd(1:n1);
    wave(a1,:) = ifd;
end

period = fourier_factor * scale;
coi = coi * dt * [1E-5,1:((n1+1)/2-1),fliplr((1:(n1/2-1))),1E-5];  % COI [Sec.3g]



% -------------------------------------------------------------------------
% ENTROPY
% -------------------------------------------------------------------------
function aHy = entropy_image(data,f,WindowSize,Overlap,spf)

[k1 k2] = size(data);

winlen = WindowSize;   % window size
maxi = floor(k2/winlen);

aHy = [];     % preallocating output variables

% Entropy calculation
ovlp = Overlap;
for j = 1:k1/spf
    tHy = [];   % preallocating temporal variables
    for i = 1:maxi*ovlp-ovlp+1        % ABS LOOP
        inx1 = (i - 1) * winlen / ovlp + 1;  % Note: overlaping windows!
        inx1 = round(inx1);
        inx2 = inx1 + winlen - 1;
        inxx1 = (j - 1) * spf + 1;
        inxx2 = inxx1 + spf - 1;
        y = data(inxx1:inxx2,inx1:inx2);   % eeg wavelet "cells"
        numy = numel(y);
        bno = fix(exp(0.626+0.4*log(numy-1)));   % number of bins
        Hy = entry_sub1(data(inxx1:inxx2,:),y,bno); % MAIN
        tHy = [tHy Hy];
    end
    aHy(j,:) = tHy;
end     % end of abs loop

% -------------------------------------------------------------------------
function [Hy,Bias_HR] = entry_sub1(tfv,y,h)
%ENTRY_SUB1   Wavelet entropy calculation.
%   [hx,Bias_HR] = ENTRY_SUB1(tfv,y,h)
%   needs three inputs:
%       tfv: eeg wavelet "row"
%       y: eeg wavelet "cell"
%       h: number of bins for eeg
%   It produces two outputs:
%       Hy: eeg wavelet magnitude entropy
%       Bias_HR: bias of entropy estimation
%
%   ENTROPY_SUB1 uses zeros when deviding with zero entropy values.
%
%   See also ENTRY_SUB2, ENTRY and ENTRY_BETA1.

% Input argument check
error(nargchk(3,3,nargin))

% Histogram estimation
n = numel(y);    % size of wavelet "cell"

miny = min(min(tfv));     % minimum of eeg wavelet "row"
maxy = max(max(tfv));     % maximum of eeg wavelet "row"
binwidth = (maxy - miny) ./ h;
xx = miny + binwidth * (0:h);   % bin edges
xx(length(xx)) = maxy;
xx(1) = -inf;
x = xx(1:length(xx)-1) + binwidth / 2;     % bin halves
nbin = length(xx);

ty = y(:) - miny;
p = fix((ty-1000000*eps)/binwidth) + 1;
hyy = zeros(1,nbin-1);     % EEG distribution
for i = 1:n
    ind = min(p(i),length(hyy));
    hyy(1,ind) = hyy(1,ind) + 1;
end

% Calculation of entropies & uncertainity coefficients
N = sum(hyy);
hy = hyy / N;      % normalized distribution: eeg

Hy = 0;
for k = 1:nbin-1
    if hy(k) ~= 0
        a = hyy(k) / N;     % normalization
        Hy = Hy - (a * log2(a));   % EEG ENTROPY
    end
end

% Bias correction
Nt = n;   % total number of trials
R_bar = h;    % number of possible responses
Bias_HR = ((-1) / (2 * Nt * log(2))) * (R_bar - 1);
Hy = Hy - Bias_HR;

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
function [hx,hy,jh,Hx,Hy,Hxy,Hxcy,Hycx,Ixy,Uxy,Uyx,Bias_HR,Bias_I] = ...
    entry_sub2(tfv1,tfv2,y1,h1,y2,h2)
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