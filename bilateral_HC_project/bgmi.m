function bgmi(data,trials)
%BGMI   Wavelet mutual information on Buzsaki's EEG files.
%   BGMI(DATA,TRIALS) calculates time-resolved wavelet MI between
%   frequencies from 1 to 119 Hz from the TRIALS of DATA. A matrix of MI
%   values are saved in the result directory.
%
%   See also BGWE2.

% Results directory
global DATAPATH
resdir = [DATAPATH 'Wheel\MI\'];
mm = pwd;
cd(resdir)

% Get trial EEG
lt = length(trials);
% lt = 1;
% error_WE = [];
% futureL_WE = [];
% futureR_WE = [];
disp('From k = 30!!!')
for k = 30:lt
    k
    tS = round(trials{k}.tStart*4/5);  % indices are transformed to match new samp. rate
    tE = round(trials{k}.tEnd*4/5);
    eeg = data(10,tS:tE);   % shank: 4 or 10
    
% Filter
    sr = 1000;
    nqf = sr / 2;
    
% Wavelet
    [wavea_abs wavea_phase f spf] = eeg_wavelet(eeg); % eeg sampled on 1000 Hz
%     H1 = figure;
%     imagesc(wavea_abs)
%     b_rescaleaxis('Y',f)
        
% Wavelet MI
    MI = wavelet_mi(wavea_abs,f,1000,4,spf);
    
% Save
    if trials{k}.errorTrial
        cd rightHC_error
%         error_WE = [error_WE sum(aHy,2)];
    else
        if isequal(trials{k}.futureLR,'L')
            cd rightHC_futureL
%             futureL_WE = [futureL_WE sum(aHy,2)];
        elseif isequal(trials{k}.futureLR,'R')
            cd rightHC_futureR
%             futureR_WE = [futureR_WE sum(aHy,2)];
        else
            keyboard
        end
    end
    fn = ['g01_m13_t' num2str(k) '_MI.mat'];
    save(fn,'MI')
    cd ..
    close all
end
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
function MI = wavelet_mi(data,f,WindowSize,Overlap,spf)

[k1 k2] = size(data);

winlen = WindowSize;   % window size
maxi = floor(k2/winlen);

aHy = [];     % preallocating output variables

% Entropy calculation
ovlp = Overlap;
MI = zeros(floor(k1/spf),floor(k1/spf),maxi*ovlp-ovlp+1);
tic
for j1 = 1:k1/spf
    for j2 = j1:k1/spf
        disp([num2str(j1) ' ' num2str(j2)])
        for i = 1:maxi*ovlp-ovlp+1        % ABS LOOP
            inx1 = (i - 1) * winlen / ovlp + 1;  % Note: overlapping windows!
            inx1 = round(inx1);
            inx2 = inx1 + winlen - 1;
            inxxj11 = (j1 - 1) * spf + 1;
            inxxj12 = inxxj11 + spf - 1;
            yj1 = data(inxxj11:inxxj12,inx1:inx2);   % eeg wavelet "cells"
            inxxj21 = (j2 - 1) * spf + 1;
            inxxj22 = inxxj21 + spf - 1;
            yj2 = data(inxxj21:inxxj22,inx1:inx2);
            numy = numel(yj1);
            bno = fix(exp(0.626+0.4*log(numy-1)));   % number of bins
            [hx,hy,jh,Hx,Hy,Hxy,Hxcy,Hycx,Ixy,Uxy,Uyx,Bias_HR,Bias_I] = ...
                entry_sub1(data(inxxj11:inxxj12,:),data(inxxj21:inxxj22,:),...
                yj1,bno,yj2,bno); % MAIN
            MI(j1,j2,i) = Ixy;
            MI(j2,j1,i) = Ixy;
        end
    end
end     % end of abs loop
toc

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
        Hx = Hx - (a * log2(a));   % UNIT ENTROPY
    end
end

Hy = 0;
for k = 1:n
    if hy(k) ~= 0
        a = hyy(k) / N;     % normalization
        Hy = Hy - (a * log2(a));   % EEG ENTROPY
    end
end

Hxy = 0;
for i = 1:m
    for k = 1:n 
        if jh(i,k) ~= 0
            a = jh(i,k) / N;     % normalization
            Hxy = Hxy - a * log2(a);     % COMMON ENTROPY
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
Bias_HR = ((-1) / (2 * Nt * log(2))) * (R_bar - 1);
Bias_HRS = ((-1) / (2 * Nt * log(2))) * sum(Rs_bar-1);
Bias_I = (1 / (2 * Nt * log(2))) * (sum(Rs_bar-1) - (R_bar - 1));

Ixy = Ixy - Bias_I;
Hx = Hx - Bias_HR;
Hy = Hy - Bias_HR;
Uyx = Ixy / (Hy + eps);
Uxy = Ixy / (Hx + eps);

% -------------------------------------------------------------------------
function [Hy,Bias_HR] = entry_sub2(tfv,y,h)
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

% Calculation of entropies
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