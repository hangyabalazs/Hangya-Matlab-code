function b_2dtf
%2DTF   Partial Directed Coherence and Directed Transfer Function.
%   2DTF calculates and plots bivariate PDC and DTF for septo-hippocampal
%   data files. It calls IN3 for import and DISC for unit discrimination.
%   Default settings: 
%       segmnent length (seglen): 20000 data points
%       number of frequencies (N): 512
%       sampling rate: 10000 Hz (after downsampling (FS): 1000 Hz).
%
%   The two data set, i.e. EEG and instantenous frequency is downsampled on
%   1000 Hz, then maximal value is normalized to one.
%
%   Output: PDC and DTF plots for 2 sec. long segments (comment out when
%               working on long registrations!)
%           Integrals for PDC in theta band in the function of time
%               ('LIN'). unit->eeg: red; eeg->unit: blue.
%
%   The script is based on Baccalá and Sameshima's paper from 2001 
%   published in Biological Cybernetics.
%
%   See also ARFIT and PLOTA.

% The proram contains modified parts of the free softwares PLOTA and ARFIT.
% Modification date: 2005/10/06.

% PLOTA
% $Revision: 1.36 $
% $Id: plota.m,v 1.36 2005/09/10 20:58:23 schloegl Exp $
% Copyright (C) 1999-2004 by Alois Schloegl <a.schloegl@ieee.org>

% ARFIT
% Modified 14-Oct-00
% Authors:
% Tapio Schneider
% tapio@gps.caltech.edu
% Arnold Neumaier
% neum@cma.univie.ac.at

% Set number of frequencies
global N
N = 512;
sr = 10000;     % sampling rate

% Import and discrimination
b_in3
b_disc
b_var2ws('all','caller')

% Initialize global LIN
global LIN1
global LIN2
LIN1 = [];
LIN2 = [];

% Main
seglen = 20000;
maxi = floor(length(eeg)/seglen);
PDC = zeros(2,2,N,maxi);
for i = 1:maxi
    ind1 = (i - 1) * seglen + 1;
    ind2 = ind1 + seglen -1;
    eeg2 = eeg(ind1:ind2);
    unit2 = unit(ind1:ind2);
    vdisc2 = vdisc(find(vdisc>=ind1&vdisc<=ind2)) - ind1;
    
    if ~isempty(vdisc2)
        [PDC(:,:,:,i),f] = dircoh(eeg2,vdisc2,sr);  % partial directed coherence
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
time = [0:2:length(eeg)/10000];
plot(time(1:end-1),LIN1,'r')  % unit->eeg
hold on
plot(time(1:end-1),LIN2,'b')  % eeg->unit
% plot(LIN1-LIN2,'g')



% -------------------------------------------------------------------------
function [PDC,f] = dircoh(eeg,vdisc,sr)

% Get number of frequencies
global N

% Instantanous frequency
if vdisc(1) == 0
    vdisc = vdisc(2:end);
end
zint = ifreq(vdisc,length(eeg));
zint = zint(1:10:end);      % downsampling instantenous frequency

% Estimate optimal model order
en = eeg(1:10:end)';        % downsampling EEG
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
p = arorder(s,1,N);

% Calculate DTF and PDC
global N;    % number of frequencies
Fs = sr;  % sampling rate
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
for k = 1:K1
    Label{k} = sprintf('#%02i',k);
end

% figure
for k1 = 1:K1
    for k2 = 1:K2
%         subplot(K1,K2,k2+(k1-1)*K1);
%         area(f,squeeze(PDC(k1,k2,:)));
%         axis([0,max(f),0,1]);
%         axis([0,100,0,1]);
        sqPDC = squeeze(PDC(k1,k2,:));
%         tt1 = num2str(sum(sqPDC));      % integral
%         tt2 = num2str(max(sqPDC));      % maximum
%         text(350,0.9,tt1)
%         text(350,0.6,tt2)
%         if k2 == 1
%             ylabel(Label{k1});
%         end
%         if k1 == 1
%             title(Label{k2});
%         end
        if k1 == 1 & k2 == 2
            LIN1(end+1) = sum(sqPDC(3:6));
        end
        if k1 == 2 & k2 == 1
            LIN2(end+1) = sum(sqPDC(3:6));
        end
    end
end

% figure
% for k1 = 1:K1
%     for k2 = 1:K2
%         subplot(K1,K2,k2+(k1-1)*K1);
%         area(f,squeeze(DTF(k1,k2,:)));
%         axis([0,max(f),0,1]);
% %         axis([0,100,0,1]);
%         if k2 == 1;
%             ylabel(Label{k1});
%         end
%         if k1 == 1;
%             title(Label{k2});
%         end
%     end
% end