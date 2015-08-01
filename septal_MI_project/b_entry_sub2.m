function [hx,hy,jh,Hx,Hy,Hxy,Hxcy,Hycx,Ixy,Uxy,Uyx] = b_entry_sub2(tfv1,tfv2,y1,h1,y2,h2)
%ENTRY_SUB2   Wavelet entropy calculation.
%   [hx,hy,jh,Hx,Hy,Hxy,Hxcy,Hycx,Ixy,Uxy,Uyx] = ENTRY_SUB2(tfv1,tfv2,y1,h1,y2,h2)
%   needs six inputs:
%       tvf1: eeg wavelet "row"
%       tvf2: unit wavelet "row"
%       y1: eeg wavelet "cell"
%       y2: unit wavelet "cell"
%       h1: number of bins for eeg
%       h2: number of bins for unit
%   It produces 11 outputs:
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
%
%   ENTROPY_SUB2 uses NaNs when deviding with zero entropy values.
%
%   See also ENTRY_SUB1, ENTRY and ENTRY_BETA2.

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
    jh(p(i),q(i)) = jh(p(i),q(i)) + 1;
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

warning off
Uyx = (Hy - Hycx) / Hy;     % uncertainity coefficient
Uxy = (Hx - Hxcy) / Hx;
Ux2y = 2 * ((Hy + Hx - Hxy) / Hx + Hy);
warning backtrace

Ixy = Hx + Hy - Hxy;    % mutual information