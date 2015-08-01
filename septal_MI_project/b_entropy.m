function [hx,hy,jh,Hx,Hy,Hxy,Hxcy,Hycx,Ixy,Uxy,Uyx] = b_entropy(y1,y2)
%ENTROPY   Entropy calculation.
%   [hx,hy,jh,Hx,Hy,Hxy,Hxcy,Hycx,Ixy,Uxy,Uyx] = ENTROPY(y1,y2) needs 2
%   inputs:
%       y1: first data series
%       y2: second data series
%
%   It produces 11 outputs:
%       hx: normalized distribution of y1
%       hy: normalized distribution of y2
%       jh: common distribution of y1 and y2
%       Hx: y1 wavelet magnitude entropy
%       Hy: y2 wavelet magnitude entropy
%       Hxy: combined entropy
%       Hxcy: conditional entropy (H(y1|y2))
%       Hycx: conditional entropy (H(y2|y1))
%       Ixy: mutual information
%       Uxy: uncertainity coefficient (y2->y1)
%       Uyx: uncertainity coefficient (y1->y2)
%
%   See also ENTR2 and KWENTROPY.

% Input argument check
error(nargchk(2,2,nargin))

% Histogram estimation
n1 = size(y1,1) * size(y1,2);    % size of data
n2 = size(y2,1) * size(y2,2);

miny1 = min(min(y1));     % minimum of data1
maxy1 = max(max(y1));     % maximum of data1
h1 = fix(exp(0.626+0.4*log(n1-1)));   % number of bins
binwidth1 = (maxy1 - miny1) ./ h1;
xx1 = miny1 + binwidth1 * (0:h1);   % bin edges
xx1(length(xx1)) = maxy1;
xx1(1) = -inf;
x1 = xx1(1:length(xx1)-1) + binwidth1 / 2;     % bin halves


miny2 = min(min(y2));     % minimum of data2
maxy2 = max(max(y2));     % maximum of data2
h2 = fix(exp(0.626+0.4*log(n2-1)));   % number of bins
binwidth2 = (maxy2 - miny2) ./ h2;
xx2 = miny2 + binwidth2 * (0:h2);   % bin edges
xx2(length(xx2)) = maxy2;
xx2(1) = -inf;
x2 = xx2(1:length(xx2)-1) + binwidth2 / 2;     % bin halves

nbin1 = length(xx1);
nbin2 = length(xx2);
jh = zeros(nbin1-1,nbin2-1);

ty1 = y1 - miny1;
ty2 = y2 - miny2;
p = fix(ty1/binwidth1-10000*eps) + 1;
q = fix(ty2/binwidth2-10000*eps) + 1;
for i = 1:n1
    jh(p(i),q(i)) = jh(p(i),q(i)) + 1;
end

% Calculation of entropies & uncertainity coefficients
[m,n] = size(jh); 
N = sum(sum(jh));
hxx = sum(jh');      % marginal distribution: y1
hyy = sum(jh);     % marginal distribution: y2
N1 = sum(hxx);
N2 = sum(hyy);
hx = hxx / N1;      % normalized distribution: y1
hy = hyy / N2;      % normalized distribution: y2

Hx = 0;
for i = 1:m
    if hx(i) ~= 0
        a = hxx(i) / N;     % normalization
        Hx = Hx - (a * log(a));   % Y1 ENTROPY
    end
end

Hy = 0;
for k = 1:n
    if hy(k) ~= 0
        a = hyy(k) / N;     % normalization
        Hy = Hy - (a * log(a));   % Y2 ENTROPY
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