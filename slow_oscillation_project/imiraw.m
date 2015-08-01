function [rHx,rHy,rHxy,rIxy,rUxy,rUyx,cHx,cHy,cHxy,cIxy,cUxy,cUyx] = imiraw(eeg1,eeg2,sr,overlap)
%IMIRAW   Runs entropy on a sequence of recordings.
%   IMIRAW calculates entropy measures on 1000-point raw EEG traces. It
%   returns entropy, mutual information and uncertainty coefficient values
%   for real data and shuffled segment controls. Note, that input EEG
%   should be downsampled on 1000 Hz!
%
%   IMIRAW(EEG1,EEG2,SR,OLP) calculates mutual information between EEG1 and
%   EEG2 sampled on SR using overlapping windows. Overlap is controlled by
%   OLP (1 for non-overlapping, 2 for 50% overlapping, etc.).
%
%   Panzeri-Treves bias correction is implemented for entropy and mutual
%   information calculation.
%
%   Reference: Panzeri S, Senatore R, Montemurro MA, Petersen RS (2007)
%   Correcting for the sampling bias problem in spike train information
%   measures. Journal of neurophysiology 98:1064-1072
%
%   See also IENTRY.

% Input argument check
error(nargchk(4,4,nargin))

% Create random eeg (segment shuffle)
seglen = length(eeg2) * sr / 1000;
ic1 = fix(seglen/sr);
increm = round(seglen/ic1);     % 1-1.2 sec. long segments
ed = {};
for t = 1:ic1
    ind1 = (t - 1) * increm + 1;
    ind2 = ind1 + increm - 1;
    ed{end+1} = eeg2(ind1:ind2);
end
led = length(ed);
rp = randperm(led);
while any(rp==[1:led])
    rp = randperm(led);
end
psed = [];
for j = 1:led
    psed = [psed ed{rp(j)}];
end

% Entropy
elen = min(length(psed),seglen);
[rHx,rHy,rHxy,rIxy,rUxy,rUyx] = entropy_line(eeg1(1:elen),eeg2(1:elen),1000,overlap);    % real
[cHx,cHy,cHxy,cIxy,cUxy,cUyx] = entropy_line(eeg1(1:elen),psed(1:elen),1000,overlap);    % control



% -------------------------------------------------------------------------
% ENTROPY
% -------------------------------------------------------------------------
function [aHx,aHy,aHxy,aIxy,aUxy,aUyx] = entropy_line(data1,data2,WindowSize,Overlap)

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

    y1 = data1(inx1:inx2);   % eeg1 segment
    y2 = data2(inx1:inx2);   % eeg2 segment
    numy = numel(y1);
    bno = fix(exp(0.626+0.4*log(numy-1)));   % number of bins
    [hx,hy,jh,Hx,Hy,Hxy,Hxcy,Hycx,Ixy,Uxy,Uyx,Bias_HR,Bias_I] =...
            entry_sub1(data1,data2,y1,bno,y2,bno); % MAIN
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
n1 = length(y1);    % length of eeg segment
n2 = length(y2);

miny1 = min(tfv1);     % minimum of eeg1
maxy1 = max(tfv1);     % maximum of eeg1
binwidth1 = (maxy1 - miny1) ./ h1;
xx1 = miny1 + binwidth1 * (0:h1);   % bin edges
xx1(length(xx1)) = maxy1;
xx1(1) = -inf;
x1 = xx1(1:length(xx1)-1) + binwidth1 / 2;     % bin halves

miny2 = min(tfv2);     % minimum of eeg2
maxy2 = max(tfv2);     % maximum of eeg2
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
hxx = sum(jh);      % marginal distribution: eeg2
hyy = sum(jh');     % marginal distribution: eeg1
N1 = sum(hxx);
N2 = sum(hyy);
hx = hxx / N1;      % normalized distribution: eeg2
hy = hyy / N2;      % normalized distribution: eeg1

Hx = 0;
for i = 1:m
    if hx(i) ~= 0
        a = hxx(i) / N;     % normalization
        Hx = Hx - (a * log2(a));   % EEG2 ENTROPY
    end
end

Hy = 0;
for k = 1:n
    if hy(k) ~= 0
        a = hyy(k) / N;     % normalization
        Hy = Hy - (a * log2(a));   % EEG1 ENTROPY
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