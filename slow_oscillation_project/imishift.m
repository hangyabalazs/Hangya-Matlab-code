function [MaxLocIxy MaxIxy] = imishift(eeg1,eeg2,sr)
%IMISHIFT   Time-shifted mutual information.
%   [ML M] = IMISHIFT(EEG1,EEG2,SR) calculates mutual information (MI) on
%   1000-point non-overlapping raw EEG traces sampled on SR. MI is
%   calculated when shifting EEG2 by 0 - 1000 ms forward. Maximal MI in the
%   function of the time shift is returned in M, maximum location is
%   returned in ML.
%   Note, that input EEG should be downsampled on 1000 Hz!
%
%   Panzeri-Treves bias correction is implemented for MI calculation.
%
%   Reference: Panzeri S, Senatore R, Montemurro MA, Petersen RS (2007)
%   Correcting for the sampling bias problem in spike train information
%   measures. Journal of neurophysiology 98:1064-1072
%
%   See also IMIRAW.

% Input argument check
error(nargchk(3,3,nargin))

% MI significance levels based on shuffled segment EEG controls
siglev95 = 1.5958;      % 95% significance level
siglev99 = 1.9121;      % 99% significance level

% Entropy calculation
leneeg = length(eeg1);
winlen = 1 * sr;   % window size: 1 sec.
maxi = floor(leneeg/winlen);
bno = fix(exp(0.626+0.4*log(winlen-1)));   % number of bins for distribution estimation
ovlp = 1;       % non-overlapping windows
maxii = maxi * ovlp - ovlp + 1;
miny1 = min(eeg1);
maxy1 = max(eeg1);
miny2 = min(eeg2);
maxy2 = max(eeg2);
MaxIxy = zeros(1,maxii);
MaxLocIxy = zeros(1,maxii);
for i = 1:maxii-1
    inx1 = (i - 1) * winlen / ovlp + 1;  % Note: overlaping windows if ovlp ~= 1
    inx1 = round(inx1);
    inx2 = inx1 + winlen - 1;
    
    T = [0:sr/1000:sr];
    Ixy = zeros(size(T));
    for t = T
        data1 = eeg1(inx1:inx2);   % eeg1 segment
        data2 = eeg2(inx1+t:inx2+t);   % shifted eeg2 segment
        Ixy(t+1) = mi(data1,data2,miny1,maxy1,miny2,maxy2,bno); % MAIN
    end
    MaxIxy(i) = max(Ixy);
    fnm = find(Ixy==MaxIxy(i));
    MaxLocIxy(i) = fnm(1);
%     if MaxIxy(i) > siglev99
%         MaxLocIxy(i) = find(Ixy==MaxIxy(i));
%     else
%         MaxIxy(i) = NaN;
%         MaxLocIxy(i) = NaN;
%     end
end



% -------------------------------------------------------------------------
function Ixy = mi(y1,y2,miny1,maxy1,miny2,maxy2,bno)

% Histogram estimation
n1 = length(y1);    % length of eeg segment
n2 = length(y2);

binwidth1 = (maxy1 - miny1) ./ bno;
xx1 = miny1 + binwidth1 * (0:bno);   % bin edges
xx1(length(xx1)) = maxy1;
xx1(1) = -inf;
binwidth2 = (maxy2 - miny2) ./ bno;
xx2 = miny2 + binwidth2 * (0:bno);   % bin edges
xx2(length(xx2)) = maxy2;
xx2(1) = -inf;
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

pHx = hxx(hxx~=0) / N;     % normalization
Hx = -sum(pHx.*log2(pHx));   % EEG2 ENTROPY

pHy = hyy(hyy~=0) / N;     % normalization
Hy = -sum(pHy.*log2(pHy));   % EEG1 ENTROPY

pHxy = jh(jh~=0) / N;     % normalization
Hxy = -sum(pHxy.*log2(pHxy));     % COMMON ENTROPY

Ixy = Hx + Hy - Hxy;    % mutual information

% Bias correction
if ~isequal(n1,n2)
    error('Repair bias correction!')
end
Nt = n2;   % total number of trials
R_bar = bno;    % number of possible responses
Rs_bar = zeros(1,size(jh,1));
for k = 1:size(jh,1)
    Rs_bar(k) = bayescount(Nt,jh(k,:)/N);
end
Bias_I = (1 / (2 * Nt * log(2))) * (sum(Rs_bar-1) - (R_bar - 1));
Ixy = Ixy - Bias_I;