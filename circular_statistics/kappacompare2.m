function [p, H] = kappacompare2(ng1,ng2,degrad,str)
%KAPPACOMPARE2   Test for equality of concentration parameters of two circular distributions.
%   [P,H] = KAPPACOMPARE2(NG1,NG2,DEGRAD) tests the null hypothesis kappa1
%   = kappa2 against the alternative of different concentration parameters
%   of the von Mises samples NG1 and NG2. The third input argument (DEGRAD)
%   should contain the information whether the sample is measured in
%   degress ('deg') or radians ('rad'). p-value (P) and indicator of
%   rejection of the null hypothesis (H) are returned. A randomisation test
%   is applied with the test statistic kappa1 - kappa2 in the onesided and
%   abs(kappa1-kappa2) in the twosided case. The test is exact incase the 
%   number of all sample combinations is below 10000, otherwise a Monte
%   Carlo randomisation using 10000 random permutations is applied.
%
%   [P,H] = KAPPACOMPARE2(NG1,NG2,DEGRAD,'onesided') evaluates the null
%   hypothesis kappa1 <= kappa2 against the alternative kappa1 > kappa2.
%
%   See also WATSON and KAPPACOMPARE.

% Input argument check
error(nargchk(3,4,nargin))
if nargin == 3
    str = 'twosided';
elseif nargin == 4
    if ~isequal(str,'onesided') && ~isequal(str,'twosided')
        error('Forth input argument should be ''onesided'' or ''twosided''.')
    end
end
if ~isequal(degrad,'deg') && ~isequal(degrad,'rad')
    error('Third input argument should be ''deg'' or ''rad''.')
end
alpha = 0.05;

% Convert to radians
if isequal(degrad,'deg')
    ng1 = ng1 / 180 * pi;
    ng2 = ng2 / 180 * pi;
end

% Performing Monte Carlo randomisation test
p = randomisationtest(ng1,ng2,str);

% Output  
if p <= alpha
    H = 1;
else
    H = 0;
end

% -------------------------------------------------------------------------
function T = teststat(ng1,ng2,str)
% Test statistic.

[kappa_ML1 kappa1] = kappaest(ng1,'rad');
[kappa_ML2 kappa2] = kappaest(ng2,'rad');
if isequal(str,'onesided')
    T = kappa1 - kappa2;
else
    T = abs(kappa1-kappa2);
end

% -------------------------------------------------------------------------
function p = randomisationtest(ng1,ng2,str)

png = [ng1 ng2];      % pooled sample
png = png(:);
n1 = length(ng1);
mno = 10000;     % maximal number of test statistics calculation
warning off
nk = nchoosek(length(png),n1);
warning backtrace
if nk < mno
    N = nk;
    allcomb = nchoosek(png,n1);     % all possible combinations
else
    N = mno;
    allcomb = zeros(mno,n1);
    rand('twister', sum(100*fliplr(clock)));    % initialize the state of the random generator
    for k = 1:mno
        rp = randperm(length(png));
        allcomb(k,:) = png(rp(1:n1));
    end
    disp('Randomisation distribution is estimated.')
end
frs = zeros(1,N);
for k = 1:N
    tng1 = allcomb(k,:);
    tng2 = setdiff(png,tng1);
    frs(k) = teststat(tng1,tng2,str);
end
fr = teststat(ng1,ng2,str);
sfrs = sort(frs,'ascend');
if fr > sfrs(end)
    m = mno;
else
    m2 = find(sfrs>fr,1,'first');
    m1 = m2 - 1;
    m = m1 + (fr - sfrs(m1)) / (sfrs(m2) - sfrs(m1));
end
p = (N - m + 1) / N;