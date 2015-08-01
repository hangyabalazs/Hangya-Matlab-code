function [p, H] = kappacompare(ng1,ng2,degrad,str)
%KAPPACOMPARE   Test for equality of concentration parameters of two von Mises distributions.
%   [P,H] = KAPPACOMPARE(NG1,NG2,DEGRAD) tests the null hypothesis kappa1 =
%   kappa2 against the alternative of different concentration parameters of
%   the von Mises samples NG1 and NG2. The third input argument (DEGRAD)
%   should contain the information whether the sample is measured in
%   degress ('deg') or radians ('rad'). p-value (P) and indicator of
%   rejection of the null hypothesis (H) are returned. The test is due to
%   N. I. Fisher.
%
%   [P,H] = KAPPACOMPARE(NG1,NG2,DEGRAD,'onesided') evaluates the null
%   hypothesis kappa1 <= kappa2 against the alternative kappa1 > kappa2.
%
%   Reference: Fisher NI (1995) Statistical analysis of circular data.
%   Cambridge University Press, Cambridge UK, pp. 131-132.
%
%   See also WATSON and KAPPACOMPARE2.

% Input argument check
error(nargchk(3,4,nargin))
if nargin == 3
    str = 'twosided';
elseif nargin == 4
    if ~isequal(str,'onesided') && ~isequal(str,'twosided')
        error('Forth input argument should be ''onesided'' or ''twosided''.')
    end
end
if ~isequal(degead,'deg') && ~isequal(degrad,'rad')
    error('Third input argument should be ''deg'' or ''rad''.')
end
alpha = 0.05;

% Convert to radians
if isequal(degrad,'deg')
    ng1 = ng1 / 180 * pi;
    ng2 = ng2 / 180 * pi;
end

% Circular mean, mean resultant length and concentration parameter (kappa)
n1 = length(ng1);
n2 = length(ng2);
N = n1 + n2;
[ftm1, cm1, mvl1] = mvlmn(ng1,'rad');
[ftm2, cm2, mvl2] = mvlmn(ng2,'rad');
[kappa_ML1 kappa1] = kappaest(ng1,'rad');
[kappa_ML2 kappa2] = kappaest(ng2,'rad');

% Performing the test
kappa_wave = median([kappa1,kappa2]);
if n1 >= 10 && n2 >= 10 && kappa_wave >= 1      % exact test
    fr = teststat(ng1,ng2);
%     cv = finv(1-alpha,1,N-2);
    p = 1 - fcdf(fr,1,N-2);
else        % randomisation test
    if n1 < 10 || n2 < 10
        disp('Applying randomisation test due to small sample size.')
    elseif kappa_wave < 1
        disp('Applying randomisation test due to small concentration parameter.')
    end
    cng1 = ng1 - cm1;   % centering samples
    cng2 = ng2 - cm2;
    cng1(cng1<-pi) = cng1(cng1<-pi) + 2 * pi;
    cng1(cng1>pi) = cng1(cng1>pi) - 2 * pi;
    cng2(cng2<-pi) = cng2(cng2<-pi) + 2 * pi;
    cng2(cng2>pi) = cng2(cng2>pi) - 2 * pi;
    p = randomisationtest(cng1,cng2);
end
if isequal(str,'twosided')
    p = 2 * p;
end

% Output  
if p <= alpha
    H = 1;
else
    H = 0;
end

% -------------------------------------------------------------------------
function fr = teststat(ng1,ng2)
% Test statistic.

n1 = length(ng1);
n2 = length(ng2);
N = n1 + n2;
[ftm1, cm1, mvl1] = mvlmn(ng1,'rad');
[ftm2, cm2, mvl2] = mvlmn(ng2,'rad');
d1 = abs(sin(0.5*(ng1-cm1)));
d2 = abs(sin(0.5*(ng2-cm2)));
d1_bar = mean(d1);
d2_bar = mean(d2);
d_bar = (n1 * d1_bar + n2 * d2_bar) / N;
fr = ((N - 2) * (n1 * (d1_bar - d_bar)^2 + n2 * (d2_bar - d_bar)^2)) /...
    (sum((d1-d1_bar).^2) + sum((d2-d2_bar).^2));
% if fr < 1
%     fr
% end

% -------------------------------------------------------------------------
function p = randomisationtest(cng1,cng2)
%   Reference: Fisher NI (1995) Statistical analysis of circular data.
%   Cambridge University Press, Cambridge UK, pp. 131-132., 214-216.

png = [cng1 cng2];      % pooled sample
png = png(:);
n1 = length(cng1);
mno = 5000;     % maximal number of test statistics calculation
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
    frs(k) = teststat(tng1,tng2);
end
fr = teststat(cng1,cng2);
sfrs = sort(frs,'ascend');
if fr > sfrs(end)
    m = 5000;
else
    m2 = find(sfrs>fr,1,'first');
    m1 = m2 - 1;
    m = m1 + (fr - sfrs(m1)) / (sfrs(m2) - sfrs(m1));
end
p = (N - m + 1) / N;