function [q,h] = b_loghist_for_loghistrun(m,k,isi,nmax)
%LOGHIST_FOR_LOGHISTRUN Version of LOGHIST used by LOGHISTRUN and CONTRASTLOGHIST.
%  [Q,H] = LOGHIST_FOR_LOGHISTRUN has two output arguments: H contains the frequency counts and Q
%   is a 1:1:nmax aritmetic serie for plotting the histogram (where nmax is the number of bins).
%
%   See also HIST, LOGHIST, LOGHISTRUN and CONTRASTLOGHIST.

% Input arguments check
error(nargchk(4,4,nargin));

% Logaritmic histogram
figure;
isi = isi * 10000;
for n = 1:nmax,
    x(n) = m*(1 + k)^(n-1);
end;
[h,xout] = hist(isi,x);
q(1:nmax) = [1:nmax];
%b = bar(q,h);