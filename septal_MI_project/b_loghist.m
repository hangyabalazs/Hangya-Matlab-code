function loghist(m,k,isi,nmax)
%LOGHIST    Logaritmic interspike interval histogram.
%   LOGHIST(M,K,ISI,NMAX) computes logaritmic interspike interval histogram where M and K are constants
%   used for computing the bin borders: M determines the base of the logaritm and K determines the center
%   of the first bin. NMAX is the number of the bins and ISI is the matrix containing the interspike
%   intervals given in number of data points (the measure that DISC uses).
%
%   See also HIST, LOGHISTRUN, CONTRASTLOGHIST and DISC.

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
bar(q,h);