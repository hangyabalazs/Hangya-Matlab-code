function psvd_new = randpoisson(lvd,leneeg)
%RANDPOISSON   Poisson process.
%   PN = RANDPOISSON(LVD,LEEG) creates frequency-adjusted Poisson spike
%   train. Input arguments:
%       LVD: number of spikes.
%       LEEG: segment length.
%       SR: sampling rate
%
%   See also RANDOM.

% Create random unit (frequency-adjusted Poisson process)
frq = lvd / leneeg;
lambda = frq;

r = random('exp',1/lambda,1,150000);     % 1/lambda param. exp.!!!! (MATLAB-ban forditva van...)
s = cumsum(r);
psvd = unique(ceil(s));     % 'pseudo vdisc'
psvd_new = psvd(psvd<leneeg);   % expected value of length(psvd): leneeg*lambda
if isequal(length(psvd),length(psvd_new))
    error('Technical error 19')
end