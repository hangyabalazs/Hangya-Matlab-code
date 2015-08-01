function D = Hdisc(P,Q)
%HDISC   Hellinger's discrimination.
%   D = HDISC(P,Q) calculates Hellinger's discrimination of the two input 
%   distributions.
%
%   See also KLDIST, JSDIV, FDIV, BDIST, CHISQUAREDIV, VARDIST and 
%   HARMONICMEAN.

% Input argument check
error(nargchk(2,2,nargin))
if abs(sum(P(:))-1) > 0.00001 || abs(sum(Q(:))-1) > 0.00001
    error('Input arguments must be probability distributions.')
end
if ~isequal(size(P),size(Q))
    error('Input distributions must be of the same size.')
end

% Hellinger's discrimination
D = 0.5 * sum((sqrt(P)-sqrt(Q)).^2);