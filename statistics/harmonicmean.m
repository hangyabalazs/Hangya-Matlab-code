function D = harmonicmean(P,Q)
%HARMONICMEAN   Harmonic mean for distributions.
%   D = HARMONICMEAN(P,Q) calculates the harmonic mean of the two input 
%   distributions.
%
%   See also KLDIST, HDISC, BDIST, CHISQUAREDIV, VARDIST and JSDIV.

% Input argument check
error(nargchk(2,2,nargin))
if abs(sum(P(:))-1) > 0.00001 || abs(sum(Q(:))-1) > 0.00001
    error('Input arguments must be probability distributions.')
end
if ~isequal(size(P),size(Q))
    error('Input distributions must be of the same size.')
end

% Harmonic mean
D = sum((2*P.*Q)./(P+Q));