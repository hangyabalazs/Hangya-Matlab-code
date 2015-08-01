function D = vardist(P,Q)
%VARDIST   Variational distance for distributions.
%   D = VARDIST(P,Q) calculates the variational distance of the two input 
%   distributions.
%
%   See also KLDIST, HDISC, BDIST, CHISQUAREDIV, JSDIV and HARMONICMEAN.

% Input argument check
error(nargchk(2,2,nargin))
if abs(sum(P(:))-1) > 0.00001 || abs(sum(Q(:))-1) > 0.00001
    error('Input arguments must be probability distributions.')
end
if ~isequal(size(P),size(Q))
    error('Input distributions must be of the same size.')
end

% Variational distance
D = sum(abs(P-Q));