function D = chisquarediv(P,Q)
%CHISQUAREDIV   Chi-square divergence.
%   D = CHISQUAREDIV(P,Q) calculates the chi-square divergence of the two 
%   input distributions.
%
%   See also KLDIST, HDISC, BDIST, JSDIV, FDIV, VARDIST and HARMONICMEAN.

% Input argument check
error(nargchk(2,2,nargin))
if abs(sum(P(:))-1) > 0.00001 || abs(sum(Q(:))-1) > 0.00001
    error('Input arguments must be probability distributions.')
end
if ~isequal(size(P),size(Q))
    error('Input distributions must be of the same size.')
end

% KL-distance
P2 = P(P.*Q>0);     % restrict to the common support
Q2 = Q(P.*Q>0);
P2 = P2 / sum(P2);  % renormalize
Q2 = Q2 / sum(Q2);

D = sum(((P2-Q2).^2)./Q2);