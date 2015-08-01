function D = Bdist(P,Q)
%BDIST   Bhattacharyya-distance.
%   D = BDIST(P,Q) calculates the Bhattacharyya-distance of the two input 
%   distributions.
%
%   See also KLDIST, HDISC, FDIV, JSDIV, CHISQUAREDIV, VARDIST, 
%   BDIST_NONNAN and HARMONICMEAN.

% Input argument check
error(nargchk(2,2,nargin))
if abs(sum(P(:))-1) > 0.00001 || abs(sum(Q(:))-1) > 0.00001
    error('Input arguments must be probability distributions.')
end
if ~isequal(size(P),size(Q))
    error('Input distributions must be of the same size.')
end

% B-distance
D = sum(sqrt(P.*Q));