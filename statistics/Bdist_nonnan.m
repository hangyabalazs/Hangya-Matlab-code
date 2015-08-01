function D = Bdist_nonnan(P,Q)
%BDIST:NONNAN   Bhattacharyya-distance.
%   D = BDIST_NONNAN(P,Q) calculates the Bhattacharyya-distance of the two 
%   input distributions ignoring NaN values.
%
%   See also KLDIST, HDISC, FDIV, BDIST, JSDIV, CHISQUAREDIV, VARDIST and 
%   HARMONICMEAN.

% Input argument check
error(nargchk(2,2,nargin))
if abs(sum(P(:))-1) > 0.00001 || abs(sum(Q(:))-1) > 0.00001
    error('Input arguments must be probability distributions.')
end
if ~isequal(size(P),size(Q))
    error('Input distributions must be of the same size.')
end

% B-distance
D = 0;
for x = 1:size(P,1)
    for y = 1:size(P,2)
        if ~isnan(P(x,y)*Q(x,y))
            D = D + sqrt(P(x,y)*Q(x,y));
        end
    end
end