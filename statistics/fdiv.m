function D = fdiv(P,Q,fcn)
%FDIV   f-divergence.
%   D = FDIV(P,Q,FCN) calculates the f-divergence of the two input 
%   distributions using the convex function FCN.
%
%   See also KLDIST, HDISC, BDIST, JSDIV, CHISQUAREDIV, VARDIST and 
%   HARMONICMEAN.

% Input argument check
error(nargchk(3,3,nargin))
if abs(sum(P(:))-1) > 0.00001 || abs(sum(Q(:))-1) > 0.00001
    error('Input arguments must be probability distributions.')
end
if ~isequal(size(P),size(Q))
    error('Input distributions must be of the same size.')
end
if ~(isequal(exist(fcn),5)||isequal(exist(fcn),2))
    error('Third input argument must be a valid MATLAB function.')
end

% f-divergence
P2 = P(P.*Q>0);     % restrict to the common support
Q2 = Q(P.*Q>0);
P2 = P2 / sum(P2);  % renormalize
Q2 = Q2 / sum(Q2);

str = ['D = sum(P2.*' fcn '(P2./Q2));'];
eval(str)