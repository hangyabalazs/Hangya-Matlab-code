function C = czcmpl(P,Q,c)
%CZCMPL   Place field complementarity index.
%   C = CZCMPL(P,Q,L) calculates complementarity between place fields (not
%   rate maps! - see CZPLACEANALYSIS) P and Q, when L is the size of the
%   arena in pixels. Complementarity index is defined as the percentage of
%   the arena covered by the symmetric difference of the two place fields.
%
%   See also CZPLACEANALYSIS.

% Input argument check
error(nargchk(3,3,nargin))
if abs(sum(P(:))-1) > 0.00001 || abs(sum(Q(:))-1) > 0.00001
    error('Input arguments must be probability distributions.')
end
if ~isequal(size(P),size(Q))
    error('Input distributions must be of the same size.')
end

% Complementarity index
C = length(find((isnan(P)&~isnan(Q))|(~isnan(P)&isnan(Q)))) / c;