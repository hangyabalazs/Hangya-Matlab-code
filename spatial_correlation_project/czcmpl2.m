function [C Cb Cc] = czcmpl2(P,Q,c)
%CZCMPL2   Place field complementarity index.
%   [C CB CC] = CZCMPL(P,Q,L) calculates complementarity between place fields
%   (not rate maps! - see CZPLACEANALYSIS) P (interneuron) and Q (place
%   cell), when L is the size of the arena in pixels. Complementarity index
%   is defined as follows.
%       C = (|Q \ P| / |Q| + |P u Q| / L) / 2
%       CB = (|P u Q|) / (|P| + |Q|) + |P u Q| / L) / 2
%       CC = (|Q \ P| / |Q| + |Q \ P| / (L - |P|)) / 2
%   The first element of the addidion quantifies 'disjunctness', the second
%   stands for coverage.
%
%   See also CZPLACEANALYSIS, CZCMPL and CZPLACEFIGS.

% Input argument check
error(nargchk(3,3,nargin))
if abs(sum(P(:))-1) > 0.00001 || abs(sum(Q(:))-1) > 0.00001
    error('Input arguments must be probability distributions.')
end
if ~isequal(size(P),size(Q))
    error('Input distributions must be of the same size.')
end

% Complementarity index
C1 = length(find(isnan(P)&~isnan(Q))) / length(find(~isnan(Q)));     % disjunct
C1b = length(find(~isnan(P)|~isnan(Q))) / (length(find(~isnan(Q))) + length(find(~isnan(P))));
C2 = length(find(~isnan(P)|~isnan(Q))) / c;         % coverage
C2b = length(find(isnan(P)&~isnan(Q))) / (c - length(find(~isnan(P))));
C = (C1 + C2) / 2;
Cb = (C1b + C2) / 2;
Cc = (C1 + C2b) / 2;