function out = b_isconstant(x)
%ISCONSTANT   True for constant matrix.
%   ISCONSTANT(X) returns 1 if X is constant and 0 otherwise.

y = x - x(1);
out = isempty(find(y,1));