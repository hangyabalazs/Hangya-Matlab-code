function M = b_man_nonnan(X)
%MAX_NONNAN    Maximum for non-NaN elements.
%   M = MAX_NONNAN(X) calculates maximum for non-NaN elements of X.
%
%   See also MAX.

Y = X(find(~isnan(X)));
M = max(Y);