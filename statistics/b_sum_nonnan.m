function M = b_sum_nonnan(X)
%SUM_NONNAN    Summation for non-NaN elements.
%   M = SUM_NONNAN(X) calculates summation for non-NaN elements of X.
%
%   See also SUM.

Y = X(find(~isnan(X)));
M = sum(Y);