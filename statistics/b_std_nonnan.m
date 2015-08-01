function S = b_std_nonnan(X)
%STD_NONNAN    Standard deviation for non-NaN elements.
%   S = STD_NONNAN(X) calculates standard deviation for non-NaN elements of X.
%
%   See also STD and MEAN_NONNAN.

Y = X(find(~isnan(X)));
S = std(Y);