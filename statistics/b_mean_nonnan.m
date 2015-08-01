function M = b_mean_nonnan(X)
%MEAN_NONNAN    Mean for non-NaN elements.
%   M = MEAN_NONNAN(X) calculates mean for non-NaN elements of X.
%
%   See also MEAN and STD_NONNAN.

Y = X(find(~isnan(X)));
M = mean(Y);