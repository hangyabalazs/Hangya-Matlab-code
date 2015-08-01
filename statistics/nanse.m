function SE = nanse(X)
%NANSE   Standars error of mean, ignoring NaNs.
%   SE = NANSE(X) returns the sample standard error of mean of X, treating
%   NaNs as missing values.

nanlen = sum(~isnan(X));
SE = nanstd(X) ./ sqrt(nanlen);