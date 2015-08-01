function bootstrap_se = se_of_median(sample_dist,bno)
%SE_OF_MEDIAN   Standard error of median.
%   BOOTSTRAP_SE = SE_OF_MEDIAN(SAMPLE_DIST) calculates bootstrap standard
%   error of median (BOOTSTRAP_SE) from a bootstrap sample of 1000
%   resamplings from the empirical sample (SAMPLE_DIST).
%
%   BOOTSTRAP_SE = SE_OF_MEDIAN(SAMPLE_DIST,BNO) resamples the sample
%   distribution BNO times.
%
%   See also MEDIAN.

% Input argument check
error(nargchk(1,2,nargin))
if nargin < 2
    bno = 1000;     % number or resamplings
end

% Get rid of NaNs
sample_dist = sample_dist(~isnan(sample_dist));

% Bootstrap SE
n = numel(sample_dist);  % sample size
bootstrap_mn = nan(1,bno);
for k = 1:bno
    bootstrap_sample = sample_dist(randi(n,1,n));   % resample with replacement
    bootstrap_mn(k) = median(bootstrap_sample);
end
bootstrap_se = std(bootstrap_mn);