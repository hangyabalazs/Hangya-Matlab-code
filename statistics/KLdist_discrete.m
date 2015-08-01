function D = KLdist_discrete(X,Y)
%KLDIST_DISCRETE   Kullbach-Leibler distance for discrete variables.
%   D = KLDIST_DISCRETE(X,Y) calculates the Kullbach-Leibler distance (information
%   divergence) of the two discrete variables.
%
%   See also KLDIST.

% Input argument check
error(nargchk(2,2,nargin))

% Sort to corresponding bins
bins = unique([X; Y]);
lb = length(bins);
P_X = arrayfun(@(s)sum(X==bins(s)),1:lb)' / length(X);       % X distribution
P_Y = arrayfun(@(s)sum(Y==bins(s)),1:lb)' / length(Y);       % Y distribution

% Kullback-Leibler divergence
D = KLdist(P_X,P_Y);