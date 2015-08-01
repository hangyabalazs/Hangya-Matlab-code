function D = JSdiv_discrete(X,Y)
%JSDIV_DISCRETE   Jensen-Shannon divergence for discrete variables.
%   D = JSDIV_DISCRETE(X,Y) calculates the Jensen-Shannon divergence
%   version) of the two discrete variables. Note that sqrt(2*D) defines a
%   metric.
%
%   See also JSDIV.

% Input argument check
error(nargchk(2,2,nargin))

% Sort to corresponding bins
bins = unique([X; Y]);
lb = length(bins);
P_X = arrayfun(@(s)sum(X==bins(s)),1:lb)' / length(X);       % X distribution
P_Y = arrayfun(@(s)sum(Y==bins(s)),1:lb)' / length(Y);       % Y distribution

% Kullback-Leibler divergence
D = JSdiv(P_X,P_Y);