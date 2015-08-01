function theta = b_rrmixedvm(n)
%RRMIXEDVM   Random values from mixed random parameter von Mises distributions.
%   [THETA] = RRMIXEDVM(N) generates uniform random parameters for two von
%   Mises distribution and mixes them using uniform random mixing ratio. It
%   generates N random values from this mixture.
%
%   See also RVM and RMIXEDVM.

% Input argumnet check
error(nargchk(1,1,nargin))

% Uniform random parameter generation
m1 = rand(1) * 2 * pi;
m2 = rand(1) * 2 * pi;
k1 = rand(1) * 9.5 + 0.5;
k2 = rand(1) * 9.5 + 0.5;
pp = rand(1);

% Random sample from mixed von Mises
theta = b_rmixedvm(n,m1,m2,k1,k2,pp);