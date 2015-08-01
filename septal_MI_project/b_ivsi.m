function b_ivsi
%IVSI   Plots "inter-first-spike interval" length versus eeg theta segment minimum location distances.
%
%   See also IVSIRUN.

% Input arguments check
error(nargchk(0,0,nargin));

% Close all
close all;

% Import
b_in3;

% Ivsi
[burstloc,intervals,minimums] = b_ivsi_main;