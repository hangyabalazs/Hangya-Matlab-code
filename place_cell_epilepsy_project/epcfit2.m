function [constM k_PL0 f] = epcfit2(constM,k_PL0)
%EPCFIT2   Fits kinetic model on real data.
%   EPCFIT2 tries to fit the two-compound kinetic model realized by
%   EPCMETAMODEL2 and EPCKINETICMODEL on real data, using a linear
%   optimization with least square error.
%
%   [CONSTM K_PL0 ERR] = EPCFIT2 returns the optimal model parameters as
%   well as the corresponding fitting error.
%
%   [CONSTM K_PL0 F] = EPCFIT2(CONSTM,K_PL0) uses the input parameters as
%   initial state.
%
%   See also EPCKINETICMODEL, EPCMETAMODEL2, EPCERROR2, EPCFITCALL2 and
%   FMINSEARCHBND.

% Initial parameters
if nargin < 2
    k_PL0 = 0.00002;
end
if nargin < 1
    constM = 0.7;
end

% Bounds
LB = [0 0];   % lower bound
UB = [1 0.01];   % upper bound

% Linear programming
global COUNTER
COUNTER = 0;    % count the iterations
[x,f,exitflag,output] = fminsearchbnd('epcerror2',[constM k_PL0],LB,UB,optimset('MaxFunEvals',150));

% Parameters
constM = x(1);
k_PL0 = x(2);

% Plot
global DYNAMICS
figure
hold on
P1 = plot(DYNAMICS.T,DYNAMICS.LC,'b');
P2 = plot(DYNAMICS.T,DYNAMICS.PC,'r');
plot(DYNAMICS.REALTIME,DYNAMICS.REALDATA,'ko')
legend([P2 P1],{'place cell' 'low rate cell'})