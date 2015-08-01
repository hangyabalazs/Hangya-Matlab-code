function [constM k_PL0 f] = epcfit2b(constM,k_PL0)
%EPCFIT2B   Fits kinetic model on real data.
%   EPCFIT2B tries to fit the two-compound kinetic model realized by
%   EPCMETAMODEL2B and EPCKINETICMODEL on real data, using a linear
%   optimization with least square error.
%
%   [CONSTM K_PL0 ERR] = EPCFIT2B returns the optimal model parameters as
%   well as the corresponding fitting error.
%
%   [CONSTM K_PL0 F] = EPCFIT2B(CONSTM,K_PL0) uses the input parameters as
%   initial state.
%
%   See also EPCKINETICMODEL, EPCMETAMODEL2B, EPCERROR2B, EPCFITCALL2B and
%   FMINSEARCHBND.

% Initial parameters
if nargin < 2
    k_PL0 = 0.00002;
end
if nargin < 1
    constM = 15;
end

% Bounds
LB = [0 0];   % lower bound
UB = [100 0.01];   % upper bound

% Linear programming
global COUNTER
COUNTER = 0;    % count the iterations
[x,f,exitflag,output] = fminsearchbnd('epcerror2b',[constM k_PL0],LB,UB,optimset('MaxFunEvals',150));

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