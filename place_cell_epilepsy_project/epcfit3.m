function [constB constA k_PL0 f] = epcfit3(constB,constA,k_PL0)
%EPCFIT3   Fits kinetic model on real data.
%   EPCFIT3 tries to fit the two-compound kinetic model realized by
%   EPCMETAMODEL3 and EPCKINETICMODEL on real data, using a linear
%   optimization with least square error.
%
%   [CONSTB CONSTA K_PL0 ERR] = EPCFIT3 returns the optimal model parameters
%   as well as the corresponding fitting error.
%
%   [CONSTB CONSTA K_PL0 F] = EPCFIT3(CONSTB,CONSTA,K_PL0) uses the input
%   parameters as initial state.
%
%   See also EPCKINETICMODEL, EPCMETAMODEL3, EPCERROR3, EPCFITCALL3 and
%   FMINSEARCHBND.

% Initial parameters
if nargin < 3
    k_PL0 = 0.00002;
end
if nargin < 2
    constA = 25;
end
if nargin < 1
    constB = 0.05;
end

% Bounds
LB = [0 0 0];   % lower bound
UB = [100 100 0.1];   % upper bound

% Linear programming
global COUNTER
COUNTER = 0;    % count the iterations
[x,f,exitflag,output] = fminsearchbnd('epcerror3',[constB constA k_PL0],LB,UB,optimset('MaxFunEvals',150));

% Parameters
constB = x(1);
constA = x(2);
k_PL0 = x(3);

% Plot
global DYNAMICS
figure
hold on
P1 = plot(DYNAMICS.T,DYNAMICS.LC,'b');
P2 = plot(DYNAMICS.T,DYNAMICS.PC,'r');
plot(DYNAMICS.REALTIME,DYNAMICS.REALDATA,'ko')
legend([P2 P1],{'place cell' 'low rate cell'})