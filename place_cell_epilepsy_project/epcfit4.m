function [constB constA constM k_PL0 f] = epcfit4(constB,constA,constM,k_PL0)
%EPCFIT4   Fits kinetic model on real data.
%   EPCFIT4 tries to fit the two-compound kinetic model realized by
%   EPCMETAMODEL4 and EPCKINETICMODEL on real data, using a linear
%   optimization with least square error.
%
%   [CONSTB CONSTA CONSTM K_PL0 ERR] = EPCFIT4 returns the optimal model
%   parameters as well as the corresponding fitting error.
%
%   [CONSTB CONSTA CONSTM K_PL0 F] = EPCFIT4(CONSTB,CONSTA,CONSTM,K_PL0)
%   uses the input parameters as initial state.
%
%   See also EPCKINETICMODEL, EPCMETAMODEL4, EPCERROR4, EPCFITCALL4 and
%   FMINSEARCHBND.

% Initial parameters
if nargin < 4
    k_PL0 = 0.0002;
end
if nargin < 3
    constM = 0.7;
end
if nargin < 2
    constA = 1;
end
if nargin < 1
    constB = 0.05;
end

% Bounds
LB = [0 0.5 0 0];   % lower bound
UB = [0.9 30 1 0.01];   % upper bound

% Linear programming
global COUNTER
COUNTER = 0;    % count the iterations
[x,f,exitflag,output] = fminsearchbnd('epcerror4',[constB constA constM k_PL0],LB,UB,optimset('MaxFunEvals',150));

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