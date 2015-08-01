function [constB constA constM k_LP0 f] = epcfit6(constB,constA,constM,k_LP0)
%EPCFIT6   Fits kinetic model on real data.
%   EPCFIT6 tries to fit the two-compound kinetic model realized by
%   EPCMETAMODEL6 and EPCKINETICMODEL2 on real data, using a linear
%   optimization with least square error.
%
%   [CONSTB CONSTA CONSTM K_LP0 ERR] = EPCFIT6 returns the optimal model
%   parameters as well as the corresponding fitting error.
%
%   [CONSTB CONSTA CONSTM K_LP0 F] = EPCFIT6(CONSTB,CONSTA,CONSTM,K_LP0)
%   uses the input parameters as initial state.
%
%   See also EPCKINETICMODEL2, EPCMETAMODEL6, EPCERROR6, EPCFITCALL6 and
%   FMINSEARCHBND.

% Initial parameters
if nargin < 4
    k_LP0 = 0.00001;
end
if nargin < 3
    constM = 15;
end
if nargin < 2
    constA = 7;
end
if nargin < 1
    constB = 0.006;
end

% Bounds
LB = [0 0.5 0 0];   % lower bound
UB = [0.9 30 100 0.1];   % upper bound

% Linear programming
global COUNTER
COUNTER = 0;    % count the iterations
[x,f,exitflag,output] = fminsearchbnd('epcerror6',[constB constA constM k_LP0],LB,UB,optimset('MaxFunEvals',150));

% Parameters
constB = x(1);
constA = x(2);
constM = x(3);
k_LP0 = x(4);

% Plot
global BESTDYNAMICS
figure
hold on
P1 = plot(BESTDYNAMICS.T,BESTDYNAMICS.LC,'b');
P2 = plot(BESTDYNAMICS.T,BESTDYNAMICS.PC,'r');
plot(BESTDYNAMICS.REALTIME,BESTDYNAMICS.REALDATA,'ko')
legend([P2 P1],{'place cell' 'low rate cell'})