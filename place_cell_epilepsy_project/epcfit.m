function [constB constA k_PL0 f] = epcfit(constB,constA,k_PL0)
%EPCFIT   Fits kinetic model on real data.
%   EPCFIT tries to fit the two-compound kinetic model realized by
%   EPCMETAMODEL and EPCKINETICMODEL on real data, using a linear
%   optimization with least square error.
%
%   [CONSTB CONSTA K_PL0 ERR] = EPCFIT returns the optimal model parameters
%   as well as the corresponding fitting error.
%
%   [CONSTB CONSTA K_PL0 F] = EPCFIT(CONSTB,CONSTA,K_PL0) uses the input
%   parameters as initial state.
%
%   See also EPCKINETICMODEL, EPCMETAMODEL, EPCERROR, EPCFITCALL and
%   FMINSEARCHBND.

% Initial parameters
if nargin < 3
    k_PL0 = 0.0001;
end
if nargin < 2
    constA = 4.5;
end
if nargin < 1
    constB = 0.002;
end

% Bounds
LB = [0 0.5 0.0001];   % lower bound
UB = [0.9 30 0.007];   % upper bound

% Linear programming
global COUNTER
COUNTER = 0;    % count the iterations
[x,f,exitflag,output] = fminsearchbnd('epcerror',[constB constA k_PL0],LB,UB,optimset('MaxFunEvals',150));

% Parameters
constB = x(1);
constA = x(2);
k_PL0 = x(3);

% Plot
global BESTDYNAMICS
figure
hold on
P1 = plot(BESTDYNAMICS.T,BESTDYNAMICS.LC,'b','LineWidth',2);
P2 = plot(BESTDYNAMICS.T,BESTDYNAMICS.PC,'r','LineWidth',2);
plot(BESTDYNAMICS.REALTIME,BESTDYNAMICS.REALDATA,'ko')
set(gca,'LineWidth',2,'XTick',[])
legend([P2 P1],{'place cell' 'low rate cell'})