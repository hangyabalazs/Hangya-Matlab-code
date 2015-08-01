function [constB constA constM k_LP0 f] = epcfit5(constB,constA,constM,k_LP0)
%EPCFIT5   Fits kinetic model on real data.
%   EPCFIT5 tries to fit the two-compound kinetic model realized by
%   EPCMETAMODEL5 and EPCKINETICMODEL2 on real data, using a linear
%   optimization with least square error.
%
%   [CONSTB CONSTA CONSTM K_LP0 ERR] = EPCFIT5 returns the optimal model
%   parameters as well as the corresponding fitting error.
%
%   [CONSTB CONSTA CONSTM K_LP0 F] = EPCFIT5(CONSTB,CONSTA,CONSTM,K_LP0)
%   uses the input parameters as initial state.
%
%   See also EPCKINETICMODEL2, EPCMETAMODEL5, EPCERROR5 and FMINSEARCHBND.

% Initial parameters
if nargin < 4
    k_LP0 = 0.00001;
end
if nargin < 3
    constM = 0.65;
end
if nargin < 2
    constA = 7;
end
if nargin < 1
    constB = 0.006;
end

% Bounds
LB = [0 0.5 0 0];   % lower bound
UB = [0.9 30 1 0.1];   % upper bound

% Linear programming
global COUNTER
COUNTER = 0;    % count the iterations
[x,f,exitflag,output] = fminsearchbnd('epcerror5',[constB constA constM k_LP0],LB,UB,optimset('MaxFunEvals',150));

% Parameters
constB = x(1);
constA = x(2);
constM = x(3);
k_LP0 = x(4);

% Plot
global DYNAMICS
figure
hold on
P1 = plot(DYNAMICS.T,DYNAMICS.LC,'b');
P2 = plot(DYNAMICS.T,DYNAMICS.PC,'r');
plot(DYNAMICS.REALTIME,DYNAMICS.REALDATA,'ko')
legend([P2 P1],{'place cell' 'low rate cell'})