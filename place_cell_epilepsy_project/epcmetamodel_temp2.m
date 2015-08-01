function [PC LC T] = epcmetamodel(constB,constA,k_PL0)
%EPCMETAMODEL   Rate constant dynamics for kinetic models.
%   EPCMETAMODEL models temporal changes of place cell-to-low rate cell
%   rate constant (low rate cell-to-place cell rate constant is kept fixed.)
%   It calls EPCKINETICMODELL for each day and hour (two nested cycles).
%   EPCKINETICMODEL simulates one hour with given rate constants.
%
%   EPCMETAMODEL simulates rate constant changes as follows. Each day
%   begins with an abrupt increase, followed by an exponential decay. The
%   dynamics is of the form exp(-bt-a). 'b' is the time constant of the
%   decay and 'a' controls the extent of the increase.
%
%   [PC LC T] = EPCMETAMODEL(B,A) runs the model with the constants given 
%   as  input arguments. It returns the dynamics of place cells (PC) and
%   low rate cells (LC) with their corresponding time vector (T).
%
%   [PC LC T] = EPCMETAMODEL(B,A,K_PL0) takes the initial PC to LC rate
%   constant as third input argument.
%
%   Initial rate constants and concentrations represent a steady state
%   system, where place cells constitute 70.75% of all cells. Initial PC to
%   LC rate constant determines LC to PC rate constant based on the initial
%   steady state (kept constant afterwards). It also determines the value
%   the PC to LC rate constant converges to after every abrupt increase
%   (since convergence is to steady state).
%
%   See also EPCKINETICMODEL, EPCFIT and EPCERROR.

% calls epckineticmodel
% epckineticmodell simulates one hour with given parameters
% parameters change in epcmetamodel
% cycle for the hours and days
% rate constant changes: get abrubtly multiplied by sg, and slowly decays (two params)
% these two params we should fit to real data
% params of epckineticmodell: initial concentrations - should estimate from control data
% and steady-state rate constants: sets the noise - can be set arbitrarily?

% fitting to data: a third program should get the samples corresponding to
% the recordings
% least-square fit for the two parameters? might be VERY SLOW...

% Input argument check
error(nargchk(0,3,nargin))
if nargin < 3
    k_PL0 = 0.0001;
end
if nargin < 2
    constA = 2.5;
end
if nargin < 1
    constB = 0.1;
end

% Equilibrium rate constants
k_LP0 = 70.75 / (100 - 70.75) * k_PL0;   % 70.75% PC based on real data; steady state

% Generate rate constant dynamics
x = 1:24*60;
ex = exp(-constB*x-constA);
k_PL = k_PL0 + ex;
for dys = 2:8
    k_PL = [k_PL ex+k_PL(end)];
end
k_PL = [k_PL k_PL(end) k_PL(end)];   % we need two more elements for the loop
figure;plot(k_PL)

% First run (3 hours)
pc = 74.06;
[PC LC T] = epckineticmodel(pc,k_PL0*ones(1,61),k_LP0,0);
for phrs = 2:3
    [PCt LCt Tt] = epckineticmodel(pc,k_PL0*ones(1,61),k_LP0,0);
    PC = [PC PCt];   %#ok<AGROW>
    LC = [LC LCt];   %#ok<AGROW>
    T = [T Tt+phrs*60];   %#ok<AGROW>
end

% Loop days and hours (resolution will be minutes)
for dys = 1:8
%     disp(dys)
    for hrs = 1:24
        pc = PC(end);   % start from last concentration
        t0 = (dys - 1) * 24 * 60 + (hrs - 1) * 60 + 1;
        t1 = t0 + 61;
        [PCt LCt Tt] = epckineticmodel(pc,k_PL(t0:t1),k_LP0,0);
        PC = [PC PCt];   %#ok<AGROW>
        LC = [LC LCt];   %#ok<AGROW>
        T = [T Tt+(dys-1)*24*60+hrs*60];   %#ok<AGROW>
    end
%     close all
end

% Plot
close all
figure
hold on
P1 = plot(T,LC,'b');
P2 = plot(T,PC,'r');
legend([P2 P1],{'place cell' 'low rate cell'})
% keyboard

% Supress output
if nargout == 0
    clear
end