function [PC LC T] = epcmetamodel5(constB,constA,constM,k_LP0)
%EPCMETAMODEL5   Rate constant dynamics for kinetic models (Model4A).
%   EPCMETAMODEL5 models temporal changes of low rate cell-to-place cell
%   rate constant (place cell-to-low rate cell rate constant is kept fixed.)
%   It calls EPCKINETICMODEL2 for each day and hour (two nested cycles).
%   EPCKINETICMODEL2 simulates one hour with given rate constants.
%
%   EPCMETAMODEL5 simulates k_LP rate constant changes as follows. Each day
%   begins with an abrupt increase, followed by an exponential decay. The
%   dynamics is of the form exp(-bt-a). 'b' is the time constant of the
%   decay and 'a' controls the extent of the increase (similar to
%   EPCMETAMODEL). On the other hand, an abrupt decrease of place cell
%   percentage is implemented concurrently with the abrupt rate constant
%   change (multiplication by constant 'M', similar to EPCMETAMODEL2).
%
%   [PC LC T] = EPCMETAMODEL5(B,A,M) runs the model with the constants
%   given as  input arguments. It returns the dynamics of place cells (PC)
%   and low rate cells (LC) with their corresponding time vector (T).
%
%   [PC LC T] = EPCMETAMODEL5(B,A,M,K_LP0) takes the initial LC to PC rate
%   constant as fourth input argument.
%
%   Initial rate constants and concentrations represent a steady state
%   system, where place cells constitute 70.75% of all cells. Initial LC to
%   PC rate constant determines PC to LC rate constant based on the initial
%   steady state (kept constant afterwards). It also determines the value
%   the LC to PC rate constant converges to after every abrupt increase
%   (since convergence is to steady state).
%
%   See also EPCKINETICMODEL2, ,EPCMETAMODEL4, EPCFIT5 and EPCERROR5.

% Input argument check
error(nargchk(0,4,nargin))
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

% Equilibrium rate constants
k_PL0 = (100 - 70.75) / 70.75 * k_LP0;   % 70.75% PC based on real data; steady state

% Generate rate constant dynamics
x = 1:24*60;
ex = exp(-constB*x-constA);
k_LP = k_LP0 + ex;
for dys = 2:8
    k_LP = [k_LP ex+k_LP(end)];
end
k_LP = [k_LP k_LP(end) k_LP(end)];   % we need two more elements for the loop
figure;
plot(k_LP)

% First run
pc = 70.75;
[PC LC T] = epckineticmodel2(pc,k_PL0*ones(1,61),k_LP0*ones(1,61),0);

% Loop days and hours (resolution will be minutes)
for dys = 1:8
%     disp(dys)
    for hrs = 1:24
        if isequal(hrs,1)
            pc = PC(end) * constM;   % abrupt decrease
        else
            pc = PC(end);   % start from last concentration
        end
        t0 = (dys - 1) * 24 * 60 + (hrs - 1) * 60 + 1;
        t1 = t0 + 61;
        [PCt LCt Tt] = epckineticmodel2(pc,k_PL0*ones(1,61),k_LP(t0:t1),0);
        PC = [PC PCt];   %#ok<AGROW>
        LC = [LC LCt];   %#ok<AGROW>
        T = [T Tt+(dys-1)*24*60+hrs*60];   %#ok<AGROW>
    end
%     close all
end

% Get the series of check points
cps = 60 + [0 30 3*60 3*24*60 3*24*60+30 3*24*60+3*60 7*24*60 7*24*60+30 7*24*60+3*60];
pcs = arrayfun(@(s)PC(find(T<s,1,'last')),cps);
lcs = arrayfun(@(s)LC(find(T<s,1,'last')),cps);

% Real data
realPC = [70.75 53.85 59.80 51.79 39.81 41.67 27.66 13.48 30.12];

% Least square error
R = sum((pcs-realPC).^2);
disp(R)

% Plot
close all
figure
hold on
P1 = plot(T,LC,'b');
P2 = plot(T,PC,'r');
plot(cps,realPC,'ko')
legend([P2 P1],{'place cell' 'low rate cell'})
% keyboard

% Supress output
if nargout == 0
    clear
end