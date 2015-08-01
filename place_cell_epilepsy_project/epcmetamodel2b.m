function [PC LC T] = epcmetamodel2b(constM,k_PL0)
%EPCMETAMODEL2B   Rate constant dynamics for kinetic models (Model1B).
%   EPCMETAMODEL2B models temporal changes of place cell percentage. Both
%   rate constants (low rate cell-to-place cell and place cell-to-low rate
%   cell) are kept fixed. It calls EPCKINETICMODELL for each day and hour
%   (two nested cycles). EPCKINETICMODEL simulates one hour with given rate
%   constants and place cell percentage.
%
%   EPCMETAMODEL2B simulates PC percentage changes as follows. Each day
%   begins with an abrupt decrease. After that, percentages are only
%   changed within the kinetic model. The decrease is controled by a
%   constant ('M'), that is subtracted from the current PC percentage to
%   provide the new value.
%
%   [PC LC T] = EPCMETAMODEL2B(M) runs the model with the constant M given 
%   as  input argument. It returns the dynamics of place cells (PC) and
%   low rate cells (LC) with their corresponding time vector (T).
%
%   [PC LC T] = EPCMETAMODEL2B(M,K_PL0) takes the initial PC to LC rate
%   constant as second input argument.
%
%   Initial rate constants and concentrations represent a steady state
%   system, where place cells constitute 70.75% of all cells. Initial PC to
%   LC rate constant determines LC to PC rate constant based on the initial
%   steady state (kept constant afterwards).
%
%   See also EPCKINETICMODEL, EPCFIT2 and EPCERROR2.

% Input argument check
error(nargchk(0,2,nargin))
if nargin < 2
    k_PL0 = 0.0002;
end
if nargin < 1
    constM = 15;
end

% Equilibrium rate constants
k_LP0 = 70.75 / (100 - 70.75) * k_PL0;   % 70.75% PC based on real data; steady state

% Generate rate constant dynamics
k_PL = ones(1,8*24*60+2) * k_PL0;

% First run
pc = 70.75;
[PC LC T] = epckineticmodel(pc,k_PL0*ones(1,61),k_LP0,0);

% Loop days and hours (resolution will be minutes)
for dys = 1:8
%     disp(dys)
    for hrs = 1:24
        if isequal(hrs,1)
            pc = PC(end) - constM;   % abrupt decrease
            pc = max(pc,0);
        else
            pc = PC(end);   % start from last concentration
        end
        t0 = (dys - 1) * 24 * 60 + (hrs - 1) * 60 + 1;
        t1 = t0 + 61;
        [PCt LCt Tt] = epckineticmodel(pc,k_PL(t0:t1),k_LP0,0);
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