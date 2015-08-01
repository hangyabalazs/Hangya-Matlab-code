function epcmetamodel
%EPCMETAMODEL   Rate constant dynamics for kinetic models.
%   EPCMETAMODEL models temporal changes of place cell-to-low rate cell
%   rate constant (low rate cell-to-place cell rate constant is kept fixed.)
%   It calls EPCKINETICMODELL for each day and hour (two nested cycles).
%   EPCKINETICMODEL simulates one hour with given rate constants. Note,
%   that this way, the rate constant dynamics is sampled every hour, and it
%   is not updated throughout an hour (this will be changed in future
%   versoins).
%
%   EPCMETAMODEL simulates rate constant changes as follows. Each day
%   begins with an abrupt increase, followed by an exponential decay. The
%   dynamics is of the form exp(-bt-a). 'b' is the time constant of the
%   decay and 'a' controls the extent of the increase.
%
%   Initial rate constants and concentrations represent a steady state
%   system, where place cells constitute 60% of all cells, place
%   cell-to-low rate cell rate constant is 0.004 and low rate cell-to-place
%   cell rate constant is 0.006 (kept fixed).
%
%   See also EPCKINETICMODEL.

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

% Equilibrium rate constants
k_PL0 = 0.004;

% Generate rate constant dynamics
x = 1:24*60;
ex = exp(-0.00135*x-4.8);
k_PL = k_PL0 + ex;
for dys = 2:8
    k_PL = [k_PL ex+k_PL(end)];
end
figure;plot(k_PL)

% First run
pc = 60;
[PC LC T] = epckineticmodel(pc,k_PL0);

% Loop days and hours (resolution will be minutes)
for dys = 1:8
    disp(dys)
    for hrs = 1:24
        pc = PC(end);   % start from last concentration
        [PCt LCt Tt] = epckineticmodel(pc,k_PL((dys-1)*24*60+(hrs-1)*60+1));
        PC = [PC PCt];   %#ok<AGROW>
        LC = [LC LCt];   %#ok<AGROW>
        T = [T Tt+(dys-1)*24*60+hrs*60];   %#ok<AGROW>
    end
    close all
end

% Plot
figure
hold on
P1 = plot(T,LC,'b');
P2 = plot(T,PC,'r');
legend([P2 P1],{'place cell' 'low rate cell'})
keyboard