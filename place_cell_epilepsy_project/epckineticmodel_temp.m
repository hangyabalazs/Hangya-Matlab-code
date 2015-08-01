function [PC LC T] = epckineticmodel(pc,k_PL)
%EPCKINETICMODEL   Simulation of stochastic 'kinetic reaction'.
%   EPCKINETICMODEL simulates a simple reversible reaction with two rate
%   constants.
%   Reaction: PC <-> LC
%   It is a simulation of a continuous time discrete state Markov process
%   (equivalent to M/M/1 problem in queuing theory).
%
%   EPCKINETICMODEL(PC,K_PL) receives the initial percentage of PC and the
%   PC to LC rate constants as input arguments.
%
%   [PC LC T] = epckineticmodel(pc,k_PL) returns the concentrations (PC,
%   LC) along with their time vector (T).
%
%   See also EPCMETAMODEL.

% Setting the initial number of cells (assuming deterministic initial state)
if nargin < 1
    pc = 60;   % place cells (%)
end
lc = 100 - pc;   % low rate cells (%)
N = 100;
pc = pc * N;   % convert to sufficient quantites from percentages
lc = lc * N;

% Setting rate constants
k_LP = 0.006;
if nargin < 2
    k_PL = 0.004;
end

% Initializing time
t = 0;
st = 60;   % stopping time

% Simulation
PC = pc;
LC = lc;
T = t;
while t < st
%     disp(t)
    if pc == 0
        disp('All place cells has been converted to low rate cells.')
        break
    end
    [pc,lc,t] = onestep(pc,lc,t,k_PL,k_LP);    % calculate one simulation step
    if pc < 0 || lc < 0    % programing error control
        error('Simulation error: negative quantity.')
    end
    PC = [PC pc]; %#ok<AGROW>
    LC = [LC lc]; %#ok<AGROW>
    T = [T t]; %#ok<AGROW>
end

% Convert numbers back to percentages
PC = PC / N;
LC = LC / N;

% Plot simulation results
figure
hold on
P1 = plot(T,LC,'b');
P2 = plot(T,PC,'r');
legend([P2 P1],{'place cell' 'low rate cell'})



% -------------------------------------------------------------------------
function [pc,lc,t] = onestep(pc,lc,t,k_PL,k_LP)

% Generating random reaction and waiting time (exponential)
% rand('twister', sum(100*fliplr(clock)));    % initialize the state of the random generator
if ~isequal(pc,0) && ~isequal(lc,0)
    q_PL = k_PL * pc;    % infinitesimal transmission probabilities
    q_LP = k_LP * lc;
    Q = q_PL + q_LP;
    r = rand(1) * Q + 0.5;   % uniform random number from U(0.5,Q+0.5)
    rr = round(r);      % discrete uniform random number on {1,2,...,Q}
    if rr <= q_PL       % draw reaction with a probability proportional to its reaction velocity
        R = 'PL';
    else
        R = 'LP';
    end
%     rand('twister', sum(100*fliplr(clock)));
    dt = exprnd(1/Q);    % in Matlab, expectation value should be given as input argument, which 
                         % is the reciprocal of the parameter of the exponential distribution
elseif isequal(pc,0) && ~isequal(lc,0)
    R = 'LP';
    q_LP = k_LP * lc;
    dt = exprnd(1/q_LP);
elseif isequal(lc,0) && ~isequal(pc,0)
    R = 'PL';
    q_PL = k_PL * pc;
    dt = exprnd(1/q_PL);
end

% Modify the quantities due to the selected reaction step
switch R
    case 'PL'
        pc = pc - 1;
        lc = lc + 1;
    case 'LP'
        pc = pc + 1;
        lc = lc - 1;
end
t = t + dt;     % increase time