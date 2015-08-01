function [PC LC T] = epckineticmodel(pc,k_PL,k_LP,isplot)
%EPCKINETICMODEL   Simulation of stochastic 'kinetic reaction'.
%   EPCKINETICMODEL simulates a simple reversible reaction with two rate
%   constants.
%   Reaction: PC <-> LC
%   It is a simulation of a continuous time discrete state Markov process
%   (equivalent to M/M/1 problem in queuing theory).
%
%   EPCKINETICMODEL(PC,K_PL,K_LP) receives the initial percentage of PC as
%   well as the PC to LC and LC to PC rate constants as input arguments. In
%   this implementation K_PL can change in time (it should be a vector),
%   whereas K_LP is kept constant throughout a run.
%
%   [PC LC T] = EPCKINETICMODEL(PC,K_PL,K_LP) returns the concentrations
%   (PC, LC) along with their time vector (T).
%
%   [PC LC T] = EPCKINETICMODEL(PC,K_PL,K_LP,ISPLOT) takes a fourth input
%   argument. If ISPLOT = 1, the dynamics are displayed (default).
%
%   See also EPCMETAMODEL.

% Input argument check
error(nargchk(0,4,nargin))
if nargin < 4
    isplot = 1;
end

% Setting the initial number of cells (assuming deterministic initial state)
if nargin < 1
    pc = 70.75;   % place cells (%)
end
lc = 100 - pc;   % low rate cells (%)
N = 100;
pc = round(pc*N);   % convert to sufficient quantites from percentages
lc = round(lc*N);

% Setting rate constants
if nargin < 3
    k_LP = 0.0001;
end
if nargin < 2
    pk_PL = lc / pc * k_LP;    % default: steady state
    k_PL = ones(1,60) * pk_PL;
end

% Initializing time
t = 0;
st = 60;   % stopping time

% Simulation
PC = pc;
LC = lc;
T = eps * 1000;    % numerically, k < k + eps * 1000
while t < st
%     disp(t)
    if pc == 0
        disp('All place cells has been converted to low rate cells.')
        break
    end
    [pc,lc,t] = onestep(pc,lc,t,k_PL(floor(t)+1),k_LP);    % calculate one simulation step
    if pc < 0 || lc < 0    % programing error control
        error('Simulation error: negative quantity.')
    end
    PC = [PC pc]; %#ok<AGROW>
    LC = [LC lc]; %#ok<AGROW>
    T = [T t]; %#ok<AGROW>
end
if length(PC) > 1   % keep initial values, if ther is nothing else
    PC = PC(1:end-1);   % last time instance > st
    LC = LC(1:end-1);
    T = T(1:end-1);
end

% Convert numbers back to percentages
PC = PC / N;
LC = LC / N;

% Plot simulation results
if isplot
    figure
    hold on
    P1 = plot(T,LC,'b');
    P2 = plot(T,PC,'r');
    legend([P2 P1],{'place cell' 'low rate cell'})
end



% -------------------------------------------------------------------------
function [pc,lc,t] = onestep(pc,lc,t,k_PL,k_LP)

% Generating random reaction and waiting time (exponential)
% rand('twister', sum(100*fliplr(clock)));    % initialize the state of the random generator
if ~isequal(pc,0) && ~isequal(lc,0)
    q_PL = k_PL * pc;    % infinitesimal transmission probabilities
    q_LP = k_LP * lc;
    Q = q_PL + q_LP;
    cst = 10 ^ 10;
    r = rand(1) * Q * cst + 0.5;   % uniform random number from U(0.5,Q*cst+0.5); cst multiplyer is to get rid of rounding problems
    rr = round(r);      % discrete uniform random number on {1,2,...,Q*cst}
    if rr <= q_PL * cst      % draw reaction with a probability proportional to its reaction velocity
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

if pc < 0 || lc < 0    % programing error control
    error('Simulation error: negative quantity.')
end