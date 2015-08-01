function ucmmsim
%UCMMSIM   Simulation of stochastic Michaelis-Menten reaction.

% Michaelis-Menten reaction: E + S <-> C -> E + P

% Simulation of a continuous time discrete state Markov process (equivalent to
% M/M/1 problem in queuing theory).

% Setting the initial number of molecules (assuming deterministic initial state)
e = 100;
s = 100;
c = 0;
p = 0;
if e + c == 0
    error('No enzyme in the system.')
elseif s + c == 0
    error('No substrate in the system.')
end

% Setting rate constants
k_SC = 1;
k_CS = 1;
k_CP = 1;

% Initializing time
t = 0;
st = 100;   % stopping time

% Simulation
E = e;
S = s;
C = c;
P = p;
T = t;
while t < st
    if s + c == 0
        disp('All substrate has been converted to product.')
        break
    end
    [e,s,c,p,t] = onestep(e,s,c,p,t,k_SC,k_CS,k_CP);    % calculate one simulation step
    if e < 0 | s < 0 | c < 0 | p < 0    % programing error control
        error('Simulation error: negative quantity.')
    end
    E = [E e];
    S = [S s];
    C = [C c];
    P = [P p];
    T = [T t];
end

% Plot simulation results
figure
hold on
P1 = plot(T,E,'b');
P2 = plot(T,S,'r');
P3 = plot(T,C,'y');
P4 = plot(T,P,'k');
legend([P1 P2 P3 P4],{'enzyme' 'substrate' 'complex' 'product'})

% Uncertainty coefficients
[hx,hy,jh,Hx,Hy,Hxy,Hxcy,Hycx,Ixy,Uxy,Uyx] = b_entropy(C(50:end),P(50:end));    % Uxy: P->E; Uyx: E->P
disp(['P->C ' num2str(Uxy)])
disp(['C->P ' num2str(Uyx)])
1;

% -------------------------------------------------------------------------
function [e,s,c,p,t] = onestep(e,s,c,p,t,k_SC,k_CS,k_CP)

% Generating random reaction and waiting time (exponential)
% rand('twister', sum(100*fliplr(clock)));    % initialize the state of the random generator
if ~isequal(s,0) && ~isequal(c,0) && ~isequal(e,0)
    q_SC = k_SC * s * e;    % infinitesimal transmission probabilities
    q_CS = k_CS * c;
    q_CP = k_CP * c;
    Q = q_SC + q_CS + q_CP;
    r = rand(1) * Q + 0.5;   % uniform random number from U(0.5,Q+0.5)
    rr = round(r);      % discrete uniform random number on {1,2,...,Q}
    if rr <= q_SC       % draw reaction with a probability proportional to its reaction velocity
        R = 'SC';
    elseif rr <= q_SC + q_CS
        R = 'CS';
    else
        R = 'CP';
    end
%     rand('twister', sum(100*fliplr(clock)));
    dt = exprnd(1/Q);    % in Matlab, expectation value should be given as input argument, which 
                         % is the reciprocal of the parameter of the exponential distribution
elseif (isequal(s,0) | isequal(e,0)) && ~isequal(c,0)   % lack of enzyme or substrate
    q_CS = k_CS * c;
    q_CP = k_CP * c;
    Q = q_CS + q_CP;
    r = rand(1) * Q + 0.5;      % draw from C->E+S and C->E+P reactions
    rr = round(r);
    if rr <= q_CS
        R = 'CS';
    else
        R = 'CP';
    end
%     rand('twister', sum(100*fliplr(clock)));
    dt = exprnd(1/Q);
elseif isequal(c,0) && ~isequal(e,0) && ~isequal(s,0)   % lack of complex
    R = 'SC';      % E+S->C
%     rand('twister', sum(100*fliplr(clock)));
    q_SC = k_SC * s * e;
    dt = exprnd(1/q_SC);
end

% Modify the quantities due to the selected reaction step
switch R
    case 'SC'
        e = e - 1;
        s = s - 1;
        c = c + 1;
    case 'CS'
        e = e + 1;
        s = s + 1;
        c = c - 1;
    case 'CP'
        e = e + 1;
        c = c - 1;
        p = p + 1;
end
t = t + dt;     % increase time